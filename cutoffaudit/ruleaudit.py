"""Audit a date-based selection *rule*, not just one benchmark's items.

Rebuilds the full candidate pool a rule implies -- every reviewed UniProtKB
entry with `date_created` on or after a cutoff -- and dates every sequence
(UniParc first seen) and a random sample of their UniRef50 clusters.

Rules audited (results/rules.json):
  liveproteinbench   reviewed, created >= 2025-01-01, 50-1022 aa   (arXiv:2512.22257)
  vibeproteinbench   reviewed, created >= 2025-09-01, 50-1022 aa   (arXiv:2605.10978,
                     "UniProt entries registered after September 1, 2025"; data not released)

  python -m cutoffaudit.ruleaudit
"""

from __future__ import annotations

import concurrent.futures as cf
import datetime as dt
import json
import random

from . import uniprot

RULES = {"liveproteinbench": "2025-01-01", "vibeproteinbench": "2025-09-01"}
CLUSTER_SAMPLE = 600


def audit(cutoff: str, seed: int = 0) -> dict:
    q = f"(reviewed:true) AND (date_created:[{cutoff} TO *]) AND (length:[50 TO 1022]) AND (fragment:false)"
    rows = uniprot.stream(q, "accession,date_created,length")
    accs = [r["Entry"] for r in rows]
    fs = {r["From"]: r["First seen"] for r in uniprot.idmap(accs, "UniParc", "first_seen")}
    u50 = {r["From"]: r["Cluster ID"] for r in uniprot.idmap(accs, "UniRef50", "id")}
    dated = [fs[a] for a in accs if a in fs]
    cut = dt.date.fromisoformat(cutoff)
    year_before = (cut - dt.timedelta(days=365)).isoformat()

    sample = random.Random(seed).sample(accs, min(CLUSTER_SAMPLE, len(accs)))
    before = (cut - dt.timedelta(days=1)).isoformat()

    def old_members(a):
        c = u50.get(a)
        return uniprot.count(f"uniref_cluster_50:{c} AND date_created:[* TO {before}]") if c else -1

    with cf.ThreadPoolExecutor(4) as ex:
        counts = [c for c in ex.map(old_members, sample) if c >= 0]
    new_seq = [a for a in accs if fs.get(a, "") >= cutoff]
    return {
        "cutoff": cutoff, "entries": len(accs), "dated": len(dated),
        "sequence_first_seen_before_cutoff": sum(f < cutoff for f in dated) / len(dated),
        "sequence_first_seen_more_than_1y_before": sum(f < year_before for f in dated) / len(dated),
        "sequence_first_seen_before_2021_04": sum(f < "2021-04-01" for f in dated) / len(dated),
        "median_first_seen": sorted(dated)[len(dated) // 2],
        "new_sequences": len(new_seq),
        "cluster_sample": len(counts),
        "cluster_has_member_before_cutoff": sum(c > 0 for c in counts) / len(counts),
    }


def run() -> dict:
    out = {name: audit(c) for name, c in RULES.items()}
    (uniprot.ROOT / "results" / "rules.json").write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
