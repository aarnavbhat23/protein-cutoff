"""Which LiveProteinBench questions survive increasingly strict decontamination?

  python -m cutoffaudit.subset   # results/clean_subset.tsv (question ids + flags), results/clean_subset.json

Tiers (each implies the previous):
  t1_answer_not_in_own_entry   correct option not supported by the protein's
                               own last pre-cutoff UniProtKB entry
  t2_no_relative_answer        and, for the seven annotation tasks, no
                               Swiss-Prot 2024_06 relative (>= 30% identity)
                               supports any option
  t3_sequence_new              and the sequence was first seen on/after the cutoff
  t4_family_new                and its UniRef50 cluster has no member created
                               before the cutoff
"""

from __future__ import annotations

import csv
import json
from collections import Counter

from .audit import RES, read_tsv
from . import benchmark


def run() -> dict:
    items = read_tsv(RES / "items.tsv")
    prot = {r["accession"]: r for r in read_tsv(RES / "proteins.tsv")}
    hom = {r["id"]: r for r in read_tsv(RES / "homolog_items.tsv")}
    cut = benchmark.CUTOFF
    rows = []
    for i in items:
        p = prot[i["accession"]]
        h = hom.get(i["id"])
        t1 = i["own_status"] != "leaked" and i["own_status"] != "ambiguous"
        t2 = t1 and (h is None or h["answered"] == "False")
        t3 = t2 and p["first_seen"] >= cut
        t4 = t3 and int(p["cluster_members_before_cutoff"] or -1) == 0
        rows.append({"id": i["id"], "task": i["task"], "accession": i["accession"],
                     "t1_answer_not_in_own_entry": t1, "t2_no_relative_answer": t2, "t3_sequence_new": t3, "t4_family_new": t4})
    with open(RES / "clean_subset.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    tiers = ["t1_answer_not_in_own_entry", "t2_no_relative_answer", "t3_sequence_new", "t4_family_new"]
    out = {"questions": len(rows), "tiers": {t: sum(r[t] for r in rows) for t in tiers},
           "by_task": {task: {"questions": sum(r["task"] == task for r in rows),
                              **{t: sum(r[t] for r in rows if r["task"] == task) for t in tiers}}
                       for task in sorted({r["task"] for r in rows})}}
    (RES / "clean_subset.json").write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
