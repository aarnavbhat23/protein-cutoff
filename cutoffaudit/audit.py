"""Was a 'post-cutoff' benchmark really unseen before its cutoff?

Three questions per protein / question:

  dates     When was the sequence first seen anywhere (UniParc)? When did the
            accession first appear in UniProtKB, and when did it reach
            Swiss-Prot? Does its UniRef50 cluster hold any member created
            before the cutoff?
  own       Using the protein's last UniProtKB entry version released before
            the cutoff (UniSave), does that entry already support the correct
            option, and only the correct one?
  current   The same check against today's entry. This calibrates the
            matchers: the correct option should nearly always be supported now.

  python -m cutoffaudit.audit dates
  python -m cutoffaudit.audit own
  python -m cutoffaudit.audit summary
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import datetime as dt
import json
from collections import Counter, defaultdict
from pathlib import Path

from . import benchmark, uniprot
from .flatfile import parse
from .match import classify

ROOT = uniprot.ROOT
RES = ROOT / "results"
CUT = dt.date.fromisoformat(benchmark.CUTOFF)


def pmap(fn, xs, workers: int = 4):
    with cf.ThreadPoolExecutor(workers) as ex:
        return list(ex.map(fn, xs))


# ----------------------------------------------------------------- dates
def protein_dates(accs: list[str], cut: dt.date = CUT) -> dict[str, dict]:
    up = {r["From"]: (r["Entry"], r["First seen"]) for r in uniprot.idmap(accs, "UniParc", "upi,first_seen")}
    u50 = {r["From"]: r["Cluster ID"] for r in uniprot.idmap(accs, "UniRef50", "id,count")}

    def history(acc):
        h = uniprot.unisave_history(acc)
        if not h:
            return acc, {}
        dates = sorted((uniprot.release_date(e), e) for e in h)
        sp = [d for d, e in dates if e["database"] == "Swiss-Prot"]
        seqv = sorted({e["sequenceVersion"] for e in h})
        last_seq_change = min((d for d, e in dates if e["sequenceVersion"] == seqv[-1]), default=None)
        return acc, {"uniprot_first": dates[0][0].isoformat(), "swissprot_first": sp[0].isoformat() if sp else "",
                     "sequence_versions": len(seqv), "current_sequence_since": last_seq_change.isoformat() if last_seq_change else ""}

    hist = dict(pmap(history, accs))

    def cluster_old(cid):
        n = uniprot.count(f"uniref_cluster_50:{cid} AND date_created:[* TO {(cut - dt.timedelta(days=1)).isoformat()}]")
        return cid, n

    clusters = dict(pmap(cluster_old, sorted({c for c in u50.values() if c})))
    out = {}
    for a in accs:
        upi, fs = up.get(a, ("", ""))
        cid = u50.get(a, "")
        out[a] = {"upi": upi, "first_seen": fs, "uniref50": cid,
                  "cluster_members_before_cutoff": clusters.get(cid, -1), **hist.get(a, {})}
    return out


def stage_dates() -> dict:
    items = benchmark.load()
    accs = sorted({it.accession for it in items})
    d = protein_dates(accs)
    RES.mkdir(exist_ok=True)
    cols = ["accession", "upi", "first_seen", "uniprot_first", "swissprot_first", "sequence_versions",
            "current_sequence_since", "uniref50", "cluster_members_before_cutoff"]
    with open(RES / "proteins.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for a in accs:
            w.writerow([a] + [d[a].get(c, "") for c in cols[1:]])
    return d


# ----------------------------------------------------------------- answer leakage
def stage_own() -> list[dict]:
    items = benchmark.load()
    accs = sorted({it.accession for it in items})

    def fetch(acc):
        meta, text = uniprot.last_version_before(acc, CUT)
        old = parse(text) if text else None
        cur = parse(uniprot.entry_text(acc))
        return acc, (meta, old, cur)

    entries = dict(pmap(fetch, accs))
    rows = []
    for it in items:
        meta, old, cur = entries[it.accession]
        o = classify(it.task, it.choices, it.answer, old)
        c = classify(it.task, it.choices, it.answer, cur)
        rows.append({"id": it.id, "task": it.task, "accession": it.accession,
                     "old_version": meta["entryVersion"] if meta else "", "old_db": meta["database"] if meta else "",
                     "old_release": meta["firstReleaseDate"] if meta else "",
                     "own_status": o["status"], "own_supported": ",".join(o["supported"]), "own_lookup_acc": o["lookup_acc"],
                     "current_status": c["status"], "current_supported": ",".join(c["supported"]),
                     "sequence_matches_old": (old.sequence == it.sequence) if old else ""})
    with open(RES / "items.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return rows


# ----------------------------------------------------------------- summary
def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def summarise() -> dict:
    prot = {r["accession"]: r for r in read_tsv(RES / "proteins.tsv")}
    items = read_tsv(RES / "items.tsv")
    cut = benchmark.CUTOFF
    n = len(prot)
    s = {"cutoff": cut, "benchmark_commit": benchmark.COMMIT, "proteins": n, "questions": len(items)}
    fs = [r["first_seen"] for r in prot.values() if r["first_seen"]]
    s["sequence_first_seen_before_cutoff"] = sum(f < cut for f in fs) / len(fs)
    for y in ("2023", "2021", "2015", "2010"):
        s[f"sequence_first_seen_before_{y}"] = sum(f < f"{y}-01-01" for f in fs) / len(fs)
    s["median_first_seen"] = sorted(fs)[len(fs) // 2]
    s["accession_in_uniprotkb_before_cutoff"] = sum(r["uniprot_first"] and r["uniprot_first"] < cut for r in prot.values()) / n
    s["swissprot_first_after_cutoff"] = sum(r["swissprot_first"] >= cut for r in prot.values() if r["swissprot_first"]) / n
    s["cluster_has_member_before_cutoff"] = sum(int(r["cluster_members_before_cutoff"] or 0) > 0 for r in prot.values()) / n
    s["sequence_changed_since_cutoff"] = sum(r["current_sequence_since"] >= cut for r in prot.values() if r["current_sequence_since"]) / n
    first_seen = {a: r["first_seen"] for a, r in prot.items()}
    s["questions_on_sequences_before_cutoff"] = sum(first_seen.get(i["accession"], "9") < cut for i in items) / len(items)

    by = defaultdict(Counter)
    acc_sum = defaultdict(float)
    for i in items:
        by[i["task"]]["questions"] += 1
        by[i["task"]][f"own_{i['own_status']}"] += 1
        by[i["task"]][f"current_{i['current_status']}"] += 1
        acc_sum[i["task"]] += float(i["own_lookup_acc"])
    tasks = {}
    for t, c in sorted(by.items()):
        q = c["questions"]
        tasks[t] = {"questions": q,
                    "leaked": c["own_leaked"] / q, "ambiguous": c["own_ambiguous"] / q,
                    "misleading": c["own_misleading"] / q, "absent": c["own_absent"] / q, "no_entry": c["own_no_entry"] / q,
                    "lookup_accuracy": acc_sum[t] / q,
                    "calibration_current_leaked": c["current_leaked"] / q,
                    "calibration_current_misleading": c["current_misleading"] / q}
    s["tasks"] = tasks
    tot = Counter(i["own_status"] for i in items)
    s["leaked"] = tot["leaked"] / len(items)
    s["lookup_accuracy"] = sum(float(i["own_lookup_acc"]) for i in items) / len(items)
    s["calibration_current_leaked"] = sum(i["current_status"] == "leaked" for i in items) / len(items)
    (RES / "summary.json").write_text(json.dumps(s, indent=2))
    return s


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["dates", "own", "summary", "all"])
    a = ap.parse_args()
    if a.stage in ("dates", "all"):
        stage_dates()
    if a.stage in ("own", "all"):
        stage_own()
    if a.stage in ("summary", "all"):
        print(json.dumps(summarise(), indent=1))


if __name__ == "__main__":
    main()
