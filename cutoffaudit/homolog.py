"""Could the questions be answered from Swiss-Prot as it stood before the cutoff?

A retrieval baseline with no learning: search each benchmark protein against
UniProtKB/Swiss-Prot release 2024_06 (27 Nov 2024, the last release before
2025-01-01) and answer each multiple-choice question with the option best
supported by the annotations of its close relatives there.

  index   parse uniprot_sprot.dat.gz (2024_06) into FASTA + annotation index
  search  MMseqs2 (-s 7.5), all hits with E < 1e-3
  score   per question: options supported by hits with >= MIN_ID identity over
          >= MIN_COV of the query, weighted by identity; argmax, ties split
          evenly, no support = uniform guess. Positional tasks (active site,
          motif, transmembrane) would need alignment-mapped coordinates and
          pH/temperature are rarely annotated, so this baseline covers the
          seven annotation tasks only.

  python -m cutoffaudit.homolog index
  python -m cutoffaudit.homolog search
  python -m cutoffaudit.homolog score
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from . import benchmark
from .flatfile import Entry, parse
from .match import supports
from .uniprot import ROOT

SPROT = ROOT / "data" / "raw" / "sprot_2024_06"
IDX = ROOT / "data" / "cache" / "sprot_2024_06"
TASKS = ["GO_molecular_function", "GO_biological_process", "GO_cellular_component", "EC_number",
         "cofactor", "catalytic_activity", "pathway"]
MIN_ID, MIN_COV, TOP = 0.30, 0.50, 10


def entries(path: Path):
    buf = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            buf.append(line)
            if line.startswith("//"):
                yield "".join(buf)
                buf = []


def build_index() -> dict:
    IDX.mkdir(parents=True, exist_ok=True)
    ann, n = {}, 0
    with open(IDX / "sprot.fasta", "w") as fa:
        for text in entries(SPROT / "uniprot_sprot.dat.gz"):
            e = parse(text)
            acc = e.accessions[0]
            fa.write(f">{acc}\n{e.sequence}\n")
            ann[acc] = {"go": e.go, "ec": sorted(e.ec), "cc": {k: v for k, v in e.cc.items()
                        if k in ("COFACTOR", "CATALYTIC ACTIVITY", "PATHWAY")}, "kw": e.keywords}
            n += 1
    with gzip.open(IDX / "annotations.json.gz", "wt") as fh:
        json.dump(ann, fh)
    meta = {"release": "2024_06", "entries": n}
    (IDX / "meta.json").write_text(json.dumps(meta))
    return meta


def as_entry(a: dict) -> Entry:
    e = Entry()
    e.go = [tuple(g) for g in a["go"]]
    e.ec = set(a["ec"])
    e.cc = a["cc"]
    e.keywords = a["kw"]
    return e


def run_search(threads: int = 4) -> Path:
    items = benchmark.load()
    seqs = {}
    for it in items:
        seqs.setdefault(it.accession, it.sequence)
    out = IDX / "hits.m8"
    exe = shutil.which("mmseqs")
    with tempfile.TemporaryDirectory(prefix="mmseqs_") as tmp:   # mmseqs breaks on paths with spaces
        with open(Path(tmp) / "q.fasta", "w") as fh:
            for a, s in seqs.items():
                fh.write(f">{a}\n{s}\n")
        shutil.copy(IDX / "sprot.fasta", Path(tmp) / "t.fasta")
        args = ["easy-search", "q.fasta", "t.fasta", "hits.m8", "tmp", "-s", "7.5", "-e", "1e-3", "--max-seqs", "300",
                "--format-output", "query,target,fident,qcov,tcov,evalue,bits", "--threads", str(threads)]
        subprocess.run([exe] + args, check=True, cwd=tmp, capture_output=True)
        shutil.copy(Path(tmp) / "hits.m8", out)
    return out


def score(max_id: float = 1.01, write: bool = True) -> dict:
    with gzip.open(IDX / "annotations.json.gz", "rt") as fh:
        ann = json.load(fh)
    hits = defaultdict(list)
    for line in (IDX / "hits.m8").read_text().splitlines():
        q, t, fid, qcov, tcov, ev, bits = line.split("\t")
        hits[q].append((float(fid), float(qcov), float(tcov), t))
    for q in hits:
        hits[q].sort(reverse=True)
    rows = []
    for it in benchmark.load():
        if it.task not in TASKS:
            continue
        close = [(fid, t) for fid, qc, tc, t in hits.get(it.accession, [])
                 if MIN_ID <= fid < max_id and qc >= MIN_COV and t in ann][:TOP]
        identical = any(fid >= 0.999 and qc >= 0.999 and tc >= 0.999 for fid, qc, tc, _ in hits.get(it.accession, []))
        support = {k: 0.0 for k in it.choices}
        for fid, t in close:
            e = as_entry(ann[t])
            for k, v in it.choices.items():
                if supports(it.task, v, e):
                    support[k] = max(support[k], fid)
        best = max(support.values())
        top = [k for k, v in support.items() if v == best]
        acc = (1.0 / len(top) if it.answer in top else 0.0) if best > 0 else 1.0 / len(it.choices)
        rows.append({"id": it.id, "task": it.task, "accession": it.accession, "n_close_hits": len(close),
                     "best_identity": close[0][0] if close else 0.0, "identical_in_sprot_2024_06": identical,
                     "answered": best > 0, "expected_accuracy": acc})
    if write:
        with open(ROOT / "results" / "homolog_items.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    out = {"swissprot_release": "2024_06", "min_identity": MIN_ID, "max_identity": max_id,
           "min_query_coverage": MIN_COV, "tasks": {}}
    for t in TASKS:
        r = [x for x in rows if x["task"] == t]
        out["tasks"][t] = {"questions": len(r), "accuracy": sum(x["expected_accuracy"] for x in r) / len(r),
                           "answered": sum(x["answered"] for x in r) / len(r),
                           "accuracy_when_answered": (sum(x["expected_accuracy"] for x in r if x["answered"]) /
                                                      max(1, sum(x["answered"] for x in r)))}
    out["overall_accuracy"] = sum(x["expected_accuracy"] for x in rows) / len(rows)
    out["questions"] = len(rows)
    out["proteins_with_identical_sequence_in_sprot_2024_06"] = len({x["accession"] for x in rows if x["identical_in_sprot_2024_06"]})
    out["proteins_with_close_hit"] = len({x["accession"] for x in rows if x["n_close_hits"]})
    if write:
        (ROOT / "results" / "homolog.json").write_text(json.dumps(out, indent=2))
    return out


def sensitivity() -> dict:
    """Same baseline using only relatives below an identity ceiling (distant homologs only)."""
    res = {}
    for cap in (0.5, 0.7, 0.9, 1.01):
        r = score(max_id=cap, write=False)
        res[f"identity<{cap:.2f}"] = {"overall": r["overall_accuracy"],
                                      **{t: round(v["accuracy"], 4) for t, v in r["tasks"].items()}}
    (ROOT / "results" / "homolog_sensitivity.json").write_text(json.dumps(res, indent=2))
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["index", "search", "score", "sensitivity", "all"])
    a = ap.parse_args()
    if a.stage in ("index", "all"):
        print(build_index())
    if a.stage in ("search", "all"):
        print(run_search())
    if a.stage in ("score", "all"):
        print(json.dumps(score(), indent=1))
    if a.stage in ("sensitivity", "all"):
        print(json.dumps(sensitivity(), indent=1))


if __name__ == "__main__":
    main()
