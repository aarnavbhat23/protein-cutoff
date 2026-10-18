"""Does sequence age change how 'predictable' a post-cutoff protein looks to a protein LM?

ESM-2 was trained on UniRef from release 2021_04. Pseudo-perplexity (mask one
position at a time) is computed for every benchmark protein and for NovelFam
proteins (sequence and whole 50% cluster first seen after 2023-01-01, no
Swiss-Prot homolog; built in next-token-bio). Groups:

  before_esm2   benchmark sequence first seen before 2021-04-01
  after_esm2    benchmark sequence first seen on/after 2021-04-01
  novelfam      sequence and cluster new, no Swiss-Prot relative

Uses fair-esm (`pip install fair-esm`); weights come from Meta's official URLs.

  python -m cutoffaudit.plm --model esm2_t6_8M_UR50D   # larger models: same command on a GPU
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import torch

from . import benchmark
from .uniprot import ROOT

ESM2_DATA = "2021-04-01"
NOVELFAM = ROOT / "data" / "novelfam" / "novelfam_v1.tsv"   # copied from next-token-bio (UniProt data, CC BY 4.0)


def device():
    import os
    if os.environ.get("CUTOFFAUDIT_DEVICE"):          # e.g. cpu, cuda
        return torch.device(os.environ["CUTOFFAUDIT_DEVICE"])
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


SAMPLE = 32    # masked positions per protein (seeded); an unbiased estimate of the per-residue mean


@torch.no_grad()
def pseudo_perplexity(model, alphabet, seq: str, dev, batch: int = 64, max_len: int = 1022,
                      sample: int = SAMPLE, seed: int = 0) -> tuple[float, int]:
    """exp(mean of -log P(x_i | x_-i)) over up to `sample` positions, one masked at a time."""
    seq = seq[:max_len]
    _, _, toks = alphabet.get_batch_converter()([("p", seq)])
    toks = toks.to(dev)
    L = len(seq)
    rng = np.random.default_rng(seed)
    positions = np.sort(rng.choice(L, size=min(L, sample), replace=False))
    nll = 0.0
    for s in range(0, len(positions), batch):
        pos = torch.from_numpy(positions[s: s + batch]).to(dev)
        x = toks.repeat(len(pos), 1)
        x[torch.arange(len(pos), device=dev), pos + 1] = alphabet.mask_idx  # +1 for <cls>
        logits = model(x)["logits"][torch.arange(len(pos), device=dev), pos + 1]
        lp = torch.log_softmax(logits.float(), -1)
        nll -= lp.gather(1, toks[0, pos + 1][:, None])[:, 0].sum().item()
    return math.exp(nll / len(positions)), L


def proteins() -> list[dict]:
    prot = {r["accession"]: r for r in csv.DictReader(open(ROOT / "results" / "proteins.tsv"), delimiter="\t")}
    seqs = {}
    for it in benchmark.load():
        seqs.setdefault(it.accession, it.sequence)
    rows = [{"accession": a, "sequence": s, "first_seen": prot.get(a, {}).get("first_seen", ""),
             "group": "before_esm2" if prot.get(a, {}).get("first_seen", "9") < ESM2_DATA else "after_esm2"}
            for a, s in seqs.items()]
    if NOVELFAM.exists():
        for r in csv.DictReader(open(NOVELFAM), delimiter="\t"):
            rows.append({"accession": r["Entry"], "sequence": r["Sequence"], "first_seen": r["first_seen"], "group": "novelfam"})
    return rows


def bootstrap_diff(a, b, B=4000, seed=0):
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a), np.asarray(b)
    d = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(B)]
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def run(name: str) -> dict:
    import esm
    model, alphabet = getattr(esm.pretrained, name)()
    dev = device()
    model = model.eval().to(dev)
    out_tsv = ROOT / "results" / f"plm_{name}.tsv"
    done = {}
    if out_tsv.exists():
        done = {r["accession"]: r for r in csv.DictReader(open(out_tsv), delimiter="\t")}
    rows = proteins()
    with open(out_tsv, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "group", "first_seen", "length", "pseudo_perplexity"], delimiter="\t")
        if not done:
            w.writeheader()
        for i, r in enumerate(rows):
            if r["accession"] in done:
                continue
            ppl, L = pseudo_perplexity(model, alphabet, r["sequence"], dev)
            w.writerow({"accession": r["accession"], "group": r["group"], "first_seen": r["first_seen"], "length": L,
                        "pseudo_perplexity": f"{ppl:.4f}"})
            fh.flush()
            if i % 50 == 0:
                print(f"{i}/{len(rows)}", flush=True)
    res = [r for r in csv.DictReader(open(out_tsv), delimiter="\t")]
    groups = {}
    for g in ("before_esm2", "after_esm2", "novelfam"):
        v = sorted(float(r["pseudo_perplexity"]) for r in res if r["group"] == g)
        if v:
            groups[g] = {"n": len(v), "median": v[len(v) // 2], "q25": v[len(v) // 4], "q75": v[3 * len(v) // 4]}
    val = lambda g: [float(r["pseudo_perplexity"]) for r in res if r["group"] == g]
    s = {"model": name, "esm2_training_data": f"UniRef {ESM2_DATA[:7]}", "groups": groups}
    if groups.get("before_esm2") and groups.get("novelfam"):
        s["median_diff_before_minus_novelfam_ci95"] = bootstrap_diff(val("before_esm2"), val("novelfam"))
    if groups.get("before_esm2") and groups.get("after_esm2"):
        s["median_diff_before_minus_after_ci95"] = bootstrap_diff(val("before_esm2"), val("after_esm2"))
    (ROOT / "results" / f"plm_{name}.json").write_text(json.dumps(s, indent=2))
    return s


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="esm2_t6_8M_UR50D")
    print(json.dumps(run(ap.parse_args().model), indent=1))
