"""LiveProteinBench items, fetched at a pinned commit.

The benchmark repository has no licence, so its files are downloaded into
data/raw/ (git-ignored) and never redistributed; this project publishes only
accession-level results.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from .uniprot import ROOT, http

REPO = "Rongdingyi/LiveProteinBench"
COMMIT = "8b9710932511f816651fc1a7e8a1bc85c3c77767"   # HEAD on 2026-10-17 (committed 2026-05-07)
TASKS = ["EC_number", "GO_biological_process", "GO_cellular_component", "GO_molecular_function", "active_site",
         "catalytic_activity", "cofactor", "motif_position", "pathway", "ph", "temperature", "transmembrane"]
CUTOFF = "2025-01-01"   # "protein records released after January 1, 2025" (arXiv:2512.22257)


@dataclass(frozen=True)
class Item:
    task: str
    idx: int
    accession: str
    choices: dict
    answer: str
    sequence: str

    @property
    def id(self) -> str:
        return f"{self.task}/{self.idx}"


def raw_dir(commit: str = COMMIT) -> Path:
    return ROOT / "data" / "raw" / "liveproteinbench" / commit[:12]


def fetch(commit: str = COMMIT) -> Path:
    d = raw_dir(commit)
    d.mkdir(parents=True, exist_ok=True)
    for t in TASKS:
        p = d / f"{t}.json"
        if not p.exists():
            body, _ = http(f"https://raw.githubusercontent.com/{REPO}/{commit}/dataset/QA/{t}.json")
            p.write_bytes(body)
    return d


def load(commit: str = COMMIT) -> list[Item]:
    d = fetch(commit)
    items = []
    for t in TASKS:
        for i, q in enumerate(json.loads((d / f"{t}.json").read_text())):
            acc = q["Protein Image"][0].split("/")[1]
            items.append(Item(t, i, acc, dict(q["Multiple Choices"]), q["Answer"], q.get("Protein Sequence", "")))
    return items
