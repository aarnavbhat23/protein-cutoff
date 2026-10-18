"""How old is the literature behind each 'post-cutoff' protein?

LiveProteinBench says it cross-validates entry creation dates with
publication dates. For every protein this collects the PubMed references cited
in its current UniProtKB entry and their publication dates (NCBI E-utilities).

  python -m cutoffaudit.literature     # results/literature.tsv, results/literature.json
"""

from __future__ import annotations

import concurrent.futures as cf
import csv
import json
import time
import urllib.parse

from . import benchmark, uniprot
from .flatfile import parse

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def pubmed_dates(pmids: list[str]) -> dict[str, str]:
    """PMID -> 'YYYY/MM/DD' (NCBI sortpubdate), batched and cached."""
    out = {}
    for i in range(0, len(pmids), 200):
        batch = pmids[i: i + 200]
        url = EUTILS + "?" + urllib.parse.urlencode({"db": "pubmed", "id": ",".join(batch), "retmode": "json"})

        def fetch(url=url):
            time.sleep(0.4)  # stay under NCBI's 3 requests/s without an API key
            return uniprot.http(url)[0]

        res = json.loads(uniprot.cached("pubmed", url, fetch))["result"]
        for p in batch:
            if p in res and "sortpubdate" in res[p]:
                out[p] = res[p]["sortpubdate"][:10]
    return out


def run() -> dict:
    items = benchmark.load()
    accs = sorted({it.accession for it in items})
    with cf.ThreadPoolExecutor(4) as ex:
        refs = dict(zip(accs, ex.map(lambda a: parse(uniprot.entry_text(a)).pubmed, accs)))
    dates = pubmed_dates(sorted({p for v in refs.values() for p in v}))
    cut = benchmark.CUTOFF.replace("-", "/")
    rows, newest = [], []
    for a in accs:
        ds = sorted(dates[p] for p in refs[a] if p in dates)
        rows.append({"accession": a, "n_pubmed": len(refs[a]), "oldest": ds[0] if ds else "", "newest": ds[-1] if ds else "",
                     "n_before_cutoff": sum(d < cut for d in ds), "n_after_cutoff": sum(d >= cut for d in ds)})
        if ds:
            newest.append(ds[-1])
    with open(uniprot.ROOT / "results" / "literature.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    cited = [r for r in rows if r["n_pubmed"]]
    s = {"proteins": len(rows), "with_pubmed_citations": len(cited),
         "all_citations_before_cutoff": sum(r["n_after_cutoff"] == 0 for r in cited) / len(cited),
         "any_citation_before_cutoff": sum(r["n_before_cutoff"] > 0 for r in cited) / len(cited),
         "median_newest_citation": sorted(newest)[len(newest) // 2] if newest else None,
         "median_oldest_citation": sorted(r["oldest"] for r in cited)[len(cited) // 2] if cited else None}
    (uniprot.ROOT / "results" / "literature.json").write_text(json.dumps(s, indent=2))
    return s


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
