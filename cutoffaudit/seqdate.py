"""Date UniProt proteins against a training cutoff.

Labels each accession:
  old_sequence             UniParc first saw the sequence before the cutoff
  new_sequence_old_family  sequence new, but its UniRef50 cluster has a
                           UniProtKB member created before the cutoff
  new                      sequence and its 50% neighbourhood are both new
  unknown                  no UniParc / UniRef50 mapping

UniProt's own 'date created' is not used: for Swiss-Prot entries it is the
curation date, typically years to decades after the sequence appeared.

  python -m cutoffaudit.seqdate --cutoff 2025-01-01 accessions.txt > dated.tsv
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys

from .audit import protein_dates

COLS = ["accession", "verdict", "first_seen", "uniprot_first", "swissprot_first", "uniref50",
        "cluster_members_before_cutoff", "upi"]


def verdict(r: dict, cutoff: str) -> str:
    if not r.get("first_seen"):
        return "unknown"
    if r["first_seen"] < cutoff:
        return "old_sequence"
    n = int(r.get("cluster_members_before_cutoff", -1))
    if n < 0:
        return "unknown"
    return "new_sequence_old_family" if n > 0 else "new"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accessions", help="file with one UniProtKB accession per line ('-' for stdin)")
    ap.add_argument("--cutoff", required=True, help="YYYY-MM-DD, e.g. a model's training-data date")
    a = ap.parse_args(argv)
    cut = dt.date.fromisoformat(a.cutoff)
    src = sys.stdin if a.accessions == "-" else open(a.accessions)
    accs = sorted({l.strip() for l in src if l.strip() and not l.startswith("#")})
    d = protein_dates(accs, cut)
    print("\t".join(COLS))
    for acc in accs:
        r = d[acc]
        print("\t".join([acc, verdict(r, a.cutoff)] + [str(r.get(c, "")) for c in COLS[2:]]))


if __name__ == "__main__":
    main()
