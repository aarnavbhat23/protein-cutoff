"""How well does UniProt's 'date created' date a sequence, in general?

For each year, a sample of UniProtKB entries whose `date_created` falls in that
year (reviewed and unreviewed separately), compared with the date UniParc first
saw the sequence. For a reviewed entry, `date_created` is the date it was
integrated into Swiss-Prot, which can be decades after the sequence appeared.

  python -m cutoffaudit.survey     # results/survey.json
"""

from __future__ import annotations

import datetime as dt
import json
import random

from . import uniprot

YEARS = list(range(2016, 2027))
PER_MONTH = 50   # first page of each month's results


def sample(reviewed: bool, year: int) -> list[dict]:
    rows = []
    for m in range(1, 13):
        a = dt.date(year, m, 1)
        b = (dt.date(year + 1, 1, 1) if m == 12 else dt.date(year, m + 1, 1)) - dt.timedelta(days=1)
        if a > dt.date.today():
            break
        q = f"(reviewed:{str(reviewed).lower()}) AND (date_created:[{a.isoformat()} TO {b.isoformat()}])"
        rows += uniprot.search(q, "accession,date_created", size=PER_MONTH)
    return rows


def run(seed: int = 0) -> dict:
    out = {"per_month_sample": PER_MONTH, "years": {}}
    for reviewed in (True, False):
        key = "reviewed" if reviewed else "unreviewed"
        allrows = {y: sample(reviewed, y) for y in YEARS}
        accs = sorted({r["Entry"] for rs in allrows.values() for r in rs})
        fs = {r["From"]: r["First seen"] for r in uniprot.idmap(accs, "UniParc", "first_seen")}
        for y, rs in allrows.items():
            lags = []
            for r in rs:
                f = fs.get(r["Entry"])
                if f:
                    lags.append((dt.date.fromisoformat(r["Date of creation"]) - dt.date.fromisoformat(f)).days / 365.25)
            if not lags:
                continue
            lags.sort()
            out["years"].setdefault(str(y), {})[key] = {
                "n": len(lags), "median_lag_years": lags[len(lags) // 2],
                "frac_sequence_older_than_1y": sum(l > 1 for l in lags) / len(lags),
                "frac_sequence_older_than_5y": sum(l > 5 for l in lags) / len(lags)}
    (uniprot.ROOT / "results" / "survey.json").write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
