"""Small UniProt REST client: retries, backoff and an on-disk cache.

Everything the audit needs from UniProt goes through here, so a rerun after an
interruption refetches nothing it already has.
"""

from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache"
API = "https://rest.uniprot.org"
UA = {"User-Agent": "protein-cutoff-audit (github.com/aarnavbhat23)"}


def http(url: str, data: bytes | None = None, tries: int = 7) -> tuple[bytes, dict]:
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=UA), timeout=180) as r:
                return r.read(), {k.lower(): v for k, v in r.headers.items()}
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504) or i == tries - 1:
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError):
            if i == tries - 1:
                raise
        time.sleep(min(60, 2 ** i))
    raise RuntimeError("unreachable")


def cached(kind: str, key: str, fetch) -> bytes:
    """Return fetch() for (kind, key), memoised under data/cache/<kind>/."""
    h = hashlib.sha1(key.encode()).hexdigest()
    path = CACHE / kind / h[:2] / f"{h}.gz"
    if path.exists():
        return gzip.decompress(path.read_bytes())
    body = fetch()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_bytes(gzip.compress(body))
    tmp.replace(path)
    return body


def parse_tsv(text: str) -> list[dict]:
    lines = text.rstrip("\n").split("\n")
    if len(lines) < 2:
        return []
    keys = [k.strip() for k in lines[0].split("\t")]
    return [dict(zip(keys, l.split("\t"))) for l in lines[1:]]


def search(q: str, fields: str, size: int = 500) -> list[dict]:
    url = f"{API}/uniprotkb/search?" + urllib.parse.urlencode({"query": q, "format": "tsv", "fields": fields, "size": size})
    return parse_tsv(cached("search", url, lambda: http(url)[0]).decode())


def stream(q: str, fields: str) -> list[dict]:
    url = f"{API}/uniprotkb/stream?" + urllib.parse.urlencode({"query": q, "format": "tsv", "fields": fields, "compressed": "true"})
    body = cached("stream", url, lambda: gzip.decompress(http(url)[0]))
    return parse_tsv(body.decode())


def count(q: str) -> int:
    url = f"{API}/uniprotkb/search?" + urllib.parse.urlencode({"query": q, "format": "tsv", "fields": "accession", "size": 1})
    body = cached("count", url, lambda: http(url)[1].get("x-total-results", "0").encode())
    return int(body.decode())


def idmap(ids: list[str], to: str, fields: str, frm: str = "UniProtKB_AC-ID") -> list[dict]:
    """UniProt ID mapping, batched; results cached per batch."""
    kind = {"UniRef50": "uniref", "UniParc": "uniparc", "UniProtKB": "uniprotkb"}[to]
    out = []
    for i in range(0, len(ids), 20000):
        batch = ids[i: i + 20000]

        def fetch(batch=batch):
            data = urllib.parse.urlencode({"ids": ",".join(batch), "from": frm, "to": to}).encode()
            job = json.loads(http(f"{API}/idmapping/run", data=data)[0])["jobId"]
            for _ in range(900):
                st = json.loads(http(f"{API}/idmapping/status/{job}")[0])
                if st.get("jobStatus") in (None, "FINISHED") or "results" in st:
                    break
                if st.get("jobStatus") in ("ERROR", "FAILED"):
                    raise RuntimeError(f"idmapping {job}: {st}")
                time.sleep(3)
            return http(f"{API}/idmapping/{kind}/results/stream/{job}?format=tsv&fields={fields}")[0]

        out += parse_tsv(cached("idmap", f"{frm}>{to}:{fields}:" + ",".join(batch), fetch).decode())
    return out


def uniparc_first_seen(upis: list[str]) -> dict[str, str]:
    out = {}
    for i in range(0, len(upis), 100):
        q = " OR ".join(f"upi:{u}" for u in upis[i: i + 100])
        url = f"{API}/uniparc/search?" + urllib.parse.urlencode({"query": q, "format": "tsv", "fields": "upi,first_seen", "size": 500})
        for r in parse_tsv(cached("uniparc", url, lambda: http(url)[0]).decode()):
            out[r["Entry"]] = r["First seen"]
    return out


def cluster_members(cid: str, cap: int = 5000) -> list[str] | None:
    """Member ids of a UniRef cluster; None if more than `cap`."""
    def fetch():
        url, members = f"{API}/uniref/{cid}/members?format=list&size=500", []
        while url:
            body, headers = http(url)
            members += [m for m in body.decode().split("\n") if m]
            if len(members) > cap:
                return b"__TOO_LARGE__"
            link = headers.get("link", "")
            url = link[link.find("<") + 1: link.find(">")] if 'rel="next"' in link else None
        return "\n".join(members).encode()
    body = cached("members", f"{cid}:{cap}", fetch).decode()
    return None if body == "__TOO_LARGE__" else [m for m in body.split("\n") if m]


# ----------------------------------------------------------------- UniSave (entry history)
def unisave_history(acc: str) -> list[dict]:
    """All archived versions of a UniProtKB entry, newest first."""
    url = f"{API}/unisave/{acc}?format=json"
    try:
        body = cached("unisave_hist", url, lambda: http(url)[0])
    except urllib.error.HTTPError as e:
        if e.code in (400, 404):
            return []
        raise
    return json.loads(body).get("results", [])


def release_date(entry: dict) -> dt.date:
    return dt.datetime.strptime(entry["firstReleaseDate"], "%d-%b-%Y").date()


def unisave_version_text(acc: str, version: int) -> str:
    url = f"{API}/unisave/{acc}?format=txt&versions={version}"
    return cached("unisave_txt", url, lambda: http(url)[0]).decode()


def last_version_before(acc: str, cutoff: dt.date) -> tuple[dict | None, str | None]:
    """The newest archived version first released before `cutoff`, and its flat-file text."""
    pre = [e for e in unisave_history(acc) if release_date(e) < cutoff]
    if not pre:
        return None, None
    e = max(pre, key=lambda e: e["entryVersion"])
    return e, unisave_version_text(acc, e["entryVersion"])


def entry_text(acc: str) -> str:
    """Current UniProtKB flat file (cached)."""
    url = f"{API}/uniprotkb/{acc}.txt"
    return cached("entry_txt", url, lambda: http(url)[0]).decode()
