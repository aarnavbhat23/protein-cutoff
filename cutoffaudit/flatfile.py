"""Parse the parts of a UniProtKB flat-file entry the audit compares against.

Handles both feature-table layouts: the current one (release 2019_11 onwards,
`FT   ACT_SITE        154` plus `/note=` qualifiers) and the older column
layout (`FT   ACT_SITE    154    154       Description.`).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

EVIDENCE = re.compile(r"\{ECO:[^}]*\}")


@dataclass
class Entry:
    accessions: list[str] = field(default_factory=list)
    reviewed: bool = False
    go: list[tuple[str, str, str]] = field(default_factory=list)      # (GO id, aspect F/P/C, term name)
    ec: set[str] = field(default_factory=set)
    cc: dict[str, list[str]] = field(default_factory=dict)            # topic -> texts
    keywords: list[str] = field(default_factory=list)
    features: list[tuple[str, str, str]] = field(default_factory=list)  # (type, location, note)
    pubmed: list[str] = field(default_factory=list)
    sequence: str = ""

    # ------------------------------------------------------------- derived views
    def go_names(self, aspect: str) -> set[str]:
        return {norm(n) for _, a, n in self.go if a == aspect}

    def cofactors(self) -> set[str]:
        return {norm(m) for t in self.cc.get("COFACTOR", []) for m in re.findall(r"Name=([^;]+);", t)}

    def reactions(self) -> set[str]:
        return {norm(m) for t in self.cc.get("CATALYTIC ACTIVITY", []) for m in re.findall(r"Reaction=(.+?);(?:\s|$)", t)}

    def pathways(self) -> list[str]:
        return [norm(t) for t in self.cc.get("PATHWAY", [])]

    def feature_locations(self, ftype: str) -> list[str]:
        return [loc for t, loc, _ in self.features if t == ftype]

    def optimum(self, what: str) -> list[str]:
        """Values in 'Optimum pH is ...' / 'Optimum temperature is ...' sentences."""
        out = []
        for t in self.cc.get("BIOPHYSICOCHEMICAL PROPERTIES", []):
            for m in re.finditer(rf"[Oo]ptimum {what} (?:is|are|of)\s+([^.;]*(?:\.\d[^.;]*)?)", t):
                out.append(m.group(1))
        return out


def norm(s: str) -> str:
    s = EVIDENCE.sub(" ", s)
    s = s.replace("’", "'").lower()
    s = re.sub(r"\s+", " ", s)
    return s.strip(" .;")


def parse(text: str) -> Entry:
    e = Entry()
    cc_topic, cc_buf = None, []
    ft_cur: list | None = None
    in_seq = False
    seq = []

    def flush_cc():
        nonlocal cc_topic, cc_buf
        if cc_topic:
            e.cc.setdefault(cc_topic, []).append(" ".join(cc_buf).strip())
        cc_topic, cc_buf = None, []

    for line in text.splitlines():
        code, body = line[:2], line[5:]
        if in_seq:
            if line.startswith("//"):
                in_seq = False
            else:
                seq.append(re.sub(r"[^A-Z]", "", line))
            continue
        if code == "ID":
            e.reviewed = "Reviewed;" in line
        elif code == "AC":
            e.accessions += [a.strip() for a in body.split(";") if a.strip()]
        elif code == "DE":
            e.ec.update(re.findall(r"EC=(\d+\.[\d\-]+\.[\d\-]+\.[n\d\-]+)", body))
        elif code == "CC":
            if body.startswith("-!- "):
                flush_cc()
                topic, _, rest = body[4:].partition(":")
                cc_topic, cc_buf = topic.strip(), [rest.strip()]
            elif body.startswith("---") or body.startswith("Copyrighted") or body.startswith("Distributed"):
                flush_cc()
            elif cc_topic:
                cc_buf.append(body.strip())
        elif code == "DR" and body.startswith("GO;"):
            parts = [p.strip() for p in body.split(";")]
            if len(parts) >= 3 and ":" in parts[2]:
                aspect, _, name = parts[2].partition(":")
                e.go.append((parts[1], aspect, name))
        elif code == "KW":
            e.keywords += [norm(k) for k in body.rstrip(".").split(";") if k.strip()]
        elif code == "RX":
            e.pubmed += re.findall(r"PubMed=(\d+)", body)
        elif code == "FT":
            if len(line) > 5 and line[5] != " ":       # a new feature: type, location(s), [description]
                tok = line[5:].split()
                ftype = tok[0]
                pos = re.compile(r"^[<>?]?\d+$")
                if len(tok) >= 3 and pos.match(tok[1]) and pos.match(tok[2]):   # old layout: from, to, text
                    loc = tok[1] if tok[1] == tok[2] else f"{tok[1]}..{tok[2]}"
                    note = " ".join(tok[3:])
                else:                                                         # current layout
                    loc, note = (tok[1] if len(tok) > 1 else ""), " ".join(tok[2:])
                ft_cur = [ftype, loc, note.rstrip(".")]
                e.features.append(ft_cur)
            elif ft_cur is not None:                    # qualifier / continuation
                q = line[5:].strip()
                mm = re.match(r'^/note="?(.*?)"?$', q)
                if mm:
                    ft_cur[2] = (ft_cur[2] + " " + mm.group(1)).strip()
        elif code == "SQ":
            in_seq = True
        if code != "CC" and cc_topic:
            flush_cc()
    flush_cc()
    e.ec.update(m for t in e.cc.get("CATALYTIC ACTIVITY", []) for m in re.findall(r"EC=(\d+\.[\d\-]+\.[\d\-]+\.[n\d\-]+)", t))
    e.features = [tuple(f) for f in e.features]
    e.sequence = "".join(seq)
    return e
