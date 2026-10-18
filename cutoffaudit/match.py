"""Does an entry support a multiple-choice option? One rule per benchmark task.

Each rule compares an option with the specific annotation the task is about
(the GO terms of one aspect, the EC numbers, the active-site positions, ...),
never with the entry's free text as a whole.
"""

from __future__ import annotations

import re

from .flatfile import Entry, norm

GO_ASPECT = {"GO_molecular_function": "F", "GO_biological_process": "P", "GO_cellular_component": "C"}


def _nums(s: str) -> list[float]:
    return [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", s)]


def _value_supported(option: str, texts: list[str]) -> bool:
    """A pH/temperature option (a number or a range) appears as a stated optimum."""
    opt = option.strip()
    for t in texts:
        if re.search(rf"(?<![\d.]){re.escape(opt)}(?![\d.])", t.replace(" to ", "-")):
            return True
        v = _nums(opt)
        if len(v) == 1 and v[0] in _nums(t.split(" degree")[0]):
            return True
    return False


def supports(task: str, option: str, e: Entry) -> bool:
    o = norm(option)
    if not o:
        return False
    if task in GO_ASPECT:
        return o in e.go_names(GO_ASPECT[task])
    if task == "EC_number":
        ec = re.search(r"\d+\.[\d\-]+\.[\d\-]+\.[n\d\-]+", option)
        return bool(ec) and ec.group(0) in e.ec
    if task == "cofactor":
        return o in e.cofactors()
    if task == "catalytic_activity":
        return o in e.reactions()
    if task == "pathway":
        o = re.sub(r"[\s.]+$", "", o)
        return len(o) >= 4 and (any(o in p for p in e.pathways()) or o in e.keywords)
    if task == "active_site":
        m = re.search(r"position\s+(\d+)", o)
        return bool(m) and m.group(1) in e.feature_locations("ACT_SITE")
    if task == "motif_position":
        return option.strip() in e.feature_locations("MOTIF")
    if task == "transmembrane":
        spans = {s.strip() for s in option.split(";") if s.strip()}
        return bool(spans) and spans == set(e.feature_locations("TRANSMEM"))
    if task == "ph":
        return _value_supported(option, e.optimum("pH"))
    if task == "temperature":
        return _value_supported(option, e.optimum("temperature"))
    raise ValueError(f"unknown task {task}")


def classify(task: str, choices: dict[str, str], answer: str, e: Entry | None) -> dict:
    """leaked / ambiguous / misleading / absent, and the lookup baseline's expected accuracy."""
    if e is None:
        return {"status": "no_entry", "supported": [], "lookup_acc": 1 / max(1, len(choices))}
    sup = sorted(k for k, v in choices.items() if supports(task, v, e))
    if answer in sup and len(sup) == 1:
        status, acc = "leaked", 1.0
    elif answer in sup:
        status, acc = "ambiguous", 1 / len(sup)
    elif sup:
        status, acc = "misleading", 0.0
    else:
        status, acc = "absent", 1 / max(1, len(choices))
    return {"status": status, "supported": sup, "lookup_acc": acc}
