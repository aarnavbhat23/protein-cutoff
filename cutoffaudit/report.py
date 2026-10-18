"""Figures and the results summary.

  python -m cutoffaudit.report    # results/figures/*.png, results/summary.md
"""

from __future__ import annotations

import csv
import json
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import benchmark
from .uniprot import ROOT

RES = ROOT / "results"
FIG = RES / "figures"
# reference palette (validated): slot 1 blue, slot 2 orange; neutrals for text and grids
BLUE, ORANGE = "#2a78d6", "#eb6834"
RAMP = ["#184f95", "#2a78d6", "#86b6ef", "#c3c2b7"]   # leaked -> ambiguous -> misleading (validated ordinal ramp), absent = gray
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
    "axes.titlesize": 11, "axes.titleweight": "bold", "legend.frameon": False, "lines.linewidth": 2,
})
LABEL = {"GO_molecular_function": "GO molecular function", "GO_biological_process": "GO biological process",
         "GO_cellular_component": "GO cellular component", "EC_number": "EC number", "cofactor": "Cofactor",
         "catalytic_activity": "Catalytic activity", "pathway": "Pathway", "active_site": "Active site",
         "motif_position": "Motif position", "transmembrane": "Transmembrane", "ph": "Optimal pH",
         "temperature": "Optimal temperature"}


def load(name):
    p = RES / name
    if not p.exists():
        return None
    return json.loads(p.read_text()) if p.suffix == ".json" else list(csv.DictReader(open(p), delimiter="\t"))


def fig_first_seen(prot):
    years = [int(r["first_seen"][:4]) for r in prot if r["first_seen"]]
    fig, ax = plt.subplots(figsize=(8, 3.3))
    lo, hi = min(years), max(years)
    cnt = Counter(years)
    xs = np.arange(lo, hi + 1)
    ax.bar(xs, [cnt.get(x, 0) for x in xs], width=0.8, color=[BLUE if x < 2025 else ORANGE for x in xs])
    ax.axvline(2024.5, color=INK2, lw=1, ls="--")
    pre = sum(y < 2025 for y in years) / len(years)
    ax.text(2024.3, ax.get_ylim()[1] * 0.92, f"benchmark cutoff\n{pre:.0%} of sequences are older",
            ha="right", va="top", fontsize=9, color=INK2)
    ax.set_xlabel("year the sequence was first seen anywhere (UniParc)")
    ax.set_ylabel("benchmark proteins")
    ax.set_title(f"{len(years)} benchmark proteins, all 'created' in 2025: when their sequences first appeared")
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    fig.savefig(FIG / "first_seen.png", dpi=160)
    plt.close(fig)


def fig_own_entry(summary):
    tasks = sorted(summary["tasks"], key=lambda t: -summary["tasks"][t]["leaked"])
    parts = [("leaked", "answer, and only the answer, in the pre-cutoff entry"),
             ("ambiguous", "answer plus a wrong option"), ("misleading", "only a wrong option"),
             ("absent", "no option / no pre-cutoff entry")]
    fig, ax = plt.subplots(figsize=(9, 4.4))
    y = np.arange(len(tasks))[::-1]
    left = np.zeros(len(tasks))
    for (k, lab), col in zip(parts, RAMP):
        v = np.array([summary["tasks"][t][k] + (summary["tasks"][t]["no_entry"] if k == "absent" else 0) for t in tasks]) * 100
        ax.barh(y, v, left=left, color=col, height=0.7, label=lab, edgecolor=SURF, linewidth=1.5)
        left += v
    for yi, t in zip(y, tasks):
        v = summary["tasks"][t]["leaked"] * 100
        inside = v > 8
        ax.text(v - 1.2 if inside else v + 1.2, yi, f"{v:.0f}%", va="center", fontsize=8,
                color=SURF if inside else INK2, ha="right" if inside else "left")
    ax.set_yticks(y, [LABEL[t] for t in tasks])
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of questions")
    ax.set_title("Was the answer already in the protein's own UniProt entry before 2025?")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=2, fontsize=8)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(FIG / "own_entry_leakage.png", dpi=160)
    plt.close(fig)


def fig_lookup(hom, sens, reported):
    tasks = list(hom["tasks"])
    models = reported["models"]
    best = {t: max((m[t], n) for n, m in models.items() if t in m) for t in tasks}
    fig, ax = plt.subplots(figsize=(8.5, 4))
    y = np.arange(len(tasks))[::-1]
    for yi, t in zip(y, tasks):
        a, b, r = hom["tasks"][t]["accuracy"] * 100, sens["identity<0.50"][t] * 100, best[t][0]
        ax.plot([min(a, b, r), max(a, b, r)], [yi, yi], color=GRID, lw=2, zorder=1)
        ax.scatter(a, yi, s=60, color=BLUE, zorder=3, label="lookup: Swiss-Prot 2024_06 relatives" if yi == y[0] else None)
        ax.scatter(b, yi, s=60, facecolor=SURF, edgecolor=BLUE, linewidth=2, zorder=3,
                   label="lookup: relatives < 50% identity only" if yi == y[0] else None)
        ax.scatter(r, yi, s=60, marker="D", color=ORANGE, zorder=3, label="best model in the paper" if yi == y[0] else None)
    ax.axvline(25, color=MUTED, lw=1, ls=":")
    ax.text(25.8, y[-1] - 0.45, "chance (4 options)", fontsize=8, color=MUTED, va="center")
    ax.set_ylim(y[-1] - 0.8, y[0] + 0.6)
    ax.set_yticks(y, [f"{LABEL[t]}\nbest model: {best[t][1]}" for t in tasks], fontsize=8)
    ax.set_xlim(0, 100)
    ax.set_xlabel("accuracy (%) on LiveProteinBench")
    ax.set_title("A no-learning lookup in pre-cutoff Swiss-Prot vs the best reported model")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.38), ncol=3, fontsize=8)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(FIG / "lookup_vs_models.png", dpi=160)
    plt.close(fig)


def fig_survey(survey):
    years = sorted(survey["years"])
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    for key, col, lab in (("reviewed", BLUE, "Swiss-Prot (reviewed)"), ("unreviewed", ORANGE, "TrEMBL (unreviewed)")):
        xs = [int(y) for y in years if key in survey["years"][y]]
        ys = [survey["years"][str(x)][key]["frac_sequence_older_than_1y"] * 100 for x in xs]
        ax.plot(xs, ys, color=col, marker="o", ms=4, label=lab)
        ax.text(xs[-1] + 0.2, ys[-1], lab.split(" ")[0], color=INK2, fontsize=8, va="center")
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of entries whose sequence is\n> 1 year older than 'date created'")
    ax.set_xlabel("year of UniProt 'date created'")
    ax.set_title("UniProt's creation date is a curation date, not a sequence date")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "survey.png", dpi=160)
    plt.close(fig)


def md_table(rows, header):
    return ("| " + " | ".join(header) + " |\n|" + "---|" * len(header) + "\n" +
            "".join("| " + " | ".join(map(str, r)) + " |\n" for r in rows))


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    prot, summary = load("proteins.tsv"), load("summary.json")
    hom, sens, reported = load("homolog.json"), load("homolog_sensitivity.json"), load("liveproteinbench_reported.json")
    survey, lit = load("survey.json"), load("literature.json")
    lines = ["# Results (generated by `python -m cutoffaudit.report`)\n",
             f"Benchmark: LiveProteinBench @ `{benchmark.COMMIT[:12]}`, cutoff {benchmark.CUTOFF}.\n"]
    if prot:
        fig_first_seen(prot)
    if summary:
        fig_own_entry(summary)
        lines.append("## Dates\n")
        keys = ["proteins", "questions", "sequence_first_seen_before_cutoff", "sequence_first_seen_before_2023",
                "sequence_first_seen_before_2015", "median_first_seen", "accession_in_uniprotkb_before_cutoff",
                "swissprot_first_after_cutoff", "cluster_has_member_before_cutoff", "questions_on_sequences_before_cutoff"]
        lines.append(md_table([[k, f"{summary[k]:.1%}" if isinstance(summary[k], float) else summary[k]] for k in keys], ["measure", "value"]))
        lines.append("\n## Own pre-cutoff entry (UniSave), per task\n")
        rows = [[LABEL[t], v["questions"], f"{v['leaked']:.1%}", f"{v['ambiguous']:.1%}", f"{v['misleading']:.1%}",
                 f"{v['absent'] + v['no_entry']:.1%}", f"{v['lookup_accuracy']:.1%}", f"{v['calibration_current_leaked']:.1%}"]
                for t, v in summary["tasks"].items()]
        lines.append(md_table(rows, ["task", "Q", "leaked", "ambiguous", "misleading", "absent/none", "own-entry lookup acc.",
                                     "calibration: current entry leaked"]))
        lines.append(f"\nAll questions: leaked {summary['leaked']:.1%}; own-entry lookup accuracy {summary['lookup_accuracy']:.1%}; "
                     f"calibration on current entries {summary['calibration_current_leaked']:.1%}.\n")
    if hom and sens and reported:
        fig_lookup(hom, sens, reported)
        lines.append("\n## Lookup in Swiss-Prot 2024_06 vs reported models\n")
        models = reported["models"]
        rows = []
        for t, v in hom["tasks"].items():
            b = max((m[t], n) for n, m in models.items())
            rows.append([LABEL[t], v["questions"], f"{v['accuracy']:.1%}", f"{sens['identity<0.50'][t]:.1%}", f"{b[0]:.1f}% ({b[1]})"])
        lines.append(md_table(rows, ["task", "Q", "lookup", "lookup, < 50% identity", "best reported model"]))
    if survey:
        fig_survey(survey)
        lines.append("\n## UniProt-wide: date created vs first seen\n")
        rows = [[y, *(f"{v[k]['frac_sequence_older_than_1y']:.0%} (n={v[k]['n']}, median lag {v[k]['median_lag_years']:.1f} y)"
                      if k in v else "-" for k in ("reviewed", "unreviewed"))] for y, v in sorted(survey["years"].items())]
        lines.append(md_table(rows, ["year created", "reviewed: sequence > 1 y older", "unreviewed: sequence > 1 y older"]))
    if lit:
        lines.append("\n## Literature cited in the current entries\n")
        lines.append(md_table([[k, f"{v:.1%}" if isinstance(v, float) else v] for k, v in lit.items()], ["measure", "value"]))
    for name in sorted(p.name for p in RES.glob("plm_*.json")):
        p = load(name)
        lines.append(f"\n## {p['model']}: pseudo-perplexity by sequence age\n")
        lines.append(md_table([[g, v["n"], f"{v['median']:.2f}", f"[{v['q25']:.2f}, {v['q75']:.2f}]"] for g, v in p["groups"].items()],
                              ["group", "n", "median", "IQR"]))
        for k in ("median_diff_before_minus_novelfam_ci95", "median_diff_before_minus_after_ci95"):
            if k in p:
                lines.append(f"\n{k}: [{p[k][0]:.2f}, {p[k][1]:.2f}]\n")
    (RES / "summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
