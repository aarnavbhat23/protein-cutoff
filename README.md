# protein-cutoff-audit

[![tests](https://github.com/aarnavbhat23/protein-cutoff-audit/actions/workflows/tests.yml/badge.svg)](https://github.com/aarnavbhat23/protein-cutoff-audit/actions/workflows/tests.yml)

Temporal cutoffs are becoming the standard defence against contamination in scientific benchmarks: evaluate only on data "released after" a model's training date. For proteins, the date that gets used is UniProt's `date created`. **That date says when an entry was curated, not when the protein became known.**

This repository measures what that does to a benchmark built on it. [LiveProteinBench](https://arxiv.org/abs/2512.22257) describes itself as contamination-free because every protein was "first publicly released after January 1, 2025". We audit all **1,797 questions on 651 proteins**, pinned to commit [`8b97109`](https://github.com/Rongdingyi/LiveProteinBench/tree/8b9710932511f816651fc1a7e8a1bc85c3c77767), against UniProt's own archives.

![When the benchmark's proteins were first seen](results/figures/first_seen.png)

## What we found

| Check | Result |
|---|---:|
| Proteins whose first Swiss-Prot version is ≥ 2025 (the date the benchmark used) | **100%** |
| Sequence first seen anywhere (UniParc) before the cutoff | **95.4%** |
| … before 2015 / before 2010 | 68.5% / 50.7% (median: July 2009) |
| Same accession already public in UniProtKB/TrEMBL before the cutoff | **89.1%** |
| UniRef50 cluster contains a protein created before the cutoff | **92.9%** (median 22 such members) |
| Every PubMed paper cited by the current entry predates the cutoff | **93.2%** (median newest citation: Aug 2021) |
| Questions whose correct option, and no wrong option, is in the protein's **own pre-2025 entry** | **36.1%** (649 / 1,797) |
| Accuracy of a no-learning lookup in **Swiss-Prot 2024_06** (7 annotation tasks, 1,370 questions; chance 25%) | **70.7%** |
| … using only relatives below 50% identity | 66.8% |
| Questions on proteins whose sequence *and* 50% family are new at the cutoff | **4 (0.2%)** |

The date that defines the benchmark is when each protein reached Swiss-Prot curation. In 2025 that was typically 15–20 years after the sequence, its unreviewed entry, its annotations and its literature were public.

![Answer already in the protein's own pre-2025 entry](results/figures/own_entry_leakage.png)

**The answer was often already in the protein's own entry.** For each protein we fetch the last UniProtKB version released before 2025-01-01 from UniSave, UniProt's entry archive. For Q8MVC4, the November 2024 TrEMBL entry already carried the benchmark's correct answer (`GO:0004867 serine-type endopeptidase inhibitor activity`), and its function paper dates from 2004. Leakage is highest for GO molecular function (76%) and lowest for the experimental pH and temperature optima (0%), which TrEMBL does not carry.

![Lookup vs reported models](results/figures/lookup_vs_models.png)

**Pre-cutoff relatives answer most questions.** A lookup with no learning searches each protein against Swiss-Prot release 2024_06 (27 Nov 2024) with MMseqs2. It picks the option best supported by the annotations of relatives at ≥ 30% identity. It beats the best model in the paper on EC number (71.5% vs 31.0%), pathway (76.6% vs 30.5%) and cofactor (81.2% vs 59.1%). It roughly ties GPT-5 on the GO tasks, and trails it on catalytic activity, where our exact-reaction matching is strict. When a relative supports an option, that option is correct 97–100% of the time.

Annotation transfer from homologs is a legitimate method, not cheating. The point is narrower: **a date cutoff does not make the questions new.** The benchmark largely measures recall of family-level knowledge that predates its cutoff.

![UniProt-wide](results/figures/survey.png)

**It is not specific to this benchmark.** In samples of Swiss-Prot entries created in each year from 2016 to 2026, 93–100% have sequences more than a year older than their creation date. The median gap has grown from 13 years (2016) to 22 years (2026). Any "created after" filter on reviewed UniProt entries inherits this.

**Not yet measured: whether this changes protein-LM scores.** ESM-2 was trained on UniRef 2021_04, and 83.1% of the benchmark's sequences were first seen before then. `notebooks/plm_colab.ipynb` compares ESM-2 pseudo-perplexity for benchmark proteins first seen before vs after April 2021, and for the NovelFam set of genuinely new proteins (`data/novelfam/`). It is not run here: it needs a GPU, and the laptop used for this audit stalled on it. No result is claimed.

## A cleaner subset

`results/clean_subset.tsv` flags every question (by id) through four nested filters:

| Filter | Questions left |
|---|---:|
| All | 1,797 |
| Answer not supported by the protein's own pre-cutoff entry | 1,142 |
| … and no Swiss-Prot 2024_06 relative supports any option (annotation tasks) | 780 |
| … and the sequence was first seen after the cutoff | 28 |
| … and its UniRef50 cluster has no member created before the cutoff | 4 |

A benchmark that wants genuinely unseen proteins has to be built from new sequences and new families. As our [NovelFam](https://github.com/aarnavbhat23/next-token-bio) build shows, that pool is small: 64 of 10,127 recently created UniProt entries survived the same filters.

## Use it on your own data

```bash
printf "Q8MVC4\nP0DXN9\n" | python -m cutoffaudit.seqdate --cutoff 2025-01-01 -
```

| verdict | meaning |
|---|---|
| `old_sequence` | UniParc first saw the sequence before the cutoff |
| `new_sequence_old_family` | new sequence, but its UniRef50 cluster has a member created before the cutoff |
| `new` | sequence and its 50% neighbourhood are new |

## Method

- **Dates.** UniParc `first_seen` (when the sequence first appeared in any source database), UniSave entry histories (first UniProtKB version, first Swiss-Prot version, sequence-version changes), and a count of UniRef50 cluster members created before the cutoff.
- **Own-entry leakage.** One matcher per task compares each option with the specific annotation the question is about. The annotations are: GO terms of the right aspect, EC numbers, `COFACTOR` names, `CATALYTIC ACTIVITY` reactions, `PATHWAY` text or keywords, `ACT_SITE` / `MOTIF` / `TRANSMEM` locations, and the stated pH and temperature optima. A question counts as leaked only if the correct option is supported and no wrong option is.
- **Calibration.** Run on the *current* entries, the same matchers find the correct option alone for 86.9% of questions and flag a wrong option alone for 0.06%. So they are precise, and their recall bounds the leakage numbers from below.
- **Relatives baseline.** MMseqs2 `-s 7.5`, E < 1e-3, ≥ 30% identity over ≥ 50% of the query, top 10 relatives, options scored by the best supporting identity, ties split, no support = uniform guess. Positional and numeric tasks are excluded: active site, motif and transmembrane would need coordinates mapped through the alignment, and pH and temperature optima are rarely annotated.
- **Literature.** Publication dates of every PubMed reference in the current entry (NCBI E-utilities).

## Limits

- This shows that the information **was public** before the cutoff. It cannot show that any particular model **was trained on** it; model training sets are not disclosed.
- The matchers are exact. Paraphrased or renamed annotations are missed, so leakage is a lower bound (calibration recall: 69% for GO biological process, 72% for catalytic activity).
- The survey samples the first page of each month's results (≤ 50 entries), so the TrEMBL series in particular is noisy.
- UniParc dates the sequence, not the protein's characterisation. A new sequence can still belong to a well-studied family, which is why the cluster check is reported separately.

## Prior work and credit

- **LiveProteinBench:** Rong et al., [arXiv:2512.22257](https://arxiv.org/abs/2512.22257), [code and data](https://github.com/Rongdingyi/LiveProteinBench). This audit uses their released questions and reported accuracies (their Table 2) and does not redistribute their files.
- **Related concerns:** a [secondary summary](https://www.emergentmind.com/topics/contamination-free-rolling-benchmark) noted that a post-2025 timestamp "does not logically exclude earlier homologous information". Benchmark audits such as the [biomolecular benchmark audit](https://pubs.rsc.org/sc/article/doi/10.1039/d6sc01799a/1285850/Auditing-widely-used-biomolecular-benchmarks) and work on [leakage in protein–ligand benchmarks](https://www.biorxiv.org/content/10.64898/2026.06.29.735309v1) study split contamination.

  We found no prior quantitative audit of UniProt-date cutoffs using sequence first-seen dates and entry histories. If one exists, please open an issue.
- **Data:** [UniProt](https://www.uniprot.org) (UniProtKB, UniParc, UniRef, UniSave; CC BY 4.0), NCBI PubMed E-utilities, [MMseqs2](https://github.com/soedinglab/MMseqs2).

## Reproduce

Python 3.12+ with NumPy and matplotlib; `mmseqs` for the relatives baseline; `fair-esm` and PyTorch for the optional protein-LM check. All HTTP responses are cached under `data/cache/`.

```bash
python -m unittest discover -s tests -v
python -m cutoffaudit.audit all        # dates + own-entry leakage + calibration  -> results/proteins.tsv, items.tsv, summary.json
python -m cutoffaudit.literature       # -> results/literature.*
python -m cutoffaudit.survey           # -> results/survey.json
./scripts/fetch_sprot_2024_06.sh       # ~1.7 GB download, keeps uniprot_sprot.dat.gz (680 MB)
python -m cutoffaudit.homolog all      # index, MMseqs2 search, score, sensitivity
python -m cutoffaudit.subset
python -m cutoffaudit.plm              # optional; on a GPU (notebooks/plm_colab.ipynb)
python -m cutoffaudit.report           # figures + results/summary.md
```

Code: MIT. Built with AI assistance. Every number comes from the files in `results/`. There is no independent review or adoption claim.
