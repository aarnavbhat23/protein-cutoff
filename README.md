# Protein Cutoff Audit

Of the 651 proteins in LiveProteinBench, 95.4% had sequences in a public archive before the cutoff. Swiss-Prot’s curation date records when the entry was curated, which can be years after the sequence first appeared.

![Earliest public archive dates for the benchmark proteins](results/figures/first_seen.png)

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

I started by tracing when the sequences appeared. Then I checked whether the old entries already gave the answer. For 649 of 1,797 questions, the entry supported the correct option and none of the alternatives.

## Limits

Archive dates show what was public by the cutoff. The audit has no access to individual training sets, and exact annotation matching can miss paraphrases.

[Methods, results and reproduction commands](docs/reproduce.md)

Built with AI assistance.
