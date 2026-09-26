# Pooled CRISPRi screening — NGS deconvolution and library representation QC

Sequencing-based deconvolution and quality control of a genome-wide CRISPRi
screening library (Broad GPP Human Dolcetto, Addgene
[pooled library](https://www.addgene.org/pooled-library/broadgpp-human-crispri-dolcetto/),
`pXPR_050` backbone, used with lenti-KRAB-dCas9).

A pooled screen is only interpretable if every guide in the pool can be counted
accurately from sequencing reads. This repository contains the workflow that does
the counting — FASTQ in, per-guide and per-gene counts plus representation
statistics out — and applies it to both amplified library sets before they are
used in a screen.

Set A: 57,050 sgRNAs / 18,901 genes. Set B: 57,011 sgRNAs / 18,899 genes.
Three guides per gene per set, 500 non-targeting controls, ~148M reads total.

## What the workflow does

| Step | Script | Output |
| --- | --- | --- |
| 1. Deconvolute reads | `extract_grna.py` | 20-nt spacer counts from FASTQ, with per-read classification (missing 5' `CACCG`, missing 3' `GTTT`, truncated, below Phred 20, valid) |
| 2. Assign to library | `map_grna.py` | Counts mapped to the reference library, exact and 1-mismatch (unique vs ambiguous), plus non-reference reads |
| 3. Representation statistics | `map_grna.py` | Mean/median reads per guide, zero-count guides, 90th/10th percentile skew, Gini coefficient, log2 fold-change vs uniform |
| 4. Gene-level power | `gene_level_coverage.py`, `compare_thresholds.py` | Per gene, how many of its 3 guides clear a read threshold (300x / 200x / 100x) |
| 5. Plots and report | `visualize_qc.py`, `generate_report.py` | Count distribution, Lorenz curve, rank-abundance, gene coverage; markdown QC report |

Step 1 and 2 are the part that carries over to a screen readout: the same
anchor-based extraction and mismatch-tolerant assignment turns amplicon reads
from a selected cell pool into a count matrix. Hit calling from those counts
(log2 fold-change of treated vs plasmid/early time point, gene-level
aggregation) is a downstream step and is not included here.

Guide structure used for extraction:

```
5' ... CACCG [20-nt spacer] GTTT ... 3'
```

## Results

| Set | Reads | Exact-mapped | Guides detected | Mean reads/guide | Zero-count guides | Log-count Gini | Genes with >=2 guides at 300x | QC |
| :-- | --: | --: | --: | --: | --: | --: | --: | :-: |
| A | 80,411,815 | 89.29% | 57,023 / 57,050 (99.95%) | 1,259 | 27 (0.05%) | 0.052 | 96.90% | PASS |
| B | 67,346,721 | 90.39% | 56,996 / 57,011 (99.97%) | 1,068 | 15 (0.03%) | 0.058 | 94.73% | PASS |

Both sets clear the thresholds used here (>=65% exact-mapped, >=300 mean reads
per guide, <=1% zero-count guides, log-count Gini <=0.10, and >=90% of genes
with at least two guides above 300x). The raw-count skew ratio (~14 for both
sets) sits above the flag value of 10, which is typical for an amplified pooled
library and does not by itself compromise screening power — the gene-level
coverage numbers are the more informative measure.

Dropout at 300x is small and asymmetric: 48 genes in Set A and 103 in Set B have
no guide above threshold, with three genes (`ATP6V1E2`, `SETBP1`, `ZNF804A`)
dropping out in both. Relaxing the threshold to 100x removes nearly all of it
(>=2 guides for 99.58% of Set A genes and 99.24% of Set B genes), which is the
trade-off to weigh when setting screen sequencing depth.

Full numbers, per-gene tables and figures are in
[`2026_07_hpc_analysis/`](2026_07_hpc_analysis/README.md).

## Repository layout

```
2026_07_hpc_analysis/     Current pipeline (Python, SLURM)
  scripts/                Extraction, mapping, coverage, plotting
  crispri_analysis.slurm  Batch job for both sets
  results/ figures/       QC metrics, reports, plots
2025_06_old_analysis/     Earlier R implementation, kept for comparison
LIbrary_QC_Novogene/      Vendor sequencing reports (FASTQs not tracked)
Library target genes/     Dolcetto Set A / Set B reference files
```

## Running it

On a SLURM cluster, both sets end to end (~40 min, 4 CPUs, 64 GB):

```bash
sbatch 2026_07_hpc_analysis/crispri_analysis.slurm
```

One set locally:

```bash
cd 2026_07_hpc_analysis
./scripts/run_analysis.sh SetB \
  "../LIbrary_QC_Novogene/01.RawData/sg_set_B/sg_set_B_..._1.fq.gz" \
  "../Library target genes/broadgpp-dolcetto-targets-setb.txt"
```

Requires Python 3.7+ with pandas, numpy, matplotlib and tabulate. FASTQs are not
in the repository (~9 GB); paths are set at the top of the SLURM script.

## Notes

- Only R1 is used for counting; the spacer sits entirely within the first read.
- 1-mismatch matches are reported separately from exact matches so sequencing
  error is not silently folded into the counts, and ambiguous 1-mismatch reads
  (matching more than one reference guide) are tracked rather than assigned.
- The Lorenz curve and Gini coefficient measure how unevenly reads are spread
  across guides; log-count Gini is used for the pass/fail call because raw-count
  Gini is dominated by the long tail.
