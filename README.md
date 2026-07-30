# CRISPRi Library QC Analysis

Quality control analysis of the Broad GPP Dolcetto CRISPRi screening library (Addgene: [broadgpp-human-crispri-dolcetto](https://www.addgene.org/pooled-library/broadgpp-human-crispri-dolcetto/)).

The library consists of two sets (A and B), each with ~57,000 gRNAs targeting ~18,700 human genes. This repository contains the analysis pipeline, sequencing data, and QC results.

## Repository Structure

```
CRISPRi/
├── 2026_07_hpc_analysis/     # Current HPC pipeline (SLURM-based)
│   ├── scripts/              # Python analysis scripts
│   ├── results/              # QC metrics and reports
│   └── figures/              # QC plots (Set A & Set B)
├── 2025_06_old_analysis/     # Previous R-based analysis
├── LIbrary_QC_Novogene/      # Novogene sequencing reports
│   ├── 01.RawData/           # FASTQ files (not in git, ~9GB)
│   └── 02.Report_*/          # Sequencing QC reports
└── Library target genes/     # gRNA library reference files
```

## Quick Start

### Run on HPC (SLURM)

```bash
sbatch 2026_07_hpc_analysis/crispri_analysis.slurm
```

### Run Locally

```bash
cd 2026_07_hpc_analysis
./scripts/run_analysis.sh SetB \
  "../LIbrary_QC_Novogene/01.RawData/sg_set_B/sg_set_B_CKDL260013755-1A_23KG5MLT4_L3_1.fq.gz" \
  "../Library target genes/broadgpp-dolcetto-targets-setb.txt"
```

## Dependencies

- Python 3.7+
- pandas, numpy, matplotlib
- SLURM (for HPC batch jobs)

## Results Summary

| Set | Total Reads | % In Library | Coverage | Mean Reads/Guide | Status |
|:---:|------------:|-------------:|---------:|-----------------:|:------:|
| A | 80,411,815 | 89.29% | 99.95% | 1,259 | PASS |
| B | 67,346,721 | 90.39% | 99.97% | 1,068 | PASS |

Both sets exceed QC thresholds (>94% of genes have >=2 gRNAs at 300x coverage).

See [`2026_07_hpc_analysis/README.md`](2026_07_hpc_analysis/README.md) for detailed analysis documentation.
