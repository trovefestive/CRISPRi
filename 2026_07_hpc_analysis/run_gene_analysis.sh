#!/bin/bash
cd /sfs/weka/scratch/uwy2ak/CRISPRi
python3 scripts/gene_level_coverage.py \
    --set-a-library "Library target genes/broadgpp-dolcetto-targets-seta.txt" \
    --set-b-library "Library target genes/broadgpp-dolcetto-targets-setb.txt" \
    --set-a-mapped results/SetA/SetA_mapped_counts.txt \
    --set-b-mapped results/SetB/SetB_mapped_counts.txt \
    -t 300 \
    -o results/gene_level_coverage_report.md
