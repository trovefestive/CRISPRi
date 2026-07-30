#!/bin/bash
# Master pipeline for CRISPRi library QC analysis

set -e

# Input parameters
SET_NAME=$1
FASTQ_FILE=$2
LIBRARY_FILE=$3

if [ -z "$SET_NAME" ] || [ -z "$FASTQ_FILE" ] || [ -z "$LIBRARY_FILE" ]; then
    echo "Usage: ./run_analysis.sh <set_name> <fastq_file> <library_file>"
    echo "Example: ./run_analysis.sh SetA path/to/setA.fq.gz 'Library target genes/broadgpp-dolcetto-targets-seta.txt'"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="results/${SET_NAME}"
FIGURES_DIR="figures/${SET_NAME}"

mkdir -p "$RESULTS_DIR" "$FIGURES_DIR"

echo "=========================================="
echo "CRISPRi Library QC Analysis: $SET_NAME"
echo "=========================================="
echo ""

# Step 1: Extract gRNA counts
echo "Step 1: Extracting gRNA counts from FASTQ..."
python3 "$SCRIPT_DIR/extract_grna.py" \
    "$FASTQ_FILE" \
    -o "$RESULTS_DIR/grna_counts.txt" \
    --min-quality 20

echo ""
echo "Step 2: Mapping to library reference..."
python3 "$SCRIPT_DIR/map_grna.py" \
    "$RESULTS_DIR/grna_counts.txt" \
    "$LIBRARY_FILE" \
    -o "$RESULTS_DIR/mapped_counts.txt" \
    -m "$RESULTS_DIR/metrics.txt"

echo ""
echo "Step 3: Generating visualizations..."
python3 "$SCRIPT_DIR/visualize_qc.py" \
    "$RESULTS_DIR/mapped_counts.txt" \
    -o "$FIGURES_DIR" \
    -n "$SET_NAME"

echo ""
echo "=========================================="
echo "Analysis complete for $SET_NAME!"
echo "=========================================="
echo "Results: $RESULTS_DIR/"
echo "Figures: $FIGURES_DIR/"
