#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=standard
#SBATCH --mem=64G
#SBATCH -J "crispri_library_qc"
#SBATCH --mail-user=subash@virginia.edu
#SBATCH --mail-type=ALL
#SBATCH --account=mohilab

# Enable strict error handling
set -euo pipefail

# Load required modules
module purge
module load miniforge

# --- Configuration ---
PROJECT_DIR="/sfs/weka/scratch/uwy2ak/CRISPRi"
SCRIPT_DIR="${PROJECT_DIR}/scripts"
RESULTS_DIR="${PROJECT_DIR}/results"
FIGURES_DIR="${PROJECT_DIR}/figures"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Use system Python - packages already installed
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"

which python3
python3 --version
python3 -c "import pandas; import numpy; import matplotlib; print('All required packages available')"

# Library reference files
SETA_LIBRARY="${PROJECT_DIR}/Library target genes/broadgpp-dolcetto-targets-seta.txt"
SETB_LIBRARY="${PROJECT_DIR}/Library target genes/broadgpp-dolcetto-targets-setb.txt"

# FASTQ files
SETA_FASTQ="${PROJECT_DIR}/LIbrary_QC_Novogene/01.RawData/sg_set_A/sg_set_A_CKDL260013754-1A_23KG5MLT4_L3_1.fq.gz"
SETB_FASTQ="${PROJECT_DIR}/LIbrary_QC_Novogene/01.RawData/sg_set_B/sg_set_B_CKDL260013755-1A_23KG5MLT4_L3_1.fq.gz"

echo "========================================"
echo "CRISPRi Library QC Analysis Pipeline"
echo "Start time: $(date)"
echo "========================================"
echo ""
echo "Configuration:"
echo "  Project Dir: ${PROJECT_DIR}"
echo "  Results Dir: ${RESULTS_DIR}"
echo "  Figures Dir: ${FIGURES_DIR}"
echo "  Threads: ${THREADS}"
echo ""

# Create output directories
mkdir -p "${RESULTS_DIR}/SetA" "${RESULTS_DIR}/SetA_lowqual"
mkdir -p "${RESULTS_DIR}/SetB" "${RESULTS_DIR}/SetB_lowqual"
mkdir -p "${FIGURES_DIR}/SetA" "${FIGURES_DIR}/SetB"

# Function to run gRNA extraction
run_extraction() {
    local set_name="$1"
    local fastq_file="$2"
    local output_dir="$3"
    local min_quality="$4"
    local counts_output="${output_dir}/${set_name}_grna_counts.txt"

    echo ""
    echo "========================================"
    echo "Step 1: Extracting gRNA from ${set_name}"
    echo "  FASTQ: ${fastq_file}"
    echo "  Output: ${counts_output}"
    echo "  Min Quality: ${min_quality}"
    echo "========================================"

    if [ ! -f "${fastq_file}" ]; then
        echo "ERROR: FASTQ file not found: ${fastq_file}"
        return 1
    fi

    python3 "${SCRIPT_DIR}/extract_grna.py" \
        "${fastq_file}" \
        -o "${counts_output}" \
        --min-quality "${min_quality}"

    if [ ! -f "${counts_output}" ]; then
        echo "ERROR: gRNA extraction failed for ${set_name}"
        return 1
    fi

    echo "Extraction complete for ${set_name}"
    ls -lh "${counts_output}"
}

# Function to run gRNA mapping
run_mapping() {
    local set_name="$1"
    local counts_file="$2"
    local library_file="$3"
    local output_dir="$4"
    local mapped_output="${output_dir}/${set_name}_mapped_counts.txt"
    local metrics_output="${output_dir}/${set_name}_metrics.txt"

    echo ""
    echo "========================================"
    echo "Step 2: Mapping ${set_name} to library"
    echo "  Counts: ${counts_file}"
    echo "  Library: ${library_file}"
    echo "  Output: ${mapped_output}"
    echo "========================================"

    if [ ! -f "${counts_file}" ]; then
        echo "ERROR: Counts file not found: ${counts_file}"
        return 1
    fi

    if [ ! -f "${library_file}" ]; then
        echo "ERROR: Library file not found: ${library_file}"
        return 1
    fi

    python3 "${SCRIPT_DIR}/map_grna.py" \
        "${counts_file}" \
        "${library_file}" \
        -o "${mapped_output}" \
        -m "${metrics_output}"

    if [ ! -f "${metrics_output}" ]; then
        echo "ERROR: Mapping failed for ${set_name}"
        return 1
    fi

    echo "Mapping complete for ${set_name}"
    echo "Metrics:"
    cat "${metrics_output}"
}

# Function to run visualization
run_visualization() {
    local set_name="$1"
    local mapped_file="$2"
    local fig_dir="$3"

    echo ""
    echo "========================================"
    echo "Step 3: Generating visualizations for ${set_name}"
    echo "  Input: ${mapped_file}"
    echo "  Output: ${fig_dir}"
    echo "========================================"

    if [ ! -f "${mapped_file}" ]; then
        echo "ERROR: Mapped counts not found: ${mapped_file}"
        return 1
    fi

    python3 "${SCRIPT_DIR}/visualize_qc.py" \
        "${mapped_file}" \
        -o "${fig_dir}" \
        -n "${set_name}"

    echo "Visualizations complete for ${set_name}"
    ls -lh "${fig_dir}/${set_name}"*.png 2>/dev/null || echo "No PNG files found"
}

# Function to run complete analysis for one set
run_full_analysis() {
    local set_name="$1"
    local fastq_file="$2"
    local library_file="$3"
    local results_subdir="$4"
    local fig_subdir="$5"
    local min_quality="${6:-20}"

    echo ""
    echo "############################################################"
    echo "# Running Full Analysis: ${set_name}"
    echo "############################################################"

    run_extraction "${set_name}" "${fastq_file}" "${results_subdir}" "${min_quality}"
    run_mapping "${set_name}" "${results_subdir}/${set_name}_grna_counts.txt" "${library_file}" "${results_subdir}"
    run_visualization "${set_name}" "${results_subdir}/${set_name}_mapped_counts.txt" "${fig_subdir}"

    echo ""
    echo "############################################################"
    echo "# ${set_name} Analysis COMPLETE"
    echo "############################################################"
}

# Main execution
echo "Starting CRISPRi Library QC Analysis..."
echo "Python: $(which python3)"
echo "Python version: $(python3 --version)"

# Check scripts exist
for script in extract_grna.py map_grna.py visualize_qc.py; do
    if [ ! -f "${SCRIPT_DIR}/${script}" ]; then
        echo "ERROR: Required script not found: ${SCRIPT_DIR}/${script}"
        exit 1
    fi
done

# Run Set A analysis (high quality filter)
if [ -f "${SETA_FASTQ}" ]; then
    run_full_analysis "SetA" "${SETA_FASTQ}" "${SETA_LIBRARY}" "${RESULTS_DIR}/SetA" "${FIGURES_DIR}/SetA" 20
else
    echo "WARNING: Set A FASTQ not found, skipping: ${SETA_FASTQ}"
fi

# Run Set B analysis (high quality filter)
if [ -f "${SETB_FASTQ}" ]; then
    run_full_analysis "SetB" "${SETB_FASTQ}" "${SETB_LIBRARY}" "${RESULTS_DIR}/SetB" "${FIGURES_DIR}/SetB" 20
else
    echo "WARNING: Set B FASTQ not found, skipping: ${SETB_FASTQ}"
fi

# Generate summary report
echo ""
echo "========================================"
echo "Generating Summary Report"
echo "========================================"

REPORT_FILE="${RESULTS_DIR}/analysis_summary.txt"
{
    echo "CRISPRi Library QC Analysis Summary"
    echo "==================================="
    echo "Generated: $(date)"
    echo ""

    for set in SetA SetB; do
        echo "${set} Results:"
        echo "---------------"
        METRICS="${RESULTS_DIR}/${set}/${set}_metrics.txt"
        if [ -f "${METRICS}" ]; then
            cat "${METRICS}"
        else
            echo "Metrics not found"
        fi
        echo ""
    done

    echo ""
    echo "Output Files:"
    echo "-------------"
    find "${RESULTS_DIR}" -name "*.txt" -type f | sort
    echo ""
    echo "Figures:"
    echo "--------"
    find "${FIGURES_DIR}" -name "*.png" -type f | sort

} > "${REPORT_FILE}"

echo ""
echo "========================================"
echo "CRISPRi Library QC Analysis COMPLETE!"
echo "========================================"
echo ""
echo "End time: $(date)"
echo ""
echo "Results Summary:"
echo "  Set A Results: ${RESULTS_DIR}/SetA/"
echo "  Set B Results: ${RESULTS_DIR}/SetB/"
echo "  Set A Figures: ${FIGURES_DIR}/SetA/"
echo "  Set B Figures: ${FIGURES_DIR}/SetB/"
echo "  Summary Report: ${REPORT_FILE}"
echo ""
echo "To check results:"
echo "  cat ${REPORT_FILE}"
echo ""
