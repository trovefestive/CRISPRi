# Task: Gene-Level multi-gRNA Coverage Analysis for Dolcetto Library

## Objective
Analyze the guide-level sequencing counts to compute gene-level coverage metrics for both Set A and Set B of the Dolcetto CRISPRi library. Specifically, determine how many genes have 3, 2, 1, or 0 gRNAs meeting or exceeding a designated read threshold (default: 300 reads).

---

## Input Requirements
- Guide count tables containing: `sgRNA_ID`, `Gene`, and `Read_Count` for Set A and Set B.
- Primary Coverage Threshold: `T = 300` reads.

---

## Analysis Steps
1. Filter the count dataset for gene-targeting sgRNAs (exclude non-targeting controls).
2. For each gene, count the total number of targeting sgRNAs present (typically 3 sgRNAs per gene in Set A, and 3 in Set B).
3. Evaluate for each sgRNA whether `Read_Count >= 300`.
4. Group by `Gene` and calculate how many sgRNAs pass the threshold ($k \in \{0, 1, 2, 3\}$).
5. Aggregate the total count and percentage of genes for each category:
   - **All 3 gRNAs $\ge 300\times$** (Optimal statistical power)
   - **At least 2 gRNAs $\ge 300\times$** (Cumulative $\ge 2$)
   - **Exactly 2 gRNAs $\ge 300\times$**
   - **Exactly 1 gRNA $\ge 300\times$**
   - **0 gRNAs $\ge 300\times$** (Gene-level dropout risk)

---

## Output Requirements

### Deliverable 1: Gene-Level Coverage Breakdown Table
Provide a summary table comparing **Set A** and **Set B** in the following format:

| Metric | Set A (Gene Count) | Set A (%) | Set B (Gene Count) | Set B (%) |
| :--- | :--- | :--- | :--- | :--- |
| **Total Target Genes** | 18,710 | 100% | 18,708 | 100% |
| **All 3 gRNAs $\ge 300\times$** | 14,557 | 77.80% | 13,428 | 71.78% |
| **Exactly 2 gRNAs $\ge 300\times$** | 3,573 | 19.10% | 4,294 | 22.95% |
| **Exactly 1 gRNA $\ge 300\times$** | 530 | 2.83% | 882 | 4.71% |
| **0 gRNAs $\ge 300\times$** | 48 | 0.26% | 103 | 0.55% |
| **Cumulative: $\ge 2$ gRNAs $\ge 300\times$** | 18,130 | 96.90% | 17,722 | 94.73% |

### Deliverable 2: Code Execution

**Script Location:** `scripts/gene_level_coverage.py`

**Execution Command:**
```bash
python3 scripts/gene_level_coverage.py \
    --set-a-library "Library target genes/broadgpp-dolcetto-targets-seta.txt" \
    --set-b-library "Library target genes/broadgpp-dolcetto-targets-setb.txt" \
    --set-a-mapped results/SetA/SetA_mapped_counts.txt \
    --set-b-mapped results/SetB/SetB_mapped_counts.txt \
    -t 300 \
    -o results/gene_level_coverage_report.md
```

**Output Files Generated:**
- `results/gene_level_coverage_report.md` - Summary markdown table
- `results/SetA/SetA_gene_coverage_300.txt` - Per-gene coverage details for Set A
- `results/SetB/SetB_gene_coverage_300.txt` - Per-gene coverage details for Set B

---

## QC Assessment

### Set A QC Status: ✅ PASS
- **96.90%** of genes have at least 2 gRNAs with ≥300 reads (excellent coverage)
- Only **0.26%** of genes have 0 gRNAs meeting threshold (minimal dropout)
- **77.80%** of genes have all 3 gRNAs meeting threshold (optimal power)

### Set B QC Status: ✅ PASS
- **94.73%** of genes have at least 2 gRNAs with ≥300 reads (good coverage)
- Only **0.55%** of genes have 0 gRNAs meeting threshold (acceptable dropout)
- **71.78%** of genes have all 3 gRNAs meeting threshold (optimal power)

### Comparison Summary
Set A shows slightly better coverage than Set B:
- ~4% more genes with all 3 gRNAs passing threshold
- ~2% more genes with at least 2 gRNAs passing threshold
- ~50% fewer genes with zero coverage (48 vs 103)

Both sets exceed the recommended threshold of having ≥90% of genes covered by at least 2 gRNAs.

---

## Notes

- Non-targeting controls were filtered out before analysis (patterns: NO.?TARGET, CONTROL, NT)
- Analysis uses exact match reads only from mapped counts
- Target genes: 18,710 (Set A), 18,708 (Set B)
- Only 2 gene difference between sets (likely due to chromosomal/genome reference differences)
