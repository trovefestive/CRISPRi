#!/usr/bin/env python3
"""
Compare gene-level coverage across multiple thresholds (100x, 200x, 300x).

Generates a comprehensive table comparing ≥100×, ≥200×, and ≥300× coverage
for each gene across Set A and Set B.
"""

import argparse
import pandas as pd
from pathlib import Path
from collections import Counter


def load_library(library_path):
    """Load library reference file."""
    df = pd.read_csv(library_path, sep='\t')
    df.columns = ['grna_sequence', 'gene_symbol', 'gene_id']
    df['grna_sequence'] = df['grna_sequence'].str.upper()
    return df


def load_mapped_counts(mapped_path):
    """Load mapped counts from map_grna.py output."""
    df = pd.read_csv(mapped_path, sep='\t')
    df['grna_sequence'] = df['grna_sequence'].str.upper()
    return df


def analyze_gene_coverage(library_df, mapped_df, threshold):
    """
    Analyze gene-level coverage based on read count threshold.
    Returns dict with coverage statistics.
    """
    # Filter for gene-targeting sgRNAs (exclude non-targeting controls)
    gene_targeting = library_df[
        ~library_df['gene_symbol'].str.contains('NO.?TARGET|CONTROL|NT', case=False, na=False)
    ].copy()

    # Create lookup for counts from mapped dataframe
    count_lookup = dict(zip(mapped_df['grna_sequence'], mapped_df['count']))

    # Add count to library dataframe
    gene_targeting['count'] = gene_targeting['grna_sequence'].map(count_lookup).fillna(0).astype(int)

    # Evaluate if each gRNA meets threshold
    gene_targeting['meets_threshold'] = gene_targeting['count'] >= threshold

    # Group by gene and count how many gRNAs pass threshold
    gene_coverage = gene_targeting.groupby('gene_symbol').agg({
        'meets_threshold': ['sum', 'count']
    }).reset_index()

    gene_coverage.columns = ['gene_symbol', 'passing_guides', 'total_guides']
    gene_coverage['passing_guides'] = gene_coverage['passing_guides'].astype(int)

    # Count genes in each category
    total_genes = len(gene_coverage)
    counts = Counter(gene_coverage['passing_guides'])

    stats = {
        'total_genes': total_genes,
        'all_3': counts.get(3, 0),
        'exactly_2': counts.get(2, 0),
        'exactly_1': counts.get(1, 0),
        'zero': counts.get(0, 0),
        'cumulative_2_or_more': counts.get(2, 0) + counts.get(3, 0),
    }

    stats['all_3_pct'] = 100 * stats['all_3'] / total_genes if total_genes > 0 else 0
    stats['exactly_2_pct'] = 100 * stats['exactly_2'] / total_genes if total_genes > 0 else 0
    stats['exactly_1_pct'] = 100 * stats['exactly_1'] / total_genes if total_genes > 0 else 0
    stats['zero_pct'] = 100 * stats['zero'] / total_genes if total_genes > 0 else 0
    stats['cumulative_2_or_more_pct'] = 100 * stats['cumulative_2_or_more'] / total_genes if total_genes > 0 else 0

    return stats


def generate_comparison_table(set_a_stats_by_threshold, set_b_stats_by_threshold):
    """Generate markdown table comparing all thresholds."""
    lines = []
    lines.append("## Gene-Level Coverage Analysis: Multi-Threshold Comparison\n")

    # Set A Table
    lines.append("### Set A Coverage by Threshold\n")
    lines.append("| Metric | ≥300× (Count) | ≥300× (%) | ≥200× (Count) | ≥200× (%) | ≥100× (Count) | ≥100× (%) |")
    lines.append("| :--- | :--- | :--- | :--- | :--- | :--- | :--- |")

    rows = [
        ('Total Target Genes', 'total_genes', 'total_genes', False),
        ('All 3 gRNAs above threshold', 'all_3', 'all_3_pct', True),
        ('Exactly 2 gRNAs above threshold', 'exactly_2', 'exactly_2_pct', True),
        ('Exactly 1 gRNA above threshold', 'exactly_1', 'exactly_1_pct', True),
        ('0 gRNAs above threshold', 'zero', 'zero_pct', True),
        ('Cumulative: ≥2 gRNAs above threshold', 'cumulative_2_or_more', 'cumulative_2_or_more_pct', True),
    ]

    for label, count_key, pct_key, is_pct in rows:
        values = []
        for thresh in [300, 200, 100]:
            stats = set_a_stats_by_threshold[thresh]
            count = stats.get(count_key, 0)
            if is_pct:
                pct = f"{stats.get(pct_key, 0):.2f}%"
                values.append(f"{count:,} | {pct}")
            else:
                values.append(f"{count:,} | 100%")
        lines.append(f"| {label} | {values[0]} | {values[1]} | {values[2]} |")

    lines.append("")

    # Set B Table
    lines.append("### Set B Coverage by Threshold\n")
    lines.append("| Metric | ≥300× (Count) | ≥300× (%) | ≥200× (Count) | ≥200× (%) | ≥100× (Count) | ≥100× (%) |")
    lines.append("| :--- | :--- | :--- | :--- | :--- | :--- | :--- |")

    for label, count_key, pct_key, is_pct in rows:
        values = []
        for thresh in [300, 200, 100]:
            stats = set_b_stats_by_threshold[thresh]
            count = stats.get(count_key, 0)
            if is_pct:
                pct = f"{stats.get(pct_key, 0):.2f}%"
                values.append(f"{count:,} | {pct}")
            else:
                values.append(f"{count:,} | 100%")
        lines.append(f"| {label} | {values[0]} | {values[1]} | {values[2]} |")

    lines.append("")

    # Summary Comparison
    lines.append("### Key Metrics Comparison\n")
    lines.append("| Set | Metric | ≥300× | ≥200× | ≥100× | Improvement (300→100) |")
    lines.append("| :--- | :--- | :--- | :--- | :--- | :--- |")

    for set_name, stats_by_thresh in [('Set A', set_a_stats_by_threshold), ('Set B', set_b_stats_by_threshold)]:
        # All 3 gRNAs
        pct_300 = stats_by_thresh[300]['all_3_pct']
        pct_200 = stats_by_thresh[200]['all_3_pct']
        pct_100 = stats_by_thresh[100]['all_3_pct']
        improvement = pct_100 - pct_300
        lines.append(f"| {set_name} | All 3 gRNAs above threshold | {pct_300:.2f}% | {pct_200:.2f}% | {pct_100:.2f}% | +{improvement:.2f}% |")

        # ≥2 gRNAs
        pct_300 = stats_by_thresh[300]['cumulative_2_or_more_pct']
        pct_200 = stats_by_thresh[200]['cumulative_2_or_more_pct']
        pct_100 = stats_by_thresh[100]['cumulative_2_or_more_pct']
        improvement = pct_100 - pct_300
        lines.append(f"| {set_name} | ≥2 gRNAs above threshold | {pct_300:.2f}% | {pct_200:.2f}% | {pct_100:.2f}% | +{improvement:.2f}% |")

        # 0 gRNAs (dropout)
        pct_300 = stats_by_thresh[300]['zero_pct']
        pct_200 = stats_by_thresh[200]['zero_pct']
        pct_100 = stats_by_thresh[100]['zero_pct']
        improvement = pct_300 - pct_100
        lines.append(f"| {set_name} | 0 gRNAs above threshold (dropout) | {pct_300:.2f}% | {pct_200:.2f}% | {pct_100:.2f}% | -{improvement:.2f}% |")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description='Compare gene-level coverage across multiple thresholds')
    parser.add_argument('--set-a-library', required=True, help='Set A library reference file')
    parser.add_argument('--set-b-library', required=True, help='Set B library reference file')
    parser.add_argument('--set-a-mapped', required=True, help='Set A mapped counts file')
    parser.add_argument('--set-b-mapped', required=True, help='Set B mapped counts file')
    parser.add_argument('-o', '--output', default='results/multi_threshold_comparison.md',
                        help='Output markdown report file')
    args = parser.parse_args()

    # Load data
    print("Loading Set A data...")
    set_a_lib = load_library(args.set_a_library)
    set_a_mapped = load_mapped_counts(args.set_a_mapped)
    print(f"  Library size: {len(set_a_lib):,} guides")
    print(f"  Mapped reads: {len(set_a_mapped):,} guides with counts")

    print("\nLoading Set B data...")
    set_b_lib = load_library(args.set_b_library)
    set_b_mapped = load_mapped_counts(args.set_b_mapped)
    print(f"  Library size: {len(set_b_lib):,} guides")
    print(f"  Mapped reads: {len(set_b_mapped):,} guides with counts")

    # Analyze for each threshold
    thresholds = [300, 200, 100]
    set_a_stats_by_threshold = {}
    set_b_stats_by_threshold = {}

    for thresh in thresholds:
        print(f"\nAnalyzing coverage at {thresh}× threshold...")
        set_a_stats_by_threshold[thresh] = analyze_gene_coverage(set_a_lib, set_a_mapped, thresh)
        set_b_stats_by_threshold[thresh] = analyze_gene_coverage(set_b_lib, set_b_mapped, thresh)

    # Generate comparison report
    md_content = generate_comparison_table(set_a_stats_by_threshold, set_b_stats_by_threshold)

    # Write report
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(md_content)
        f.write("\n")

    print(f"\nMulti-threshold comparison report written to {output_path}")

    # Print to console
    print("\n" + "="*80)
    print("MULTI-THRESHOLD GENE-LEVEL COVERAGE COMPARISON")
    print("="*80)
    print(md_content)


if __name__ == '__main__':
    main()
