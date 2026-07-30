#!/usr/bin/env python3
"""
Gene-Level multi-gRNA Coverage Analysis for Dolcetto Library.

Analyzes guide-level counts to compute gene-level coverage metrics.
Determines how many genes have 3, 2, 1, or 0 gRNAs meeting a designated read threshold.

Default Coverage Threshold: 300 reads
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


def analyze_gene_coverage(library_df, mapped_df, threshold=300):
    """
    Analyze gene-level coverage based on read count threshold.

    Returns dict with coverage statistics.
    """
    # Filter for gene-targeting sgRNAs (exclude non-targeting controls)
    # Non-targeting controls typically have names like 'NO-TARGET', 'NON-TARGETING', etc.
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
        'meets_threshold': ['sum', 'count']  # sum = passing guides, count = total guides
    }).reset_index()

    gene_coverage.columns = ['gene_symbol', 'passing_guides', 'total_guides']
    gene_coverage['passing_guides'] = gene_coverage['passing_guides'].astype(int)

    # Count genes in each category
    total_genes = len(gene_coverage)

    counts = Counter(gene_coverage['passing_guides'])

    # Calculate statistics
    stats = {
        'total_genes': total_genes,
        'all_3': counts.get(3, 0),
        'exactly_2': counts.get(2, 0),
        'exactly_1': counts.get(1, 0),
        'zero': counts.get(0, 0),
        'cumulative_2_or_more': counts.get(2, 0) + counts.get(3, 0),
    }

    # Calculate percentages
    stats['all_3_pct'] = 100 * stats['all_3'] / total_genes if total_genes > 0 else 0
    stats['exactly_2_pct'] = 100 * stats['exactly_2'] / total_genes if total_genes > 0 else 0
    stats['exactly_1_pct'] = 100 * stats['exactly_1'] / total_genes if total_genes > 0 else 0
    stats['zero_pct'] = 100 * stats['zero'] / total_genes if total_genes > 0 else 0
    stats['cumulative_2_or_more_pct'] = 100 * stats['cumulative_2_or_more'] / total_genes if total_genes > 0 else 0

    return stats, gene_coverage


def print_coverage_table(set_a_stats, set_b_stats, threshold):
    """Print gene-level coverage comparison table."""
    print("\n" + "="*80)
    print(f"GENE-LEVEL MULTI-gRNA COVERAGE ANALYSIS (Threshold: {threshold} reads)")
    print("="*80)

    # Header
    print(f"{'Metric':<35} {'Set A (Count)':>12} {'Set A (%)':>10} {'Set B (Count)':>12} {'Set B (%)':>10}")
    print("-"*80)

    rows = [
        ('Total Target Genes', 'total_genes', 'total_genes', False),
        ('All 3 gRNAs >= {}x'.format(threshold), 'all_3', 'all_3_pct', True),
        ('Exactly 2 gRNAs >= {}x'.format(threshold), 'exactly_2', 'exactly_2_pct', True),
        ('Exactly 1 gRNA >= {}x'.format(threshold), 'exactly_1', 'exactly_1_pct', True),
        ('0 gRNAs >= {}x'.format(threshold), 'zero', 'zero_pct', True),
        ('Cumulative: >= 2 gRNAs >= {}x'.format(threshold), 'cumulative_2_or_more', 'cumulative_2_or_more_pct', True),
    ]

    for label, count_key, pct_key, is_pct in rows:
        a_count = set_a_stats.get(count_key, 0)
        b_count = set_b_stats.get(count_key, 0)

        if is_pct:
            a_pct = set_a_stats.get(pct_key, 0)
            b_pct = set_b_stats.get(pct_key, 0)
            print(f"{label:<35} {a_count:>12,} {a_pct:>9.2f}% {b_count:>12,} {b_pct:>9.2f}%")
        else:
            print(f"{label:<35} {a_count:>12,} {'100.00%':>10} {b_count:>12,} {'100.00%':>10}")

    print("="*80)


def generate_markdown_table(set_a_stats, set_b_stats, threshold):
    """Generate markdown format table for reports."""
    lines = []
    lines.append(f"## Gene-Level Coverage Analysis (Threshold: {threshold} reads)\n")
    lines.append("| Metric | Set A (Gene Count) | Set A (%) | Set B (Gene Count) | Set B (%) |")
    lines.append("| :--- | :--- | :--- | :--- | :--- |")

    rows = [
        ('Total Target Genes', 'total_genes', 'total_genes', False),
        ('All 3 gRNAs >= {}x'.format(threshold), 'all_3', 'all_3_pct', True),
        ('Exactly 2 gRNAs >= {}x'.format(threshold), 'exactly_2', 'exactly_2_pct', True),
        ('Exactly 1 gRNA >= {}x'.format(threshold), 'exactly_1', 'exactly_1_pct', True),
        ('0 gRNAs >= {}x'.format(threshold), 'zero', 'zero_pct', True),
        ('Cumulative: >= 2 gRNAs >= {}x'.format(threshold), 'cumulative_2_or_more', 'cumulative_2_or_more_pct', True),
    ]

    for label, count_key, pct_key, is_pct in rows:
        a_count = set_a_stats.get(count_key, 0)
        b_count = set_b_stats.get(count_key, 0)

        if is_pct:
            a_pct = f"{set_a_stats.get(pct_key, 0):.2f}%"
            b_pct = f"{set_b_stats.get(pct_key, 0):.2f}%"
            lines.append(f"| {label} | {a_count:,} | {a_pct} | {b_count:,} | {b_pct} |")
        else:
            lines.append(f"| {label} | {a_count:,} | 100% | {b_count:,} | 100% |")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description='Gene-level multi-gRNA coverage analysis')
    parser.add_argument('--set-a-library', required=True, help='Set A library reference file')
    parser.add_argument('--set-b-library', required=True, help='Set B library reference file')
    parser.add_argument('--set-a-mapped', required=True, help='Set A mapped counts file')
    parser.add_argument('--set-b-mapped', required=True, help='Set B mapped counts file')
    parser.add_argument('-t', '--threshold', type=int, default=300,
                        help='Coverage threshold in reads (default: 300)')
    parser.add_argument('-o', '--output', help='Output markdown report file')
    args = parser.parse_args()

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

    print(f"\nAnalyzing Set A gene coverage (threshold: {args.threshold} reads)...")
    set_a_stats, set_a_coverage = analyze_gene_coverage(set_a_lib, set_a_mapped, args.threshold)

    print(f"Analyzing Set B gene coverage (threshold: {args.threshold} reads)...")
    set_b_stats, set_b_coverage = analyze_gene_coverage(set_b_lib, set_b_mapped, args.threshold)

    # Print results
    print_coverage_table(set_a_stats, set_b_stats, args.threshold)

    # Generate markdown report if requested
    if args.output:
        md_table = generate_markdown_table(set_a_stats, set_b_stats, args.threshold)
        with open(args.output, 'w') as f:
            f.write(md_table)
            f.write("\n")
        print(f"\nMarkdown report written to {args.output}")

    # Save detailed coverage files
    base_a = Path(args.set_a_mapped).stem.replace('_mapped_counts', '')
    base_b = Path(args.set_b_mapped).stem.replace('_mapped_counts', '')

    output_a = Path(args.set_a_mapped).parent / f"{base_a}_gene_coverage_{args.threshold}.txt"
    output_b = Path(args.set_b_mapped).parent / f"{base_b}_gene_coverage_{args.threshold}.txt"

    set_a_coverage.to_csv(output_a, sep='\t', index=False)
    set_b_coverage.to_csv(output_b, sep='\t', index=False)

    print(f"\nDetailed gene coverage written to:")
    print(f"  {output_a}")
    print(f"  {output_b}")


if __name__ == '__main__':
    main()
