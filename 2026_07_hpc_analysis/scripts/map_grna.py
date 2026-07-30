#!/usr/bin/env python3
"""
Map extracted gRNA counts to reference library.
Calculate library representation metrics including 1-mismatch tolerant matching.
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict


def hamming_distance(s1, s2):
    """Calculate Hamming distance between two sequences."""
    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def load_library(library_path):
    """Load library reference file."""
    df = pd.read_csv(library_path, sep='\t')
    # Library has columns: Barcode Sequence, Annotated Gene Symbol, Annotated Gene ID
    df.columns = ['grna_sequence', 'gene_symbol', 'gene_id']
    df['grna_sequence'] = df['grna_sequence'].str.upper()
    return df


def load_counts(counts_path):
    """Load gRNA counts from extract_grna.py output."""
    df = pd.read_csv(counts_path, sep='\t')
    df['grna_sequence'] = df['grna_sequence'].str.upper()
    return df


def build_mismatch_map(library_seqs):
    """
    Build a lookup table for 1-mismatch sequences.
    Returns a dict mapping 1-mismatch sequences to list of matching reference sequences.
    """
    mismatch_map = defaultdict(list)

    for ref_seq in library_seqs:
        ref_seq = ref_seq.upper()
        # For each position, create variants with each other base
        bases = ['A', 'C', 'G', 'T']
        for i in range(len(ref_seq)):
            original_base = ref_seq[i]
            for base in bases:
                if base != original_base:
                    variant = ref_seq[:i] + base + ref_seq[i+1:]
                    mismatch_map[variant].append(ref_seq)

    return mismatch_map


def map_grnas_with_mismatches(counts_df, library_df):
    """
    Map counts to library reference with 1-mismatch tolerance.

    Returns:
        - mapped_df: DataFrame with mapping results
        - mapping_stats: dict with classification counts
    """
    # Create library lookup
    library_map = dict(zip(library_df['grna_sequence'],
                          zip(library_df['gene_symbol'], library_df['gene_id'])))
    library_seqs = list(library_df['grna_sequence'].values)

    # Mapping statistics
    mapping_stats = {
        'total_reads': 0,
        'exact_reference_match': 0,
        '1_mismatch_unique': 0,
        '1_mismatch_ambiguous': 0,
        'non_reference': 0,
    }

    # Pre-compute 1-mismatch mappings
    print("Building 1-mismatch lookup table...")
    mismatch_map = build_mismatch_map(library_seqs)
    print(f"  Lookup table size: {len(mismatch_map):,} entries")

    results = []

    for _, row in counts_df.iterrows():
        grna = row['grna_sequence']
        count = row['count']
        mapping_stats['total_reads'] += count

        if grna in library_map:
            # Exact match
            gene_symbol, gene_id = library_map[grna]
            results.append({
                'grna_sequence': grna,
                'gene_symbol': gene_symbol,
                'gene_id': gene_id,
                'count': count,
                'match_type': 'exact'
            })
            mapping_stats['exact_reference_match'] += count
        elif grna in mismatch_map:
            # 1-mismatch match
            matching_refs = mismatch_map[grna]
            if len(matching_refs) == 1:
                # Unique 1-mismatch
                ref_seq = matching_refs[0]
                gene_symbol, gene_id = library_map[ref_seq]
                results.append({
                    'grna_sequence': grna,
                    'gene_symbol': gene_symbol,
                    'gene_id': gene_id,
                    'count': count,
                    'match_type': '1_mismatch_unique'
                })
                mapping_stats['1_mismatch_unique'] += count
            else:
                # Ambiguous 1-mismatch
                results.append({
                    'grna_sequence': grna,
                    'gene_symbol': 'AMBIGUOUS_1MM',
                    'gene_id': 'AMBIGUOUS_1MM',
                    'count': count,
                    'match_type': '1_mismatch_ambiguous'
                })
                mapping_stats['1_mismatch_ambiguous'] += count
        else:
            # Non-reference
            results.append({
                'grna_sequence': grna,
                'gene_symbol': 'NON_REFERENCE',
                'gene_id': 'NON_REFERENCE',
                'count': count,
                'match_type': 'non_reference'
            })
            mapping_stats['non_reference'] += count

    return pd.DataFrame(results), mapping_stats


def calculate_metrics(mapped_df, library_df, mapping_stats):
    """Calculate library representation metrics."""
    metrics = {}

    # Basic counts from mapping stats
    total_reads = mapping_stats['total_reads']
    library_size = len(library_df)

    metrics['total_reads'] = total_reads
    metrics['library_size'] = library_size
    metrics['exact_reference_match'] = mapping_stats['exact_reference_match']
    metrics['1_mismatch_unique'] = mapping_stats['1_mismatch_unique']
    metrics['1_mismatch_ambiguous'] = mapping_stats['1_mismatch_ambiguous']
    metrics['non_reference'] = mapping_stats['non_reference']

    # Exact Mapped (%) - reads mapping exactly / total reads
    metrics['exact_mapped_pct'] = 100 * mapping_stats['exact_reference_match'] / total_reads if total_reads > 0 else 0

    # For guides with exact matches, calculate distribution metrics
    exact_matches = mapped_df[mapped_df['match_type'] == 'exact']

    # Create full library counts (including zeros for missing guides)
    library_counts = {}
    for seq in library_df['grna_sequence']:
        library_counts[seq] = 0

    for _, row in exact_matches.iterrows():
        library_counts[row['grna_sequence']] = row['count']

    counts_array = np.array(list(library_counts.values()))
    counts_nonzero = counts_array[counts_array > 0]

    # Coverage metrics
    detected_guides = len(counts_nonzero)
    metrics['guides_detected'] = detected_guides
    metrics['coverage_pct'] = 100 * detected_guides / library_size
    metrics['missing_guides'] = library_size - detected_guides
    metrics['zero_count_guides'] = library_size - detected_guides
    metrics['zero_guides_pct'] = 100 * (library_size - detected_guides) / library_size

    if len(counts_nonzero) > 0:
        # Basic stats
        metrics['mean_count'] = np.mean(counts_nonzero)
        metrics['median_count'] = np.median(counts_nonzero)
        metrics['std_count'] = np.std(counts_nonzero)

        # Mean Reads/Guide (including zeros)
        metrics['mean_reads_per_guide'] = total_reads / library_size

        # 90th/10th Ratio (percentile-based, not mean of top/bottom 10%)
        p90 = np.percentile(counts_array, 90)
        p10 = np.percentile(counts_array, 10)
        metrics['p90'] = p90
        metrics['p10'] = p10
        metrics['p90_p10_ratio'] = p90 / p10 if p10 > 0 else float('inf')

        # Log-Count Gini (Gini on log2(counts + 1))
        log_counts = np.log2(counts_array + 1)
        metrics['log_count_gini'] = gini_coefficient(log_counts)

        # Traditional Gini (on raw counts for comparison)
        metrics['gini'] = gini_coefficient(counts_nonzero)

        # Original skew ratio (top 10% mean / bottom 10% mean)
        sorted_counts = np.sort(counts_nonzero)
        n = len(sorted_counts)
        if n >= 10:
            top10 = sorted_counts[-max(1, n//10):].mean()
            bottom10 = sorted_counts[:max(1, n//10)].mean()
            metrics['skew_ratio'] = top10 / bottom10 if bottom10 > 0 else float('inf')
        else:
            metrics['skew_ratio'] = float('nan')

    return metrics


def gini_coefficient(x):
    """Calculate Gini coefficient."""
    x = np.array(x, dtype=float)
    x = np.sort(x)
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    cumsum = np.cumsum(x)
    return (2 * np.sum((np.arange(1, n+1) * x))) / (n * cumsum[-1]) - (n + 1) / n


def determine_verdict(metrics):
    """Determine PASS/FAIL verdict based on thresholds."""
    thresholds = {
        'exact_mapped_pct': ('>=', 65, '%'),
        'mean_reads_per_guide': ('>=', 300, ''),
        'zero_guides_pct': ('<=', 1.0, '%'),
        'p90_p10_ratio': ('<', 10, ''),
        'log_count_gini': ('<=', 0.10, ''),
    }

    verdicts = {}
    overall_pass = True

    for metric, (operator, threshold, unit) in thresholds.items():
        value = metrics.get(metric, 0)

        if operator == '>=':
            passed = value >= threshold
        elif operator == '<=':
            passed = value <= threshold
        elif operator == '<':
            passed = value < threshold
        else:
            passed = False

        verdicts[metric] = {
            'value': value,
            'threshold': threshold,
            'operator': operator,
            'passed': passed,
            'unit': unit
        }

        if not passed:
            overall_pass = False

    return verdicts, overall_pass


def print_mapping_summary(mapping_stats):
    """Print read mapping classification summary."""
    total = mapping_stats['total_reads']

    print("\n" + "="*60)
    print("READ MAPPING CLASSIFICATION")
    print("="*60)
    print(f"{'Category':<35} {'Count':>12} {'Percent':>10}")
    print("-"*60)

    categories = [
        ('Total Reads', 'total_reads'),
        ('Exact Reference Match', 'exact_reference_match'),
        ('1-Mismatch (Unique)', '1_mismatch_unique'),
        ('1-Mismatch (Ambiguous)', '1_mismatch_ambiguous'),
        ('Non-Reference', 'non_reference'),
    ]

    for label, key in categories:
        count = mapping_stats[key]
        pct = 100 * count / total if total > 0 else 0
        print(f"{label:<35} {count:>12,} {pct:>9.2f}%")

    print("="*60)


def print_metrics_report(metrics, verdicts, overall_pass):
    """Print comprehensive metrics report."""
    print("\n" + "="*60)
    print("LIBRARY REPRESENTATION METRICS")
    print("="*60)

    print("\n--- Coverage Metrics ---")
    print(f"  Total reads:                 {metrics['total_reads']:>15,}")
    print(f"  Library size:                {metrics['library_size']:>15,}")
    print(f"  Guides detected:             {metrics['guides_detected']:>15,} ({metrics['coverage_pct']:.2f}%)")
    print(f"  Missing guides:              {metrics['missing_guides']:>15,}")
    print(f"  Zero-count guides (%):       {metrics['zero_guides_pct']:>15.2f}%")

    print("\n--- Mapping Classification ---")
    print(f"  Exact mapped (%):            {metrics['exact_mapped_pct']:>15.2f}%")
    print(f"  1-mismatch unique:           {metrics['1_mismatch_unique']:>15,}")
    print(f"  1-mismatch ambiguous:        {metrics['1_mismatch_ambiguous']:>15,}")
    print(f"  Non-reference:               {metrics['non_reference']:>15,}")

    print("\n--- Distribution Metrics ---")
    print(f"  Mean reads/guide:            {metrics['mean_reads_per_guide']:>15.1f}")
    print(f"  Median count (non-zero):     {metrics['median_count']:>15.1f}")
    print(f"  Std dev:                     {metrics['std_count']:>15.1f}")

    print("\n--- Quality Metrics ---")
    print(f"  90/10 Ratio:                 {metrics['p90_p10_ratio']:>15.2f}")
    print(f"  Log-Count Gini:              {metrics['log_count_gini']:>15.3f}")
    print(f"  Raw Count Gini:              {metrics['gini']:>15.3f}")
    print(f"  Skew ratio (top/bot 10%):    {metrics['skew_ratio']:>15.2f}")

    print("\n" + "="*60)
    print("QUALITY VERDICT")
    print("="*60)

    for metric, result in verdicts.items():
        status = "PASS" if result['passed'] else "FAIL"
        symbol = "✓" if result['passed'] else "✗"
        print(f"  {symbol} {metric:25} {result['value']:.3f} (threshold: {result['operator']} {result['threshold']}{result['unit']}) [{status}]")

    print("-"*60)
    overall_status = "PASS" if overall_pass else "FAIL"
    overall_symbol = "✓✓✓" if overall_pass else "✗✗✗"
    print(f"  {overall_symbol} OVERALL VERDICT: {overall_status}")
    print("="*60)


def main():
    parser = argparse.ArgumentParser(description='Map gRNA counts to library')
    parser.add_argument('counts', help='Input counts file from extract_grna.py')
    parser.add_argument('library', help='Library reference file')
    parser.add_argument('-o', '--output', required=True, help='Output mapped counts file')
    parser.add_argument('-m', '--metrics', required=True, help='Output metrics file')
    args = parser.parse_args()

    print("Loading library reference...")
    library_df = load_library(args.library)
    print(f"  Library size: {len(library_df):,} guides")

    print("\nLoading gRNA counts...")
    counts_df = load_counts(args.counts)
    print(f"  Unique gRNAs in sample: {len(counts_df):,}")
    print(f"  Total reads: {counts_df['count'].sum():,}")

    print("\nMapping to library (with 1-mismatch tolerance)...")
    mapped_df, mapping_stats = map_grnas_with_mismatches(counts_df, library_df)

    print_mapping_summary(mapping_stats)

    # Calculate metrics
    print("\nCalculating metrics...")
    metrics = calculate_metrics(mapped_df, library_df, mapping_stats)

    # Determine verdict
    verdicts, overall_pass = determine_verdict(metrics)

    # Print report
    print_metrics_report(metrics, verdicts, overall_pass)

    # Write outputs
    mapped_df.to_csv(args.output, sep='\t', index=False)
    print(f"\nMapped counts written to {args.output}")

    # Write metrics to TSV
    with open(args.metrics, 'w') as f:
        f.write("metric\tvalue\n")
        for key, value in metrics.items():
            f.write(f"{key}\t{value}\n")
    print(f"Metrics written to {args.metrics}")


if __name__ == '__main__':
    main()
