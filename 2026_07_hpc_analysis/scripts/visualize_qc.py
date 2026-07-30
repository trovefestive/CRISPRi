#!/usr/bin/env python3
"""
Generate QC visualizations for CRISPRi library screening.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path


def load_mapped_counts(path):
    """Load mapped counts file."""
    return pd.read_csv(path, sep='\t')


def plot_count_distribution(df, output_path, set_name):
    """Plot log2 count distribution histogram."""
    in_lib = df[(df['in_library'] == True) & (df['count'] > 0)]
    counts = in_lib['count'].values

    # Calculate log2 fold change from expected (assuming even distribution)
    total_reads = df['count'].sum()
    library_size = len(df[df['in_library'] == True])
    expected = total_reads / library_size
    log2fc = np.log2(counts / expected)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Raw counts histogram
    ax1 = axes[0]
    ax1.hist(np.log10(counts), bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.set_xlabel('Log₁₀ Count', fontsize=12)
    ax1.set_ylabel('Number of Guides', fontsize=12)
    ax1.set_title(f'{set_name}: Raw Count Distribution', fontsize=14)
    ax1.axvline(np.log10(expected), color='red', linestyle='--', linewidth=2,
                label=f'Expected (even dist): {expected:.0f}')
    ax1.legend()
    ax1.grid(alpha=0.3)

    # Log2 fold change histogram
    ax2 = axes[1]
    ax2.hist(log2fc, bins=50, edgecolor='black', alpha=0.7, color='forestgreen')
    ax2.set_xlabel('Log₂ Fold-Change (Observed/Expected)', fontsize=12)
    ax2.set_ylabel('Number of Guides', fontsize=12)
    ax2.set_title(f'{set_name}: Representation vs Expected', fontsize=14)
    ax2.axvline(0, color='red', linestyle='--', linewidth=2, label='Expected')
    ax2.axvline(np.median(log2fc), color='orange', linestyle='--', linewidth=2,
                label=f'Median: {np.median(log2fc):.2f}')
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Distribution plot saved to {output_path}")


def plot_lorenz_curve(df, output_path, set_name):
    """Plot Lorenz curve for inequality visualization."""
    in_lib = df[(df['in_library'] == True) & (df['count'] > 0)]
    counts = np.sort(in_lib['count'].values)

    # Calculate cumulative distribution
    n = len(counts)
    cumsum = np.cumsum(counts)
    cumsum_norm = cumsum / cumsum[-1]
    x = np.arange(1, n + 1) / n

    # Gini coefficient
    gini = 2 * np.sum((np.arange(1, n + 1) * counts)) / (n * cumsum[-1]) - (n + 1) / n

    fig, ax = plt.subplots(figsize=(8, 8))

    # Perfect equality line
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Perfect Equality')

    # Lorenz curve
    ax.plot(x, cumsum_norm, linewidth=2, color='steelblue', label=f'Observed (Gini = {gini:.3f})')

    # Fill area between curves
    ax.fill_between(x, cumsum_norm, x, alpha=0.3, color='steelblue')

    ax.set_xlabel('Cumulative Fraction of Guides', fontsize=12)
    ax.set_ylabel('Cumulative Fraction of Reads', fontsize=12)
    ax.set_title(f'{set_name}: Lorenz Curve (Library Inequality)', fontsize=14)
    ax.legend(loc='upper left')
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Lorenz curve saved to {output_path}")


def plot_rank_abundance(df, output_path, set_name):
    """Plot rank-abundance curve."""
    in_lib = df[(df['in_library'] == True) & (df['count'] > 0)]
    counts = np.sort(in_lib['count'].values)[::-1]  # Sort descending
    ranks = np.arange(1, len(counts) + 1)

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(ranks, counts, linewidth=1.5, color='steelblue')
    ax.set_xlabel('Guide Rank', fontsize=12)
    ax.set_ylabel('Read Count', fontsize=12)
    ax.set_title(f'{set_name}: Rank-Abundance Plot', fontsize=14)
    ax.set_yscale('log')
    ax.grid(alpha=0.3)

    # Mark top guides
    top_10_idx = len(counts) // 10
    ax.axvline(top_10_idx, color='red', linestyle='--', alpha=0.5, label=f'Top 10% boundary (rank {top_10_idx:,})')
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Rank-abundance plot saved to {output_path}")


def plot_gene_coverage(df, output_path, set_name):
    """Plot coverage statistics by gene."""
    in_lib = df[df['in_library'] == True].copy()

    # Group by gene
    gene_stats = in_lib.groupby('gene_symbol').agg({
        'count': ['sum', 'mean', 'std', 'count']
    }).reset_index()
    gene_stats.columns = ['gene_symbol', 'total_counts', 'mean_count', 'std_count', 'n_guides']

    # Calculate guides detected per gene
    gene_detected = in_lib[in_lib['count'] > 0].groupby('gene_symbol').size().reset_index(name='detected_guides')
    gene_stats = gene_stats.merge(gene_detected, on='gene_symbol', how='left')
    gene_stats['detected_guides'] = gene_stats['detected_guides'].fillna(0).astype(int)
    gene_stats['coverage_pct'] = 100 * gene_stats['detected_guides'] / gene_stats['n_guides']

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Guides per gene histogram
    ax1 = axes[0, 0]
    ax1.hist(gene_stats['n_guides'], bins=30, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.set_xlabel('Number of Guides per Gene')
    ax1.set_ylabel('Number of Genes')
    ax1.set_title('Guide Count Distribution per Gene')
    ax1.grid(alpha=0.3)

    # Coverage % histogram
    ax2 = axes[0, 1]
    ax2.hist(gene_stats['coverage_pct'], bins=30, edgecolor='black', alpha=0.7, color='forestgreen')
    ax2.set_xlabel('% Guides Detected per Gene')
    ax2.set_ylabel('Number of Genes')
    ax2.set_title('Gene Coverage Distribution')
    ax2.axvline(gene_stats['coverage_pct'].mean(), color='red', linestyle='--',
                label=f'Mean: {gene_stats["coverage_pct"].mean():.1f}%')
    ax2.legend()
    ax2.grid(alpha=0.3)

    # Total counts per gene histogram
    ax3 = axes[1, 0]
    ax3.hist(np.log10(gene_stats['total_counts'] + 1), bins=50, edgecolor='black', alpha=0.7, color='coral')
    ax3.set_xlabel('Log₁₀ Total Counts per Gene')
    ax3.set_ylabel('Number of Genes')
    ax3.set_title('Read Count Distribution per Gene')
    ax3.grid(alpha=0.3)

    # Scatter: coverage % vs mean count
    ax4 = axes[1, 1]
    scatter = ax4.scatter(gene_stats['coverage_pct'], gene_stats['mean_count'],
                          alpha=0.5, c=gene_stats['n_guides'], cmap='viridis', s=20)
    ax4.set_xlabel('% Guides Detected per Gene')
    ax4.set_ylabel('Mean Count per Guide')
    ax4.set_title('Coverage vs Abundance')
    ax4.set_yscale('log')
    plt.colorbar(scatter, ax=ax4, label='Guides per Gene')
    ax4.grid(alpha=0.3)

    plt.suptitle(f'{set_name}: Gene-Level QC Statistics', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Gene coverage plot saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Generate QC visualizations')
    parser.add_argument('mapped_counts', help='Mapped counts file from map_grna.py')
    parser.add_argument('-o', '--output-dir', default='figures', help='Output directory for figures')
    parser.add_argument('-n', '--name', default='Sample', help='Sample/set name for titles')
    args = parser.parse_args()

    Path(args.output_dir).mkdir(exist_ok=True)

    print(f"Loading mapped counts from {args.mapped_counts}...")
    df = load_mapped_counts(args.mapped_counts)

    print(f"\nGenerating plots for {args.name}...")

    # Generate all plots
    plot_count_distribution(df, f"{args.output_dir}/{args.name}_distribution.png", args.name)
    plot_lorenz_curve(df, f"{args.output_dir}/{args.name}_lorenz.png", args.name)
    plot_rank_abundance(df, f"{args.output_dir}/{args.name}_rank_abundance.png", args.name)
    plot_gene_coverage(df, f"{args.output_dir}/{args.name}_gene_coverage.png", args.name)

    print("\nAll plots generated successfully!")


if __name__ == '__main__':
    main()
