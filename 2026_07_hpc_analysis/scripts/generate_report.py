#!/usr/bin/env python3
"""
Generate consolidated QC report for CRISPRi library analysis.
Creates markdown report with executive summary, classification breakdown, and quality verdict.
"""

import argparse
import json
import pandas as pd
from pathlib import Path
from datetime import datetime


def load_metrics(metrics_path):
    """Load metrics TSV file into dict."""
    df = pd.read_csv(metrics_path, sep='\t', header=None, names=['metric', 'value'])
    return dict(zip(df['metric'], df['value']))


def load_classification(classification_path):
    """Load classification JSON file."""
    with open(classification_path, 'r') as f:
        return json.load(f)


def determine_verdict(metrics):
    """Determine PASS/FAIL verdict based on thresholds."""
    thresholds = {
        'exact_mapped_pct': ('>=', 65),
        'mean_reads_per_guide': ('>=', 300),
        'zero_guides_pct': ('<=', 1.0),
        'p90_p10_ratio': ('<', 10),
        'log_count_gini': ('<=', 0.10),
    }

    verdicts = {}
    overall_pass = True

    for metric, (operator, threshold) in thresholds.items():
        value = float(metrics.get(metric, 0))

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
            'passed': passed
        }

        if not passed:
            overall_pass = False

    return verdicts, overall_pass


def format_number(value, decimals=1):
    """Format number with commas and specified decimals."""
    if isinstance(value, (int, float)):
        if decimals == 0:
            return f"{int(value):,}"
        return f"{value:,.{decimals}f}"
    return str(value)


def generate_executive_summary(set_a_metrics, set_b_metrics):
    """Generate executive summary table."""
    rows = []

    for set_name, metrics in [('A', set_a_metrics), ('B', set_b_metrics)]:
        if not metrics:
            rows.append({
                'Set': set_name,
                'Status': 'N/A',
                'Total Reads': 'N/A',
                'Mapped (%)': 'N/A',
                'Mean Reads': 'N/A',
                'Median Reads': 'N/A',
                'Zero Guides (%)': 'N/A',
                '90/10 Ratio': 'N/A',
                'Log Gini': 'N/A'
            })
            continue

        verdicts, overall_pass = determine_verdict(metrics)
        status = 'PASS' if overall_pass else 'FAIL'

        rows.append({
            'Set': set_name,
            'Status': status,
            'Total Reads': format_number(metrics.get('total_reads', 0), 0),
            'Mapped (%)': format_number(metrics.get('exact_mapped_pct', 0), 2),
            'Mean Reads': format_number(metrics.get('mean_reads_per_guide', 0), 1),
            'Median Reads': format_number(metrics.get('median_count', 0), 1),
            'Zero Guides (%)': format_number(metrics.get('zero_guides_pct', 0), 2),
            '90/10 Ratio': format_number(metrics.get('p90_p10_ratio', 0), 2),
            'Log Gini': format_number(metrics.get('log_count_gini', 0), 3)
        })

    return pd.DataFrame(rows)


def generate_read_breakdown(set_a_class, set_b_class):
    """Generate detailed read classification breakdown table."""
    categories = [
        ('Total Reads', 'total_reads'),
        ('Valid 20-nt Insert', 'valid_20nt_insert'),
        ('Missing Prefix (no CACCG)', 'missing_prefix'),
        ('Missing Downstream Suffix (no GTTT)', 'missing_suffix'),
        ('Too Short After Prefix', 'too_short'),
        ('Low Quality', 'low_quality'),
    ]

    rows = []
    for label, key in categories:
        row = {'Category': label}

        for set_name, class_data in [('Set A', set_a_class), ('Set B', set_b_class)]:
            if class_data and key in class_data:
                count = class_data[key]
                total = class_data.get('total_reads', 1)
                pct = 100 * count / total if total > 0 else 0
                row[f'{set_name} Count'] = format_number(count, 0)
                row[f'{set_name} %'] = f"{pct:.2f}%"
            else:
                row[f'{set_name} Count'] = 'N/A'
                row[f'{set_name} %'] = 'N/A'

        rows.append(row)

    return pd.DataFrame(rows)


def generate_mapping_breakdown(set_a_metrics, set_b_metrics):
    """Generate mapping classification breakdown table."""
    categories = [
        ('Total Classified Reads', 'total_reads'),
        ('Exact Reference Match', 'exact_reference_match'),
        ('1-Mismatch (Unique)', '1_mismatch_unique'),
        ('1-Mismatch (Ambiguous)', '1_mismatch_ambiguous'),
        ('Non-Reference', 'non_reference'),
    ]

    rows = []
    for label, key in categories:
        row = {'Category': label}

        for set_name, metrics in [('Set A', set_a_metrics), ('Set B', set_b_metrics)]:
            if metrics and key in metrics:
                count = metrics[key]
                total = metrics.get('total_reads', 1)
                pct = 100 * count / total if total > 0 else 0
                row[f'{set_name} Count'] = format_number(count, 0)
                row[f'{set_name} %'] = f"{pct:.2f}%"
            else:
                row[f'{set_name} Count'] = 'N/A'
                row[f'{set_name} %'] = 'N/A'

        rows.append(row)

    return pd.DataFrame(rows)


def generate_verdict_section(set_a_metrics, set_b_metrics):
    """Generate quality verdict and next steps section."""
    sections = []

    for set_name, metrics in [('Set A', set_a_metrics), ('Set B', set_b_metrics)]:
        if not metrics:
            sections.append(f"## Set {set_name[-1]}: No Data Available\n")
            continue

        verdicts, overall_pass = determine_verdict(metrics)
        status = "PASS" if overall_pass else "FAIL"
        status_icon = "✅" if overall_pass else "❌"

        sections.append(f"## Set {set_name[-1]}: {status_icon} {status}\n")

        sections.append("### Metric Evaluation\n")
        sections.append("| Metric | Value | Threshold | Status |")
        sections.append("|--------|-------|-----------|--------|")

        metric_labels = {
            'exact_mapped_pct': 'Exact Mapped (%)',
            'mean_reads_per_guide': 'Mean Reads/Guide',
            'zero_guides_pct': 'Zero Guides (%)',
            'p90_p10_ratio': '90/10 Ratio',
            'log_count_gini': 'Log-Count Gini'
        }

        for metric, result in verdicts.items():
            label = metric_labels.get(metric, metric)
            value = format_number(result['value'], 3)
            threshold = f"{result['operator']} {result['threshold']}"
            m_status = "✅ PASS" if result['passed'] else "❌ FAIL"
            sections.append(f"| {label} | {value} | {threshold} | {m_status} |")

        sections.append("")

        # Warnings and recommendations
        warnings = []

        em_pct = float(metrics.get('exact_mapped_pct', 0))
        if em_pct < 65:
            warnings.append(f"- **Low mapping rate** ({em_pct:.1f}% < 65%): Consider sequencing depth or library quality")

        zero_pct = float(metrics.get('zero_guides_pct', 0))
        if zero_pct > 1.0:
            warnings.append(f"- **High dropout rate** ({zero_pct:.2f}% > 1%): Significant guide loss detected")

        ratio = float(metrics.get('p90_p10_ratio', 0))
        if ratio >= 10:
            warnings.append(f"- **High 90/10 ratio** ({ratio:.2f} >= 10): Uneven guide representation")

        log_gini = float(metrics.get('log_count_gini', 0))
        if log_gini > 0.10:
            warnings.append(f"- **High inequality** (Log Gini {log_gini:.3f} > 0.10): Library is skewed")

        if warnings:
            sections.append("### ⚠️ Warnings\n")
            sections.extend(warnings)
        else:
            sections.append("### ✅ No Critical Issues\n")
            sections.append("All quality metrics are within acceptable thresholds.")

        sections.append("")

    return "\n".join(sections)


def generate_full_report(set_a_metrics, set_b_metrics, set_a_class, set_b_class):
    """Generate full markdown report."""
    report = []

    # Header
    report.append("# CRISPRi Library QC Analysis Report")
    report.append("")
    report.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("")

    # Executive Summary
    report.append("---")
    report.append("")
    report.append("## Executive Summary\n")

    summary_df = generate_executive_summary(set_a_metrics, set_b_metrics)
    report.append(summary_df.to_markdown(index=False))
    report.append("")

    # Mapping Classification Breakdown
    report.append("---")
    report.append("")
    report.append("## Detailed Mapping Classification Breakdown\n")

    mapping_df = generate_mapping_breakdown(set_a_metrics, set_b_metrics)
    report.append(mapping_df.to_markdown(index=False))
    report.append("")

    # Read Classification Breakdown
    report.append("---")
    report.append("")
    report.append("## Raw Read Classification Breakdown\n")
    report.append("(From FASTQ processing)\n")

    read_df = generate_read_breakdown(set_a_class, set_b_class)
    report.append(read_df.to_markdown(index=False))
    report.append("")

    # Quality Verdict
    report.append("---")
    report.append("")
    report.append("## Quality Verdict & Next Steps\n")

    verdict_section = generate_verdict_section(set_a_metrics, set_b_metrics)
    report.append(verdict_section)

    # Footer
    report.append("---")
    report.append("")
    report.append("## Notes\n")
    report.append("- **Exact Mapped (%)**: Percentage of total reads that map exactly to reference library")
    report.append("- **90/10 Ratio**: Ratio of 90th percentile to 10th percentile guide counts (lower is more even)")
    report.append("- **Log-Count Gini**: Gini coefficient computed on log2(counts + 1) (lower is more even)")
    report.append("- **Flanking sequences**: CACCG (5') and GTTT (3') per Dolcetto requirements")
    report.append("")

    return "\n".join(report)


def main():
    parser = argparse.ArgumentParser(description='Generate consolidated QC report')
    parser.add_argument('--set-a-metrics', help='Set A metrics TSV file')
    parser.add_argument('--set-b-metrics', help='Set B metrics TSV file')
    parser.add_argument('--set-a-class', help='Set A classification JSON file')
    parser.add_argument('--set-b-class', help='Set B classification JSON file')
    parser.add_argument('-o', '--output', required=True, help='Output markdown report file')
    args = parser.parse_args()

    # Load metrics
    set_a_metrics = load_metrics(args.set_a_metrics) if args.set_a_metrics else None
    set_b_metrics = load_metrics(args.set_b_metrics) if args.set_b_metrics else None

    # Load classifications
    set_a_class = load_classification(args.set_a_class) if args.set_a_class else None
    set_b_class = load_classification(args.set_b_class) if args.set_b_class else None

    # Generate report
    report = generate_full_report(set_a_metrics, set_b_metrics, set_a_class, set_b_class)

    # Write report
    with open(args.output, 'w') as f:
        f.write(report)

    print(f"Report written to {args.output}")


if __name__ == '__main__':
    main()
