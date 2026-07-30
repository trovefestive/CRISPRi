#!/usr/bin/env python3
"""
Extract gRNA spacer sequences from FASTQ files.
Looks for gRNA sequence between 5' constant (CACCG) and 3' constant (GTTT).

Per requirements: Flanking sequences are CACCG (5') and GTTT (3'),
with exactly 20 nucleotides between them.
"""

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path


def open_fastq(path):
    """Open FASTQ file (gzipped or not)."""
    path = Path(path)
    if path.suffix == '.gz':
        return gzip.open(path, 'rt')
    return open(path, 'r')


def classify_read(sequence, upstream="CACCG", downstream="GTTT", length=20):
    """
    Classify a read based on presence of flanking sequences.

    Returns a tuple: (classification, extracted_sequence)
    - classification: one of 'valid', 'missing_prefix', 'missing_suffix', 'too_short'
    - extracted_sequence: the 20-nt spacer if valid, else None
    """
    # Look for upstream constant region
    up_pos = sequence.find(upstream)
    if up_pos == -1:
        return 'missing_prefix', None

    # gRNA starts after upstream constant
    grna_start = up_pos + len(upstream)
    grna_end = grna_start + length

    # Check if we have enough sequence
    if grna_end > len(sequence):
        return 'too_short', None

    spacer = sequence[grna_start:grna_end]

    # Verify downstream constant region is present
    downstream_seq = sequence[grna_end:grna_end + len(downstream)]
    if downstream_seq != downstream:
        return 'missing_suffix', None

    return 'valid', spacer


def count_grnas(fastq_path, min_quality=20):
    """
    Count gRNA spacers from FASTQ file with detailed classification.

    Returns:
        - counter: Counter of gRNA sequences (only valid, quality-passing reads)
        - classification: dict with read classification counts
    """
    counter = Counter()

    # Read classification tracking
    classification = {
        'total_reads': 0,
        'missing_prefix': 0,
        'too_short': 0,
        'missing_suffix': 0,
        'low_quality': 0,
        'valid_20nt_insert': 0,
    }

    with open_fastq(fastq_path) as f:
        while True:
            # Read 4 lines (one FASTQ record)
            header = f.readline()
            if not header:
                break
            sequence = f.readline().strip()
            plus = f.readline()
            quality = f.readline().strip()

            classification['total_reads'] += 1

            # Classify the read based on flanking sequences
            read_class, spacer = classify_read(sequence)

            if read_class == 'missing_prefix':
                classification['missing_prefix'] += 1
                continue
            elif read_class == 'too_short':
                classification['too_short'] += 1
                continue
            elif read_class == 'missing_suffix':
                classification['missing_suffix'] += 1
                continue

            # At this point we have a valid 20-nt insert
            # Apply quality filter (Phred+33 encoding)
            grna_start = sequence.find("CACCG") + len("CACCG")
            grna_end = grna_start + 20
            grna_quality = quality[grna_start:grna_end]

            if min_quality > 0:
                q_scores = [ord(q) - 33 for q in grna_quality]
                if min(q_scores) < min_quality:
                    classification['low_quality'] += 1
                    continue

            # Valid, quality-passing read
            classification['valid_20nt_insert'] += 1
            counter[spacer] += 1

            if classification['total_reads'] % 1000000 == 0:
                valid_pct = 100 * classification['valid_20nt_insert'] / classification['total_reads']
                print(f"Processed {classification['total_reads']:,} reads, "
                      f"{classification['valid_20nt_insert']:,} valid ({valid_pct:.1f}%)")

    return counter, classification


def print_classification_summary(classification):
    """Print read classification summary."""
    total = classification['total_reads']

    print("\n" + "="*60)
    print("READ CLASSIFICATION SUMMARY")
    print("="*60)
    print(f"{'Category':<35} {'Count':>12} {'Percent':>10}")
    print("-"*60)

    categories = [
        ('Total Reads', 'total_reads'),
        ('Missing Prefix (no CACCG)', 'missing_prefix'),
        ('Too Short After Prefix', 'too_short'),
        ('Missing Downstream Suffix (no GTTT)', 'missing_suffix'),
        ('Low Quality (< min Phred)', 'low_quality'),
        ('Valid 20-nt Insert', 'valid_20nt_insert'),
    ]

    for label, key in categories:
        count = classification[key]
        pct = 100 * count / total if total > 0 else 0
        print(f"{label:<35} {count:>12,} {pct:>9.2f}%")

    # Verification
    calculated_total = (classification['missing_prefix'] +
                       classification['too_short'] +
                       classification['missing_suffix'] +
                       classification['low_quality'] +
                       classification['valid_20nt_insert'])

    print("-"*60)
    print(f"{'Verification (should equal Total)':<35} {calculated_total:>12,}")
    print("="*60)


def main():
    parser = argparse.ArgumentParser(description='Extract gRNA counts from FASTQ')
    parser.add_argument('fastq', help='Input FASTQ file (.fq or .fq.gz)')
    parser.add_argument('-o', '--output', required=True, help='Output counts file')
    parser.add_argument('-c', '--classification', help='Output classification summary file (JSON)')
    parser.add_argument('--min-quality', type=int, default=20,
                        help='Minimum base quality (Phred, default: 20)')
    args = parser.parse_args()

    print(f"Processing {args.fastq}...")
    print(f"Using flanking sequences: 5'-CACCG-[20nt]-GTTT-3'")
    print(f"Minimum quality threshold: {args.min_quality}")
    print()

    counts, classification = count_grnas(args.fastq, min_quality=args.min_quality)

    # Print classification summary
    print_classification_summary(classification)

    # Write counts to file
    with open(args.output, 'w') as f:
        f.write("grna_sequence\tcount\n")
        for grna, count in counts.most_common():
            f.write(f"{grna}\t{count}\n")

    print(f"\nResults written to {args.output}")
    print(f"  Unique gRNAs: {len(counts):,}")
    print(f"  Total valid reads: {sum(counts.values()):,}")

    # Write classification summary to JSON if requested
    if args.classification:
        with open(args.classification, 'w') as f:
            json.dump(classification, f, indent=2)
        print(f"  Classification summary: {args.classification}")


if __name__ == '__main__':
    main()
