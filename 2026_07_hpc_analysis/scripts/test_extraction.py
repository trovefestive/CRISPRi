#!/usr/bin/env python3
"""
Quick test script to verify gRNA extraction on a small FASTQ sample.
"""

import sys
sys.path.insert(0, '.')
from extract_grna import extract_grna, open_fastq

def test_extraction(fastq_path, n_reads=10000):
    """Test gRNA extraction on first N reads."""
    print(f"Testing extraction on first {n_reads} reads from {fastq_path}...")

    extracted = []
    total = 0

    with open_fastq(fastq_path) as f:
        while total < n_reads:
            header = f.readline()
            if not header:
                break
            sequence = f.readline().strip()
            plus = f.readline()
            quality = f.readline()

            total += 1
            grna = extract_grna(sequence)
            if grna:
                extracted.append(grna)

            if total <= 10:
                print(f"\nRead {total}:")
                print(f"  Sequence: {sequence[:80]}...")
                print(f"  gRNA: {grna}")

    print(f"\n{'='*50}")
    print(f"Test Results:")
    print(f"  Total reads: {total}")
    print(f"  Valid gRNAs extracted: {len(extracted)} ({100*len(extracted)/total:.1f}%)")
    print(f"  Unique gRNAs: {len(set(extracted))}")

    if extracted:
        print(f"\n  First 10 extracted sequences:")
        for seq in extracted[:10]:
            print(f"    {seq}")

    return len(extracted) / total if total > 0 else 0


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Test gRNA extraction')
    parser.add_argument('fastq', help='FASTQ file to test')
    parser.add_argument('-n', '--num-reads', type=int, default=10000,
                        help='Number of reads to test (default: 10000)')
    args = parser.parse_args()

    success_rate = test_extraction(args.fastq, args.num_reads)
    print(f"\nSuccess rate: {success_rate*100:.1f}%")
