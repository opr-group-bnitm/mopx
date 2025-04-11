#!/usr/bin/env python3

import argparse
import pysam
from collections import Counter, defaultdict
import pandas as pd


class OperationCollector:

    operations = {
        0: 'match_mismatch',
        1: 'insertion',
        2: 'deletion',
        3: 'Skip',
        4: 'softclip',
        5: 'hardclip',
        6: 'padding',
        7: 'match',
        8: 'mismatch'
    }

    def __init__(self, align):
        self.readid = align.query_name
        self.max_counts = defaultdict(int)
        self.total_counts = Counter()
        self.occurences = Counter()
        if align.is_unmapped:
            self.aligntype = 'UNMAPPED'
        else:
            self.aligntype = 'PRIMARY'
            if align.is_secondary:
                self.aligntype = 'SECONDARY'
            if align.is_supplementary:
                self.aligntype = 'SUPPLEMENTARY'
            for operation, length in align.cigar:    
                self.total_counts[operation] += length
                self.occurences[operation] += 1
                self.max_counts[operation] = max(self.max_counts[operation], length)

    def table_column_dict(self):
        row = {
            'Read': self.readid,
            'Type': self.aligntype,
        }
        for operation, operation_name in self.operations.items():
            row[f'{operation_name}_total_nucleotides'] = self.total_counts[operation]
            row[f'{operation_name}_occurences'] = self.occurences[operation]
            row[f'{operation_name}_max_length'] = self.max_counts[operation]
        return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True, help='Mapped contigs in .bam format.')
    parser.add_argument('--out', required=True, help='Output file in .tsv format.')
    args = parser.parse_args()

    with pysam.AlignmentFile(args.bam, 'rb') as bam_file:
        operation_stats = [
            OperationCollector(alignment).table_column_dict()
            for alignment in bam_file.fetch()
        ]
    pd.DataFrame(operation_stats).to_csv(args.out, sep='\t')


if __name__ == "__main__":
    main()
