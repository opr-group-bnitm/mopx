#!/usr/bin/env python3

"""Turn secondary into primary alignments.
"""

import argparse
from collections import Counter
from Bio.Seq import Seq
import pysam


def orient(seq: str, qual, forward: bool) -> str:
    if forward:
        return seq, qual
    return str(Seq(seq).reverse_complement()), qual[::-1]


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', help='Input bam file.', required=True)
    parser.add_argument('--output', '-o', help='Output bam file.', required=True)
    args = parser.parse_args()

    with pysam.AlignmentFile(args.input, 'rb') as infile:
        secondary_ids = set(
            align.query_name
            for align in infile.fetch()
            if align.is_secondary
        )
        sequences = {
            align.query_name: orient(align.query_sequence, align.query_qualities, align.is_forward)
            for align in infile.fetch()
            if align.query_name in secondary_ids and not (align.is_secondary or align.is_supplementary)
        }

        align_id_tracker = Counter()

        with pysam.AlignmentFile(args.output, "wb", header=infile.header) as outfile:
            for alignment in infile.fetch():
                qname = alignment.query_name
                if alignment.is_unmapped:
                    continue
                if alignment.is_secondary:
                    align_id_tracker.update([qname])
                    seq, qual = orient(*sequences[alignment.query_name], alignment.is_forward)
                    
                    alignment.query_sequence = seq
                    alignment.query_qualities = qual
                    alignment.is_secondary = False
                    alignment.query_name = f'{qname}.{align_id_tracker[qname]}'

                    if alignment.has_tag("tp") and alignment.get_tag("tp") == "S":
                        # Change the 'tp' tag value to 'A:P'
                        alignment.set_tag("tp", "P")

                    outfile.write(alignment)
                else:
                    outfile.write(alignment)

    print(f'wrote alignments to {args.output}')


if __name__ == '__main__':
    main()
