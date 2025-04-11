#!/usr/bin/env python3

"""Merge reads mapped to terminal repeat regions separately and to the
whole genome construct.
"""

import argparse
from itertools import chain
from collections import defaultdict
import pysam


def score(align):
    """Alignment score.

    It's the length of the aligned query minus the number of mismatches and InDels
    divided by the length of the whole read.

    """
    return (align.query_alignment_length - align.get_tag('NM')) / align.query_length


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-score-diff', required=True, type=float)
    parser.add_argument('--all', required=True)
    parser.add_argument('--left', required=True)
    parser.add_argument('--right', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    aligns_all = pysam.AlignmentFile(args.all, 'r')
    aligns_left = pysam.AlignmentFile(args.left, 'r')
    aligns_right = pysam.AlignmentFile(args.right, 'r')

    out_bam = pysam.AlignmentFile(args.out, "wb", template=aligns_all)

    ref = aligns_all.references[0]
 
    # end pos is the last position + 1 (0 indexed)
    ref_left, range_str_lef = aligns_left.references[0].split(':')
    end_pos_left_terminal_repeat = int(range_str_lef.split('-')[1])

    # start pos is the first position (0 indexed)
    ref_right, range_str_right = aligns_right.references[0].split(':')
    start_pos_right_terminal_repeat = int(range_str_right.split('-')[0]) - 1

    if not ref == ref_left == ref_right:
        raise RuntimeError(f'References ({ref}, {ref_left}, {ref_right} do not match')

    interior_matches = set()
    for align in aligns_all.fetch():
        if align.is_unmapped:
            continue
        matches_between_tr = (
            align.reference_start <  start_pos_right_terminal_repeat
            and align.reference_end > end_pos_left_terminal_repeat
        )
        if matches_between_tr:
            interior_matches.add(align.query_name)
            out_bam.write(align)

    # get the best score for each alignment in aligns_left and aligns_right
    max_scores = defaultdict(float)
    for align in chain(aligns_left.fetch(), aligns_right.fetch()):
        if align.is_unmapped or align.query_name in interior_matches:
            continue
        max_scores[align.query_name] = max(
            max_scores[align.query_name],
            score(align)
        )

    for align in aligns_left.fetch():
        if align.is_unmapped or align.query_name in interior_matches:
            continue
        if score(align) >= args.max_score_diff * max_scores[align.query_name]:
            align.reference_id = out_bam.get_tid(ref)
            align.query_name = align.query_name + '.l'
            out_bam.write(align)

    for align in aligns_right.fetch():
        if align.is_unmapped or align.query_name in interior_matches:
            continue
        if score(align) >= args.max_score_diff * max_scores[align.query_name]:
            align.reference_id = out_bam.get_tid(ref)
            align.query_name = align.query_name + '.r'
            align.reference_start += start_pos_right_terminal_repeat
            out_bam.write(align)

    out_bam.close()
    aligns_all.close()
    aligns_left.close()
    aligns_right.close()


if __name__ == "__main__":
    main()
