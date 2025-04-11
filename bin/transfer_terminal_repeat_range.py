#!/usr/bin/env python3

"""Transfer the end position of the left terminal repeat and the start position
of the right terminal repeat from a reference sequence onto the query sequence.
"""

from Bio import SeqIO
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True)
    parser.add_argument('--alignment', required=True)
    parser.add_argument(
        '--trl-end', required=True, type=int,
        help='Last position of the left terminal repeat, (1-based index)'
    )
    parser.add_argument(
        '--trr-start', required=True, type=int,
        help='First position of the right terminal repeat (1-based index)'
    )
    args = parser.parse_args()

    sequences = [
        str(record.seq)
        for record in SeqIO.parse(args.alignment, "fasta")
    ]
    ref = sequences[0]
    query = sequences[1]

    ref_pos = [
        i
        for i, s in enumerate(ref)
        if s != '-'
    ]

    last_query_pos = -1
    query_pos = []
    for q in query:
        if q != '-':
            last_query_pos += 1
        query_pos.append(last_query_pos)

    left_column = ref_pos[args.trl_end - 1]
    left_pos_query = query_pos[left_column]

    right_column = ref_pos[args.trr_start]
    right_pos_query = query_pos[right_column]
    # in case of a gap, the next matching base is chosen as
    # the start of the terminal repeat region
    if query[right_column] == '-':
        right_pos_query += 1

    with open(args.out, 'w') as f_out:
        # + 1 because the output is 1 based index
        f_out.write(f'{left_pos_query + 1}\n{right_pos_query + 1}\n')


if __name__ == "__main__":
    main()
