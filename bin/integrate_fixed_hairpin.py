#!/usr/bin/env python3

import argparse
import json

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--draft',
        required=True,
        help=''
    )
    parser.add_argument(
        '--hairpin-consensus',
        required=True,
        help=''
    )
    parser.add_argument(
        '--hairpin-matches',
        required=True,
        help='JSON file with information about initial hairpin to draft genome alignment.'
    )
    parser.add_argument(
        '--out',
        required=True,
        help='Fasta filename for corrected draft.'
    )
    args = parser.parse_args()

    with open(args.hairpin_matches) as f:
        match_details =  json.load(f)

    arm_1 = next(SeqIO.parse(args.hairpin_consensus, 'fasta'))
    draft = SeqIO.read(args.draft, 'fasta')

    len_left = match_details['genome_hairpin_len_left']
    len_right = match_details['genome_hairpin_len_right']

    outseq = arm_1.seq + draft.seq[len_left:-len_right] + arm_1.seq.reverse_complement()

    SeqIO.write(
        SeqRecord(
            id=draft.id,
            description=draft.description,
            name=draft.name,
            seq=outseq
        ),
        args.out,
        'fasta'
    )


if __name__ == '__main__':
    main()
