#!/usr/bin/env python3

import argparse
import json
import os
import pandas as pd

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
        '--basecounts',
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
        help='Output directory'
    )
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    fname_seq_out = os.path.join(args.out, 'fixed_hairpin_draft.fasta')
    fname_basecounts_out = os.path.join(args.out, 'fixed_hairpin_basecounts.tsv')

    with open(args.hairpin_matches) as f:
        match_details = json.load(f)

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
        fname_seq_out,
        'fasta'
    )

    basecounts = pd.read_csv(args.basecounts, sep='\t')

    shift_pos = len(arm_1.seq) - len_left
    basecounts['position'] += shift_pos
    basecounts = basecounts[len_left:-len_right]

    basecounts.to_csv(fname_basecounts_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
