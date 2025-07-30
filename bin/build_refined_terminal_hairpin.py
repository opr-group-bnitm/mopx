#!/usr/bin/env python3

import argparse
import pysam
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description='Build a consensus FASTA from a Medaka gVCF with specific coordinates for the mpox non-complementary hairpin arms.'
    )
    parser.add_argument(
        '--gvcf',
        required=True,
        help='Input gVCF file (Medaka output, annotated with DP).',
    )
    parser.add_argument(
        '--hairpin-matches',
        required=True,
        help='Json file with coordinates of hairpin.',
    )
    parser.add_argument(
        '--min-depth',
        type=int,
        default=20,
        help='Minimum read depth (DP) to keep a site unmasked. Default: 20.',
    )
    parser.add_argument(
        '--min-qual',
        type=int,
        default=20,
        help='Minimum variant call quality (QUAL) to keep a site unmasked. Default: 20.',
    )
    parser.add_argument(
        '--out',
        default='consensus_hairpin.fasta',
        help='Output consensus FASTA path. Default: consensus.fasta',
    )
    return parser.parse_args()


def read_hairpin_match_details(fname):
    with open(fname) as f:
        return json.load(f)


def read_ref(fname):
    return SeqIO.read(fname, 'fasta')


def read_vcf(fname):
    return list(pysam.VariantFile(fname))


def filter_variants(variants, min_qual, min_depth):
    return [
        v
        for v in variants
        if v.info.get("DP", 0) >= min_depth and v.qual is not None and v.qual >= min_qual
    ]


def check_variants(variants, pos1, pos2, pos3):
    for v in variants:
        start = v.pos
        end = v.pos + len(v.ref) - 1
        if start != end:
            for pos in [pos1, pos2, pos3]:
                if start < pos and end >= pos:
                    return False
    return True


def build_consensus(variants, start, end):
    ref_called = {
        v.pos - 1: v.ref
        for v in variants
        if v.alts is None
    }
    var = {
        v.pos - 1: (len(v.ref), v.alts[0])
        for v in variants
        if v.alts is not None
    }
    cons = []
    i = start
    while i < end:
        if i in var:
            varlen, alt = var[i]
            cons.append(alt)
            i += varlen
        elif i in ref_called:
            cons.append(ref_called[i])
            i += 1
        else:
            cons.append('N')
            i += 1
    return ''.join(cons)


def write_consensus(fname, arm_1, arm_2):
    records = [
        SeqRecord(
            id='hairpin_arm_1',
            name='hairpin_arm_1',
            seq=Seq(arm_1).reverse_complement(),
            description=''
        ),
        SeqRecord(
            id='hairpin_arm_2',
            name='hairpin_arm_2',
            seq=Seq(arm_2),
            description=''
        )
    ]
    SeqIO.write(records, fname, 'fasta')


def main():

    args = parse_args()

    matches = read_hairpin_match_details(args.hairpin_matches)

    margin = matches['margin_len']
    hp_len_1 = len(matches['arm_1'])
    hp_len_2 = len(matches['arm_2'])

    pos1 = margin
    pos2 = margin + hp_len_1
    pos3 = pos2 + hp_len_2

    variants = read_vcf(args.gvcf)
    filtered = filter_variants(variants, args.min_qual, args.min_depth)

    variants_are_ok = check_variants(variants, pos1, pos2, pos3)

    if not variants_are_ok:
        print('Failed due to complex variants')
        exit(0)

    consensus_arm_1 = build_consensus(filtered, pos1, pos2)
    consensus_arm_2 = build_consensus(filtered, pos2, pos3)

    write_consensus(args.out, consensus_arm_1, consensus_arm_2)


if __name__ == '__main__':
    main()
