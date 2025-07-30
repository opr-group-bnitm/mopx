#!/usr/bin/env python3

import argparse
import os
import json

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner


def get_dna_aligner(right_gap_score=-1):
    aligner = PairwiseAligner()
    # Load a DNA-specific substitution matrix, e.g., IUB for nucleotides
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    aligner.mode = 'global'
    aligner.right_open_gap_score = right_gap_score
    aligner.right_extend_gap_score = right_gap_score
    aligner.match_score = 1
    aligner.mismatch_score = -1
    return aligner


def align_end(
        arm_1,
        arm_2,
        genome,
        arms_selfscore,
        margin
):
    aligner = get_dna_aligner(right_gap_score=0)

    alignment_arm_1 = aligner.align(arm_1, genome[:len(arm_1) * 2])[0]
    alignment_arm_2 = aligner.align(arm_2, genome[:len(arm_2) * 2])[0]

    if alignment_arm_1.score > alignment_arm_2.score:
        alignment = alignment_arm_1
        arm_len = len(arm_1)
    else:
        alignment = alignment_arm_2
        arm_len = len(arm_2)

    if alignment.score < arms_selfscore:
        raise RuntimeError(
            f'Something seems off since the alignment score of an arm ({alignment.score}) against the '
            f'genome is smaller than that of the arms against each other ({arms_selfscore})'
        )

    last_index_arm, last_index_genome = [
        (i, j)
        for i, j in zip(*alignment.indices)
        if i != -1 and j != -1
    ][-1]

    if last_index_arm + 1 != arm_len:
        raise RuntimeError('Last arm nucleotide did not match!')
    
    if last_index_genome == 0:
        raise RuntimeError('Something is off as the last reference index is 0!')
    
    matchlen_ref = last_index_genome + 1

    # get the flanking region
    flanking_region = genome[matchlen_ref:matchlen_ref+margin]

    return matchlen_ref, flanking_region


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--hairpin',
        required=True,
        help='Fasta file with two unequal hairpin arms.'
    )
    parser.add_argument(
        '--genome',
        required=True,
        help='Fasta file with genome.'
    )
    parser.add_argument(
        '--margin',
        required=True,
        type=int,
        help='Width of the flanking region.'
    )
    parser.add_argument(
        '--flanking_region_diff_tolerance',
        required=True,
        type=int,
        help='Differences tolerated between the flanking regions.'
    )
    parser.add_argument(
        '--out',
        required=True,
        help='Output directory.'
    )
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    fname_out_fasta = os.path.join(args.out, 'hairpin_region_to_refine.fasta')
    fname_out_results = os.path.join(args.out, 'hairpin_matching.json')

    genome = SeqIO.read(args.genome, 'fasta')
    arm_1, arm_2 = list(SeqIO.parse(args.hairpin, 'fasta'))

    arm_aligner = get_dna_aligner(right_gap_score=-1)
    hairpin_selfalignment = arm_aligner.align(arm_1, arm_2)[0]
    hairpin_selfscore = hairpin_selfalignment.score

    try:
        hp_length_l, flanking_l = align_end(
            arm_1.seq,
            arm_2.seq,
            genome.seq,
            hairpin_selfscore,
            args.margin
        )
        hp_length_r, flanking_r = align_end(
            arm_1.seq,
            arm_2.seq,
            genome.seq.reverse_complement(),
            hairpin_selfscore,
            args.margin
        )

        # compare the flanking regions
        flanking_regions_aligner = get_dna_aligner(right_gap_score=-1)
        flank_alignment = flanking_regions_aligner.align(flanking_r, flanking_l)[0]

        counts = flank_alignment.counts()
        cost = counts.gaps + counts.mismatches

        if cost > args.flanking_region_diff_tolerance:
            raise RuntimeError('Flanking regions differ too much')

        results = {
            'success': True,
            'genome_hairpin_len_left': int(hp_length_l),
            'genome_hairpin_len_right': int(hp_length_r),
            'margin_len': int(args.margin),
            'arm_1': str(arm_1.seq),
            'arm_2': str(arm_2.seq),
        }
        seq_out = (
            flanking_l.reverse_complement()
            + arm_1.seq.reverse_complement()
            + arm_2.seq
            + flanking_l
        )
        record_out = SeqRecord(
            id='hairpin_ref',
            description='',
            seq=seq_out
        )
        SeqIO.write(
            record_out,
            fname_out_fasta,
            'fasta'
        )

    except Exception as error:
        results = {
            'success': False,
            'genome_hairpin_len_left': 0,
            'genome_hairpin_len_right': 0,
            'margin_len': args.margin,
            'margin_diff_cost': int(cost),
            'arm_1': str(arm_1.seq),
            'arm_2': str(arm_2.seq),
            'err': str(error),
        }
        with open(fname_out_fasta, 'w') as f_out:
            f_out.write('')

    with open(fname_out_results, 'w') as f:
        json.dump(results, f, indent=4)


if __name__ == '__main__':
    main()
