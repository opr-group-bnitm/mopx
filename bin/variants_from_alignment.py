#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import pysam


BASECOUNT_COLS = ['A', 'C', 'G', 'T', 'deletions', 'reads_all']


def read_pairwise_alignment(fname: str):
    sequences = list(SeqIO.parse(fname, 'fasta'))
    if len(sequences) != 2:
        raise RuntimeError(f'Found more than two sequences in file {fname}')
    ref, query = sequences
    return ref, query


def sequence_positions_from_aligned(seq: Seq):
    # gaps always have the position of the following nucleotide
    csum = np.cumsum([int(c != "-") for c in seq])
    for i, (s, c) in enumerate(zip(csum, seq)):
        if c != '-':
            csum[i] -= 1
    return csum


def add_sentinel_if_necessary(pos, ref, alt, reference_seq):
    '''Add sentinels if necessary to convert .tsv to VCF entries.'''
    # Add sentinel to InDels
    if ref == '' or alt == '':
        if pos > 0:
            sentinel = reference_seq[pos - 1]
            ref = sentinel + ref
            alt = sentinel + alt
            pos = pos - 1
        elif len(ref) < len(reference_seq):
            sentinel = reference_seq[len(ref)]
            ref = ref + sentinel
            alt = alt + sentinel
        else:
            ref = ref if ref else '-'
            alt = alt if alt else '-'

    return pos, ref.upper(), alt.upper()


def get_variants(ref: Seq, query: Seq):
    """Get variants from aligned sequences.
    """
    ref_positions = sequence_positions_from_aligned(ref)
    query_positions = sequence_positions_from_aligned(query)

    pairs = list(zip(
        ref_positions,
        query_positions,
        ref,
        query
    ))

    # concatenate gaps
    vars = []
    for pos_ref, pos_query, seq_ref, seq_query in pairs:
        if seq_ref == seq_query:
            continue
        sref = seq_ref.replace('-', '')
        squery = seq_query.replace('-', '')
        if not vars or (pos_ref > vars[-1]['pos_ref'] and pos_query > vars[-1]['pos_query']):
            vars.append({
                'pos_ref': pos_ref,
                'pos_query': pos_query,
                'ref': sref,
                'alt': squery,
            })
        else:
            vars[-1]['ref'] += sref
            vars[-1]['alt'] += squery

    # Add sentinels for indels
    reference_genome = ref.replace('-', '')
    vars_with_sentinels = [
        add_sentinel_if_necessary(
            variant['pos_ref'],
            variant['ref'],
            variant['alt'],
            reference_genome
        )
        for variant in vars
    ]

    return vars_with_sentinels


def read_basecounts(fname):
    return pd.read_csv(fname, sep='\t', index_col='position')[BASECOUNT_COLS]


def write_vcf(variants, basecounts, fname, reference_id):

    header=pysam.VariantHeader()
    header.add_line('##fileformat=VCFv4.2')
    header.add_line(f'##contig=<ID={reference_id}>')
    header.add_line(
        '##INFO=<ID=BASECOUNTS,Number=1,Type=String,'
        + 'Description="Base counts: '
        + ' | '.join(BASECOUNT_COLS)
        + '">'
    )
    # header.add_meta(
    #     key='FORMAT',
    #     items=[('ID', 'GT'), ('Number', '1'), ('Type', 'String'), ('Description', 'Genotype')]
    # )
    sum_ins = 0
    sum_del = 0
    with pysam.VariantFile(fname, 'w', header=header) as vcf_out:        
        for pos, ref, alt in variants:
            record = vcf_out.new_record(
                contig=reference_id,
                start=pos,
                stop=pos + len(ref),
                alleles=(ref, alt)
            )
            pos_query = pos + sum_ins - sum_del
            sum_del += len(ref)
            sum_ins += len(alt)            
            # Add base count info if it's a SNP
            if len(ref) == len(alt) == 1 and pos_query in basecounts.index:
                counts = basecounts.loc[pos_query]
                basecount_str = '|'.join([f'{k}:{counts[k]}' for k in BASECOUNT_COLS])
                record.info['BASECOUNTS'] = basecount_str
            vcf_out.write(record)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--align", required=True, help="Input alignment file in FASTA format.")
    parser.add_argument("--basecounts", required=True, help="Basecounts .tsv file for the query sequence.")
    parser.add_argument("--out", default='out.vcf', help="Output file in .vcf format.")
    args = parser.parse_args()

    ref, query = read_pairwise_alignment(args.align)

    variants = get_variants(ref.seq, query.seq)
    basecounts = read_basecounts(args.basecounts)

    write_vcf(variants, basecounts, args.out, ref.id)


if __name__ == "__main__":
    main()
