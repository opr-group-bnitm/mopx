#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pysam


def seq_str(seqid, pos, seq):
     seq_len = len(seq.replace('-', ''))
     start = pos
     end = max(pos, pos + seq_len - 1)
     return f'{seqid:10s} ({start:6d}-{end:6d}): {seq.lower()}'


def sequence_strings(refname, queryname, ref_pos, query_pos, ref, query):
    if not ref:
         ref = '-' * len(query)
    elif not query:
         query = '-' * len(ref)
    return seq_str(refname, ref_pos, ref), seq_str(queryname, query_pos, query)


def nosentinel(pos, ref, query):
    if not len(ref) == len(query) == 1:
        if ref[0] == query[0]:
            ref = ref[1:]
            query = query[1:]
            pos += 1
        elif ref[-1] == query[-1]:
            ref = ref[:-1]
            query = query[:-1]
    return pos, ref, query


def freq_str(variant_row):
    reads_all = variant_row['reads_all']
    a = variant_row['A']
    c = variant_row['C']
    g = variant_row['G']
    t = variant_row['T']
    gap = variant_row['deletions']
    fa = a / reads_all
    fc = c / reads_all
    fg = g / reads_all
    ft = t / reads_all
    fgap = gap / reads_all
    return (
        f'Counts:      A={a:4d}|C={c:4d}|G={g:4d}|T={t:4d}|del={gap:4d}\n'
        f'Frequencies: A={fa:1.2f}|C={fc:1.2f}|G={fg:1.2f}|T={ft:1.2f}|del={fgap:1.2f}'
    )


def get_features_at(gb_record, positions):
    return [
        feat
        for feat in gb_record.features
        if any(pos in feat.location for pos in positions)
        if feat.type not in {'source'}
    ]


def gb_str(feat):
    tags = '_'.join(
        feat.qualifiers.get(
            'locus_tag',
            feat.qualifiers.get('gene', ['UNKNOWN'])
        )
    )
    return f'{feat.type}\t{feat.location}\t{tags}'
    # type: gene
    # location: [194331:196104](+)
    # qualifiers:
    # Key: locus_tag, Value: ['MPXV-PCH-176']


def report(variants, genbank):
    last_pos_ref = -100
    last_pos_query = -100
    sum_ins = 0
    sum_del = 0
    for i, variant in enumerate(variants):

        pos_ref, ref, query = nosentinel(variant.pos, variant.ref, variant.alts[0])
        pos_query = pos_ref + sum_ins - sum_del

        if i > 0:
            mindist = min(pos_ref - last_pos_ref, pos_query - last_pos_query)
            print(f'\nDistance {mindist}\n')

        print(f'Variant {i}\n')
        print(f'len={max(len(ref), len(query))}')
        seqstr_ref, seqstr_query = sequence_strings(
            variant.contig,
            'Query',
            pos_ref,
            pos_query,
            ref,
            query,
        )
        print(seqstr_ref)
        print(seqstr_query)
        print('')

        basecounts = variant.info.get('BASECOUNTS')
        if basecounts:
            variant_list = [
                count.split(':')
                for count in basecounts.split('|')
            ]
            variant_row = {k: int(val) for k, val in variant_list}
            print(freq_str(variant_row))

        ref_range = range(pos_ref, pos_ref + min(1, len(ref)))
        for entry in get_features_at(genbank, ref_range):
            print(gb_str(entry))

        print('')

        sum_ins += len(query) 
        sum_del += len(ref)

        last_pos_ref = pos_ref + len(ref)
        last_pos_query = pos_query + len(query)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--variants', required=True, help='Input variant file in .vcf format.')
    parser.add_argument('--gb', required=True, help='Annotation file matching the reference.')
    args = parser.parse_args()

    variants = list(pysam.VariantFile(args.variants, 'r'))
    genbank = SeqIO.read(args.gb, 'genbank')
    report(variants, genbank)


if __name__ == "__main__":
    main()
