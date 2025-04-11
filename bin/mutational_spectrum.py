#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import pysam


def read_snps_table(fname_vcf, fname_ref):
    ref = SeqIO.read(fname_ref, 'fasta').seq
    with pysam.VariantFile(fname_vcf, 'r') as vcf_file:
        # Careful! Pysam writes 0 indexed (and automatically adds 1)
        # but reads 1 indexed (so you have to subtract 1)
        return pd.DataFrame([
            {
                "context_5'": (ref[snp.pos - 2] if snp.pos - 1 > 0 else '').lower(),
                'ref': snp.ref.lower(),
                'alt': snp.alts[0].lower(),
                "context_3'": (ref[snp.pos] if snp.pos < len(ref) else '').lower(),
            }
            for snp in vcf_file
            if len(snp.ref) == len(snp.alts[0]) == 1
        ])


def revcomp(seq: str):
    return str(Seq(seq).reverse_complement())


def lexical_max(ref: str, alt: str):
    return max((ref, alt), (revcomp(ref), revcomp(alt)))


def get_context_table(ref_col, query_col):
    df = pd.DataFrame({
        'ref': ref_col,
        'alt': query_col,
    })
    df[['ref', 'alt']] = df.apply(
        lambda row: pd.Series(lexical_max(row['ref'], row['alt'])), axis=1
    )
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--ref', required=True)
    parser.add_argument('--out', required=True, help='Directory for the output files.')
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    snps = read_snps_table(args.vcf, args.ref)

    triples = get_context_table(
        (snps["context_5'"] + snps['ref'] + snps["context_3'"]),
        (snps["context_5'"] + snps['alt'] + snps["context_3'"]),
    )
    varcounts_triples = triples.groupby(['ref', 'alt']).size().reset_index(name='count')
    varcounts_triples.to_csv(os.path.join(args.out, 'triples.tsv'), sep='\t')

    pairs = get_context_table(
        pd.concat([
            (snps["context_5'"] + snps["ref"]),
            (snps["ref"] + snps["context_3'"]),
        ]),
        pd.concat([
            (snps["context_5'"] + snps["alt"]),
            (snps["alt"] + snps["context_3'"]),
        ]),
    )
    varcounts_pairs = pairs.groupby(['ref', 'alt']).size().reset_index(name='count')
    varcounts_pairs.to_csv(os.path.join(args.out, 'pairs.tsv'), sep='\t')


if __name__ == '__main__':
    main()
