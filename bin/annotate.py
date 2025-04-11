#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


class Orf:

    def __init__(self, start, end, atg_start, stop_found, prot_seq):
        self.start = start
        self.end = end
        self.atg_start = atg_start
        self.stop_found = stop_found
        self.prot_seq = prot_seq


class OrfFinder:

    def __init__(self, seq):
        self.seq = Seq(seq)

    def find_orf(self, start, end, forward):
        if forward:
             seq = self.seq[start:]
        else:
             seq = self.seq[:end+1].reverse_complement()

        # exclude partial triples
        seq = seq[0:3*(len(seq) // 3)]
        amino_acids = seq.translate(to_stop=True)

        # check if start and endpos were found
        dna_range = len(amino_acids) * 3
        starts_with_atg = str(seq[:3]).upper() == 'ATG'
        stop_codon_found = len(seq) >= dna_range + 3

        if forward:
            start_orf = start
            end_orf = start + dna_range + 3
        else:
            start_orf = end - dna_range - 2
            end_orf = end + 1

        return Orf(
            start_orf,
            end_orf,
            starts_with_atg,
            stop_codon_found,
            amino_acids
        )


def msa_seq_indices(msa_seq):
    return np.cumsum([int(s != '-') for s in msa_seq]) - 1


def get_genomes_and_mapping(alignment_file):
    alignment = AlignIO.read(alignment_file, "fasta")
    ref_seq = str(alignment[0].seq)
    query_seq = str(alignment[1].seq)
    mapping = {
        p1: p2
        for s1, s2, p1, p2 in zip(ref_seq, query_seq, msa_seq_indices(ref_seq), msa_seq_indices(query_seq))
        if s1 != '-' and s2 != '-'
    }
    return ref_seq.replace("-", ""), query_seq.replace("-", ""), mapping


def process_annotations(fname_alignment, fname_gb, outdir):

    _, query_seq, mapping = get_genomes_and_mapping(fname_alignment)
    query_orf = OrfFinder(query_seq)
    annotations = SeqIO.read(fname_gb, "genbank")
    alignment_file = os.path.join(outdir, "alignments.txt")

    os.makedirs(outdir, exist_ok=True)

    prot_aligner = PairwiseAligner()
    prot_aligner.mode = 'global'
    prot_aligner.substitution_matrix = None
    prot_aligner.open_gap_score = -3
    prot_aligner.extend_gap_score = -0.5

    transfers = []

    with open(alignment_file, "w") as out:

        for feature in annotations.features:

            # TODO: how to treat gene?

            if feature.type in ["CDS"]:

                start = feature.location.start
                end = feature.location.end
                strand = feature.location.strand
                locus_tag = feature.qualifiers.get(
                    'locus_tag',
                    feature.qualifiers.get('gene', ['UNKNOWN'])
                )[0]

                if strand not in [-1, 1]:
                    raise RuntimeError(f'Strand needs to be 1 or -1 but is {strand}')

                ref_protein_seq = feature.qualifiers.get("translation", [""])[0]
                query_start = mapping.get(start, None) if mapping.get(start) else None
                query_end = mapping.get(end-1, None) if mapping.get(start) else None

                case_results = {
                    'LocusTag': locus_tag,
                    'RefStart': start,
                    'RefEnd': end,
                    'Strand': strand,
                    'MappedStart': query_start,
                    'MappedEnd': query_end
                }

                has_anchor = (
                    (strand == 1 and query_start)
                    or (strand == -1 and query_end)
                )

                if has_anchor:

                    orf = query_orf.find_orf(query_start, query_end, strand==1)
                    prot_align = prot_aligner.align(ref_protein_seq, orf.prot_seq)[0]

                    out.write(f"Feature: {locus_tag}\n")
                    out.write("Protein Alignment:\n")
                    out.write(str(prot_align))

                    counts = prot_align.counts()

                    align_results = {
                        'QueryStart': orf.start,
                        'QueryEnd': orf.end,
                        'QueryStartATG': orf.atg_start,
                        'StopFound': orf.stop_found,
                        'mismatches': counts.mismatches,
                        'matches': counts.identities,
                        'gaps': counts.gaps,
                    }
                else:
                    align_results = {
                        'QueryStart': None,
                        'QueryEnd': None,
                        'QueryStartATG': None,
                        'StopFound': None,
                        'mismatches': None,
                        'matches': None,
                        'gaps': None,
                    }

                transfers.append(case_results | align_results)
    results_table = pd.DataFrame(transfers)
    results_table.to_csv(os.path.join(outdir, 'transfers.tsv'), sep="\t")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--alignment', required=True, help='FASTA file with alignment')
    parser.add_argument('--gb', required=True, help='Genbank file with annotations')
    parser.add_argument('--out', '-o', required=True, help='Directory to save output files')
    args =  parser.parse_args()
    process_annotations(args.alignment, args.gb, args.out)


if __name__ == "__main__":
    main()
