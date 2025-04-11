#!/usr/bin/env python3

import argparse
import os
from collections import Counter
import pandas as pd

import pysam
import pysamstats

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_basecounts_pysamstats(fname_bam, fname_ref):
    return [
        {
            "position": record["pos"],
            "A": record["A"],
            "C": record["C"],
            "G": record["G"],
            "T": record["T"],
            "deletions": record["deletions"],
            "reads_all": record["reads_all"]
        }
        for record in pysamstats.stat_variation(pysam.AlignmentFile(fname_bam, "rb"), fafile=fname_ref)
    ]


def oriented_seq(seq: str, forward: bool) -> str:
    if forward:
        return seq
    return str(Seq(seq).reverse_complement())


def get_basecounts(fname_bam, refseq):
    basecounters = [
        Counter({
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0,
            'deletions': 0,
        })
        for _ in refseq
    ]
    with pysam.AlignmentFile(fname_bam, 'rb') as bam_file:
        secondary_ids = set(
            align.query_name
            for align in bam_file.fetch()
            if align.is_secondary
        )
        sequences = {
            align.query_name: oriented_seq(align.query_sequence, align.is_forward)
            for align in bam_file.fetch()
            if align.query_name in secondary_ids and not (align.is_secondary or align.is_supplementary)
        }
        # This includes supplementary and secondary alignments
        for alignment in bam_file.fetch():
            if alignment.is_unmapped:
                continue
            if alignment.is_secondary:
                query_seq = oriented_seq(sequences[alignment.query_name], alignment.is_forward)
            else:
                query_seq = alignment.query_sequence

            pos_seq = 0
            pos_ref = alignment.reference_start
            for operation, length in alignment.cigar:
                if operation == 0:
                    # Operation '0' is a match/mismatch
                    for i, c in enumerate(query_seq[pos_seq:pos_seq+length], pos_ref):
                        basecounters[i].update(c)
                    pos_seq += length
                    pos_ref += length
                elif operation == 1:
                    # Insertion to the reference
                    pos_seq += length
                elif operation == 2:
                    basecounters[pos_ref].update(['deletions'])
                    # Deletion from the reference
                    pos_ref += length
                elif operation == 3:
                    # N
                    pos_seq += length
                    pos_ref += length
                elif operation == 4:
                    # Soft clipping
                    pos_seq += length
    return [
        {
            'position': i,
            'reads_all': bc.total(),
            **bc,
        }
        for i, bc in enumerate(basecounters)
    ]


def main():
    parser = argparse.ArgumentParser(
        description='Compute nucleotide counts and mask positions with low coverage.'
    )
    parser.add_argument(
        '--bam', '-b',
        required=True,
        help='Input BAM file'
    )
    parser.add_argument(
        '--reference', '-r',
        required=True,
        help='Reference FASTA file'
    )
    parser.add_argument(
        '--outdir', '-o',
        required=True,
        help='Output directory.'
    )
    parser.add_argument(
        '--min-coverage',
        default=20,
        help='Positions with lower coverage below are masked with N.',
        type=int
    )
    args = parser.parse_args()

    fname_bam = args.bam
    fname_ref = args.reference
    outdir = args.outdir

    ref = SeqIO.read(fname_ref, 'fasta')

    basecounts = get_basecounts(fname_bam, ref.seq)

    if len(basecounts) != len(ref.seq):
        raise RuntimeError(
            f'Reference sequence length ({len(ref.seq)}) does not match '
            f'the number of frequency columns ({len(basecounts)})'
        )

    rows = []
    for counts, nt in zip(basecounts, ref.seq):
        if counts['reads_all'] < args.min_coverage:
            nt_out = 'N'
        else:
            most_common = Counter({nt: counts[nt] for nt in 'ACGT'}).most_common(1)[0][0]
            nt_out = nt if nt == most_common else 'N'
        rows.append({
            **counts,
            'nt_ref': nt,
            'nt_corrected': nt_out,
        })
    df = pd.DataFrame(rows)

    seq_corrected = ''.join(df['nt_corrected'])  # .strip('N')
    seq_record_corrected = SeqRecord(
        Seq(seq_corrected),
        id=f'{ref.id}',
        description=ref.description
    )
    fname_consensus = os.path.join(outdir, 'consensus_masked.fasta')
    SeqIO.write(seq_record_corrected, fname_consensus, 'fasta')
    print(f"Corrected and masked consensus written to {fname_consensus}")

    fname_table = os.path.join(outdir, 'nucleotide_counts.tsv')
    df.to_csv(fname_table, sep="\t", index=False)
    print(f"Base counts with totals saved to {fname_table}")


if __name__ == "__main__":
    main()
