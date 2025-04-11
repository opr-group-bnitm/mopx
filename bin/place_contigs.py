#!/usr/bin/env python3

import subprocess
import tempfile
import os
import argparse
import itertools
import sys
import logging
import traceback

from Bio.Align import PairwiseAligner, Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import pysam


temp_rootdir = os.path.join(os.getcwd(), 'work')


def run_minimap2(fname_ref, fname_reads, fname_sam, secondary=False):
    """Run minimap2."""
    minimap2_cmd = [
        "minimap2",
        "-x", "map-ont",
        "-a",
        fname_ref,
        fname_reads
    ]
    if not secondary:
        minimap2_cmd.append("--secondary=no")
    minimap_out = subprocess.run(
        minimap2_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    if minimap_out.returncode != 0:
        raise RuntimeError(f'Minimap2 failed: {minimap_out.stderr}')
    with open(fname_sam, 'w') as f_sam:
        f_sam.write(minimap_out.stdout)
    sam = list(pysam.AlignmentFile(fname_sam, "r").fetch())
    return sam


def write_fasta(fname, sequences):
    """Write a list of sequence records to a fasta file."""
    SeqIO.write(sequences, fname, 'fasta')


def map_to_ref(
        ref: SeqRecord,
        contigs: SeqRecord | list[SeqRecord]
):
    with tempfile.TemporaryDirectory(dir=temp_rootdir) as temp_dir:
        fname_ref = os.path.join(temp_dir, 'ref.fasta')
        fname_contigs = os.path.join(temp_dir, 'contigs.fasta')
        fname_sam = os.path.join(temp_dir, 'mapped.sam')

        write_fasta(fname_ref, ref)
        write_fasta(fname_contigs, contigs)

        sam = run_minimap2(fname_ref, fname_contigs, fname_sam, secondary=False)
        return sam


def primary_alignments(sam_alignments):
    """Get all primary alignments from a list of sam alignments.

    This excludes secondary and supplementary alignments and unmapped reads
    """
    return [
        align
        for align in sam_alignments
        if not any((align.is_supplementary, align.is_secondary, align.is_unmapped))
    ]


def supplementary_alignments(sam_alignments):
    return [
        align
        for align in sam_alignments
        if align.is_supplementary and not align.is_unmapped
    ]


def sort_contigs(
        ref: SeqRecord,
        contigs: SeqRecord | list[SeqRecord]
):
    """Map contigs and sort them by distance to reference end."""
    sam = map_to_ref(ref, contigs)
    primaries = primary_alignments(sam)
    dist = [
        (align.query_name, min(align.reference_start, len(ref) - align.reference_end))
        for align in primaries
    ]
    sorted_dist = sorted(dist, key=lambda x: x[1], reverse=True)
    contigs_by_id = {c.id: c for c in contigs}
    sorted_contigs = [contigs_by_id[id] for id, _ in sorted_dist]

    return sorted_contigs


def place_contig(sam, ref):
    """Replace a the reference part where a sequence maps with the sequence.
    """
    candidates = primary_alignments(sam)
    if not candidates:
        return ref
    align = candidates[0]
    draft = ref.seq[:align.reference_start] + Seq(align.query_alignment_sequence) + ref.seq[align.reference_end:]
    return SeqRecord(Seq(draft), id='draft', description=f'contigs placed on {ref.id}')


def place_contigs_iteratively(ref, contigs):
    """Iteratively replace the reference sequence by the contigs.

    In each iteration, a contig is mapped to the current draft. The region
    where the contigs maps is then replace by the contig.
    """
    draft = ref
    for contig in contigs:
        sam = map_to_ref(draft, contig)
        draft = place_contig(sam, draft)
    return draft


def is_at_start(align, margin):
    return align.reference_start <= margin


def is_at_end(align, reflen, margin):
    return reflen - align.reference_end <= margin


def spans_start_to_end(align, reflen, margin):
    return is_at_start(align, margin) and is_at_end(align, reflen, margin)


def is_terminal(align, reflen, margin):
    return is_at_start(align, margin) or is_at_end(align, reflen, margin)


def get_overlapper_sam(candidates):
    if not candidates:
        return None
    alignment, _ = max(
        candidates,
        key=lambda aligns: min(aligns[0].query_alignment_length, aligns[1].query_alignment_length)
    )
    return alignment


def get_query_pos(align):
    """Get query positions with respect to the original query sequence
    in the original direction.
    """
    if align.cigartuples and align.cigartuples[0][0] == 5:
        hardclip_left = align.cigartuples[0][1]
    else:
        hardclip_left = 0
    if align.cigartuples and align.cigartuples[-1][0] == 5:
        hardclip_right = align.cigartuples[-1][1]
    else:
        hardclip_right = 0
    if align.is_reverse:
        qlen = len(align.query_sequence)
        qstart = qlen - align.query_alignment_end + hardclip_right
        qend = qlen - align.query_alignment_start + hardclip_right
    else:
        qstart = align.query_alignment_start + hardclip_left
        qend = align.query_alignment_end + hardclip_left

    return qstart, qend


def check_terminal_overlap_candidate(primary, supplementary, reflen, margin):
    """

    Returns
    -------
    is_match : bool
    side : 'start' or 'end'

    """
    # test match_direction
    primary_start, primary_end = get_query_pos(primary)
    supplementary_start, supplementary_end = get_query_pos(supplementary)
    same_direction = (primary.is_forward == supplementary.is_forward)

    suppl_is_prefix = supplementary_end <= primary_start
    suppl_is_suffix = primary_end <= supplementary_start

    if primary.is_reverse:
        suppl_is_prefix, suppl_is_suffix = (suppl_is_suffix, suppl_is_prefix)

    if suppl_is_suffix and is_at_end(primary, reflen, margin):
        # This is an end candidate
        # the alignments are expected to be like this
        # 
        # >>>*******<<<
        #      ------->
        #           <--
        #
        # or this
        # 
        # >>>*******<<<
        #      ------->
        # -->
        #
        # hence the supplementary sequence is at the end and in the opposite
        # direction or at the start and in the same direction (!= for bool make an xor)
        #
        matching_directions = (is_at_end(supplementary, reflen, margin) != same_direction)
        return 'end' if matching_directions else None
    elif suppl_is_prefix and is_at_start(primary, margin):
        # This is a start candidate
        # we expect alignments like this
        # 
        # >>>*******<<<
        # ------->
        # <--
        #
        # or this
        # 
        # >>>*******<<<
        # ------->
        #           -->
        #
        matching_directions = (is_at_start(supplementary, margin) != same_direction)
        return 'start' if matching_directions else None
    return None


class WrappedAlignment:

    def __init__(self, alignment: Alignment):
        self.alignment = alignment

    def __str__(self):
        if self.alignment is None:
            return 'Wrapped alignment of not existing alignment'
        return '\n'.join([
            'Wrapped alignment',
            f'Seqlen A={len(self.alignment.sequences[0])}',
            f'Seqlen B={len(self.alignment.sequences[1])}',
            f'Length={self.alignment_length}',
            f'Sequence identity={self.sequence_identity}',
            str(self.alignment)
        ])

    @property
    def alignment_length(self) -> int:
        """Length of the alignment without terminal gaps.
        """
        def terminal_gapcount(index_pairs):
            for gapcount, (i, j) in enumerate(index_pairs):
                if i != -1 and j != -1:
                    return gapcount
            return len(index_pairs)
        index_pairs = list(zip(*self.alignment.indices))
        start_gaps = terminal_gapcount(index_pairs)
        end_gaps = terminal_gapcount(reversed(index_pairs))
        return max(0, self.alignment.length - start_gaps - end_gaps)

    @property
    def sequence_identity(self) -> float:
        return self.alignment.counts().identities / self.alignment_length

    def reverse_complement(self):
        return type(self)(self.alignment.reverse_complement())


def find_turning_point(aligner, sequence):
    # We want to find an alignment like this:
    #
    # (i)
    # 0 1 2 3 x 4 5
    #     | |   | |
    #     5 4 x 3 2 1 0
    # (len - 1 - j)
    #
    # we want to separate at 'x'
    seq = Seq(sequence)
    alignment = aligner.align(seq, seq.reverse_complement())
    aln = alignment[0]
    aligned_pairs = [
        (i, len(sequence) - j - 1)
        for (start1, end1), (start2, end2) in zip(aln.aligned[0], aln.aligned[1])
        for i, j in zip(range(start1, end1), range(start2, end2))
        if i >= len(sequence) - j - 1
    ]
    closest_pair = min(aligned_pairs, key=lambda pair: pair[0] - pair[1])
    separation_pos = closest_pair[0]
    return WrappedAlignment(aln), separation_pos


class TerminusLoopResults:
    
    def __init__(
            self,
            draft: Seq,
            contig_name_5prime: str,
            valid_5prime: bool,
            alignment_5prime: WrappedAlignment,
            contig_name_3prime: str,
            valid_3prime: bool,
            alignment_3prime: WrappedAlignment,
    ):
        self.draft = draft

        self.contig_name_5prime = contig_name_5prime
        self.valid_5prime = valid_5prime
        self.alignment_5prime = alignment_5prime

        self.contig_name_3prime = contig_name_3prime
        self.valid_3prime = valid_3prime
        self.alignment_3prime = alignment_3prime

    def __str__(self):
        def info_str(terminus, contig_name, is_valid, align):
            return '\n'.join([
                f"Terminus={terminus}'",
                f"Contig={contig_name}",
                f"is_valid={is_valid}",
                f"{align}"
            ])
        info_string = '\n'.join([
            "Terminus Looping",
            info_str(
                "5'",
                self.contig_name_5prime,
                self.valid_5prime,
                self.alignment_5prime
            ),
            info_str(
                "3'",
                self.contig_name_3prime,
                self.valid_3prime,
                self.alignment_3prime
            )
        ])
        return info_string

    @property
    def found_start_overlap(self):
        return self.valid_5prime

    @property
    def found_end_overlap(self):
        return self.valid_3prime


class TerminusLooper:

    def __init__(
            self,
            align_margin: int = 10_000,
            max_dist_to_terminus: int = 200,
            min_seq_identity: float = 0.9,
            min_align_len: int = 200,
            align_match: float = 1.0,
            align_mismatch: float = -1.0,
            align_gap_open: float = -1.0,
            align_gap_extend: float = -1.0
    ):
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = align_match
        self.aligner.mismatch_score = align_mismatch
        self.aligner.open_gap_score = align_gap_open
        self.aligner.extend_gap_score = align_gap_extend
        self.aligner.query_end_gap_score = 0
        self.aligner.target_end_gap_score = 0

        self.align_margin = align_margin
        self.max_dist_to_terminus = max_dist_to_terminus
        self.min_seq_identity = min_seq_identity
        self.min_align_len = min_align_len

    def overlap_ends(
            self,
            ref: SeqRecord,
            contigs: list[SeqRecord]
    ):
        reflen = len(ref.seq)
        sam_alignments = map_to_ref(ref, contigs)
        primaries = [
            align
            for align in primary_alignments(sam_alignments)
            if is_terminal(align, reflen, self.max_dist_to_terminus)
        ]
        supplementaries = [
            align
            for align in supplementary_alignments(sam_alignments)
            if is_terminal(align, reflen, self.max_dist_to_terminus)
        ]
        candidate_pairs = [
            (primary, suppl)
            for primary, suppl in itertools.product(primaries, supplementaries)
            if primary.query_name == suppl.query_name 
        ]
        candidate_checks = [
            check_terminal_overlap_candidate(primary, supplementary, reflen, self.max_dist_to_terminus)
            for primary, supplementary in candidate_pairs
        ]
        start_candidates = [
            (primary, suppl)
            for (primary, suppl), start_or_end in zip(candidate_pairs, candidate_checks)
            if start_or_end == 'start'
        ]
        end_candidates = [
            (primary, suppl)
            for (primary, suppl), start_or_end in zip(candidate_pairs, candidate_checks)
            if start_or_end == 'end'
        ]
        start_align = get_overlapper_sam(start_candidates)
        end_align = get_overlapper_sam(end_candidates)

        valid_loopalign_start = False
        valid_loopalign_end = False
        loopalign_start = None
        loopalign_end = None

        if start_align is not None:
            #
            #    ------------> (reference)
            # ------->         (query)
            #    |
            #    query_alignment_start
            #
            start = max(0, start_align.query_alignment_start - self.align_margin)
            end = start_align.query_alignment_start + self.align_margin

            to_align = start_align.query_sequence[start:end]
            loopalign_start, separation_start = find_turning_point(self.aligner, to_align)
            start_5prime = start + separation_start
            end_5prime = start_align.query_alignment_end
    
            valid_loopalign_start = (
                start_5prime < end_5prime
                and loopalign_start.sequence_identity > self.min_seq_identity
                and loopalign_start.alignment_length > self.min_align_len
            )

        if end_align is not None:
            #
            # ------------>     (reference)
            #         ------->  (query)
            #
            start = max(0, end_align.query_alignment_end - self.align_margin)
            end = end_align.query_alignment_end + self.align_margin

            to_align = end_align.query_sequence[start:end]
            loopalign_end, separation_end = find_turning_point(self.aligner, to_align)
            start_3prime = end_align.query_alignment_start
            end_3prime = start + separation_end

            valid_loopalign_end = (
                start_3prime < end_3prime
                and loopalign_end.sequence_identity > self.min_seq_identity
                and loopalign_end.alignment_length > self.min_align_len
            )

        one_contig_spans_all = (
            valid_loopalign_start
            and valid_loopalign_end
            and start_align.query_name == end_align.query_name
        )
        if one_contig_spans_all:
            draft = ref.seq[start_5prime:end_3prime]
        else:
            prefix = ''
            start_mid = 0
            end_mid = len(ref.seq)
            suffix = ''
            if valid_loopalign_start:
                prefix = start_align.query_sequence[start_5prime:end_5prime]
                start_mid = start_align.reference_end
            if valid_loopalign_end:
                suffix = end_align.query_sequence[start_3prime:end_3prime]
                end_mid = end_align.reference_start
            draft = Seq(prefix + str(ref.seq)[start_mid:end_mid] + suffix)
            
        return TerminusLoopResults(
            draft,
            'None' if start_align is None else start_align.query_name,
            valid_loopalign_start,
            loopalign_start,
            'None' if end_align is None else end_align.query_name,
            valid_loopalign_end,
            loopalign_end,
        )


class TransferResults:

    def __init__(
            self,
            draft: Seq,
            align: WrappedAlignment = None,
            len_overhang: int = 0,
            start_to_end: bool = True,
    ):
        self.draft = draft
        self.align = align
        self.start_to_end = start_to_end
        self.len_overhang = len_overhang

    def __str__(self):
        direction_str = "5' to 3'" if self.start_to_end else "3' to 5'"
        info_str = '\n'.join([
            "Terminus transfer",
            f"Direction={direction_str}",
            f"Length overhang={self.len_overhang}",
            f"{self.align}"
        ])
        return info_str

    def reverse_complement(self):
        return type(self)(
            self.draft.reverse_complement(),
            self.align.reverse_complement() if self.align is not None else None,
            self.len_overhang,
            not self.start_to_end
        )


class TerminusTransfer:

    def __init__(
            self,
            align_margin_ref: int = 1_000,
            align_margin_query: int = 1_000,
            max_end_margin: int = 50,
            min_seq_identity: float = 0.8,
            min_align_len: int = 500,
            align_match: float = 1.0,
            align_mismatch: float = -1.0,
            align_gap_open: float = -1.0,
            align_gap_extend: float = -1.0
    ):
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'local'
        self.aligner.match_score = align_match
        self.aligner.mismatch_score = align_mismatch
        self.aligner.open_gap_score = align_gap_open
        self.aligner.extend_gap_score = align_gap_extend

        self.align_margin_ref = align_margin_ref
        self.align_margin_query = align_margin_query
        self.max_end_margin = max_end_margin
        self.min_seq_identity = min_seq_identity
        self.min_align_len = min_align_len

    def start_to_end(self, draft: Seq):
        """
        """
        # the terminal repeat region is transferred from the query
        # to the ref
        start_shift_ref = max(0, len(draft)-self.align_margin_ref)
        ref = draft[start_shift_ref:]
        query = draft[:self.align_margin_query].reverse_complement()
        if not ref or not query:
            return None
        alignments = self.aligner.align(ref, query)

        # check if there is a good alignment
        align = None
        for al in alignments:
            wrapped = WrappedAlignment(al)
            end_ref = wrapped.alignment.indices[0][-1]
            is_valid = (
                end_ref >= len(ref) - self.max_end_margin
                and wrapped.alignment_length >= self.min_align_len
                and wrapped.sequence_identity >= self.min_seq_identity
            )
            if is_valid:
                align = wrapped
                break

        if align is None:
            return TransferResults(draft)

        end_query = al.indices[1][-1]
        overhang = query[end_query+1:]
        draft_out = draft[:start_shift_ref + end_ref + 1] + overhang

        return TransferResults(
            draft_out,
            align=align,
            len_overhang=len(overhang)
        )

    def end_to_start(self, draft: Seq):
        transfer_rev = self.start_to_end(draft.reverse_complement())
        return transfer_rev.reverse_complement()


def place_terminal_repeats(
        ref,
        contigs,
        terminus_overlapper: TerminusLooper,
        terminus_transfer: TerminusTransfer,
    ):
    contig_overlaps = terminus_overlapper.overlap_ends(ref, contigs)

    if contig_overlaps.found_start_overlap and not contig_overlaps.found_end_overlap:
        transfer = terminus_transfer.start_to_end(contig_overlaps.draft)
    elif contig_overlaps.found_end_overlap and not contig_overlaps.found_start_overlap:
        transfer = terminus_transfer.end_to_start(contig_overlaps.draft)
    else:
        transfer = TransferResults(contig_overlaps.draft, None)

    return contig_overlaps, transfer


def build_draft(
        fname_ref: str,
        fname_contigs: str,
        fname_out: str,
        draft_id: str,
        terminus_looper: TerminusLooper,
        terminus_transfer: TerminusTransfer,
        logger
):
    """Determine the order and then place the contigs.
    """
    ref = next(SeqIO.parse(fname_ref, 'fasta'))
    contigs = list(SeqIO.parse(fname_contigs, 'fasta'))

    contigs_sorted = sort_contigs(ref, contigs)

    logger.info(f'Loaded {len(contigs_sorted)} contigs')

    contigs_placed_draft = place_contigs_iteratively(ref, contigs_sorted)

    logger.info('Contig placement done')

    tr_loop, tr_transfer = place_terminal_repeats(
        contigs_placed_draft,
        contigs,
        terminus_looper,
        terminus_transfer
    )
    final_draft = tr_transfer.draft

    logger.info(tr_loop)
    logger.info(tr_transfer)

    description=f'reference={ref.id}:n_contigs={len(contigs_sorted)}'
    draft_record = SeqRecord(Seq(final_draft), id=draft_id, description=description)
    write_fasta(fname_out, draft_record)


def logger_setup(fname=None):
    logger = logging.getLogger('logger')
    logger.setLevel(logging.DEBUG)
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)

    if fname:
        file_handler = logging.FileHandler(fname)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref', '-r', required=True)
    parser.add_argument('--contigs', '-c', required=True)
    parser.add_argument('--outfile', '-o', required=True)
    parser.add_argument('--log', '-l')
    parser.add_argument('--terminus_looper_align_margin', default=10_000, type=int)
    parser.add_argument('--terminus_looper_max_dist_to_terminus', default=200, type=int)
    parser.add_argument('--terminus_looper_min_seq_identity', default=0.8, type=float)
    parser.add_argument('--terminus_looper_min_align_len', default=200, type=int)
    parser.add_argument('--terminus_looper_align_match', default=1.0, type=float)
    parser.add_argument('--terminus_looper_align_mismatch', default=-1.0, type=float)
    parser.add_argument('--terminus_looper_align_gap_open', default=-1.0, type=float)
    parser.add_argument('--terminus_looper_align_gap_extend', default=-1.0, type=float)
    parser.add_argument('--terminus_transfer_align_margin_ref', default=2000, type=int)
    parser.add_argument('--terminus_transfer_align_margin_query', default=2000, type=int)
    parser.add_argument('--terminus_transfer_max_end_margin', default=100, type=int)
    parser.add_argument('--terminus_transfer_min_seq_identity', default=0.8, type=float)
    parser.add_argument('--terminus_transfer_min_align_len', default=300, type=int)
    parser.add_argument('--terminus_transfer_align_match', default=1.0, type=float)
    parser.add_argument('--terminus_transfer_align_mismatch', default=-1.0, type=float)
    parser.add_argument('--terminus_transfer_align_gap_open', default=-1.0, type=float)
    parser.add_argument('--terminus_transfer_align_gap_extend', default=-1.0, type=float)
    parser.add_argument('--work')
    parser.add_argument('--draftid', default='draft')
    args = parser.parse_args()

    logger = logger_setup(args.log)

    #if args.work:
    #    temp_rootdir = args.work

    os.makedirs(temp_rootdir, exist_ok=True)

    cmd = ' '.join(sys.argv)
    logger.info(f'Command: {cmd}')

    logger.info('Options:')
    for option, value in vars(args).items():
        logger.info(f'{option}: {value}')

    try:
        terminus_looper = TerminusLooper(
            align_margin=args.terminus_looper_align_margin,
            max_dist_to_terminus=args.terminus_looper_max_dist_to_terminus,
            min_seq_identity=args.terminus_looper_min_seq_identity,
            min_align_len=args.terminus_looper_min_align_len,
            align_match=args.terminus_looper_align_match,
            align_mismatch=args.terminus_looper_align_mismatch,
            align_gap_open=args.terminus_looper_align_gap_open,
            align_gap_extend=args.terminus_looper_align_gap_extend,
        )
        terminus_transfer = TerminusTransfer(
            align_margin_ref=args.terminus_transfer_align_margin_ref,
            align_margin_query=args.terminus_transfer_align_margin_query,
            max_end_margin=args.terminus_transfer_max_end_margin,
            min_seq_identity=args.terminus_transfer_min_seq_identity,
            min_align_len=args.terminus_transfer_min_align_len,
            align_match=args.terminus_transfer_align_match,
            align_mismatch=args.terminus_transfer_align_mismatch,
            align_gap_open=args.terminus_transfer_align_gap_open,
            align_gap_extend=args.terminus_transfer_align_gap_extend,
        )
        build_draft(
            args.ref,
            args.contigs,
            args.outfile,
            args.draftid,
            terminus_looper,
            terminus_transfer,
            logger
        )
    except Exception as err:
        logger.error(f'Failed with exception\n{err}')
        stack_trace = traceback.format_exc()
        logger.error(stack_trace)


if __name__ == '__main__':
    main()
