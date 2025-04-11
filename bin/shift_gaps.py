#!/usr/bin/env python3

import argparse
from Bio.Seq import Seq
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    seqs = list(SeqIO.parse(args.align, 'fasta'))
    if not len(seqs) == 2:
        raise RuntimeError(f'Expecting exactly two sequences but got {len(seqs)}')

    ref, query = seqs
    pairs = list(zip(ref.seq, query.seq))

    if any(nt_ref == nt_query == '-' for nt_ref, nt_query in pairs):
        raise RuntimeError('Gap in both sequences is not allowed!')

    # leave the 5' terminal gap where it is
    for i, (nt_ref, nt_query) in enumerate(pairs):
        if nt_ref != '-' and nt_query != '-':
            break

    j = -1
    while i < len(pairs):
        if j < i:
            j = i
        nt_ref, nt_query = pairs[i]
        if nt_ref == '-':
            # find the end of the insertion
            while j < len(pairs) and pairs[j][0] == '-':
                j += 1
            if j == len(pairs):
                # we are at the end, no more changes possible
                i = j
                break
            if nt_query == pairs[j][1] or nt_query == pairs[j][0]:
                # the gap placement is equivalent or it is changed to a match=> swap the ref part
                pairs[i], pairs[j] = (pairs[j][0], pairs[i][1]), (pairs[i][0], pairs[j][1])
                i += 1
            else:
                # no swap, we do not need to look at this gap any longer
                i = j
        elif nt_query == '-':
            # find the end of the deletion
            while j < len(pairs) and pairs[j][1] == '-':
                j += 1
            if j == len(pairs):
                # we are at the end, no more changes possible
                i = j
                break
            if nt_ref == pairs[j][0] or nt_query == pairs[j][1]:
                # the gap placement is equivalent or changed to a match => swap the query part
                pairs[i], pairs[j] = (pairs[i][0], pairs[j][1]), (pairs[j][0], pairs[i][1])
                i += 1
            else:
                # no swap, we do not need to look at this gap any longer
                # and we can go to its end
                i = j
        else:
            i += 1

    ref_shifted, query_shifted = zip(*pairs)
    ref.seq = Seq(''.join(ref_shifted))
    query.seq = Seq(''.join(query_shifted))
    SeqIO.write([ref, query], args.out, 'fasta')


if __name__ == '__main__':
    main()
