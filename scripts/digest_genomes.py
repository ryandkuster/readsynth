#!/usr/bin/env python

import gzip
import pandas as pd
import re

from scripts.gzip_test import test_unicode


def main(args):
    """
    process fasta sequences as restricion enzyme fragments

    input args.g is fasta

    output 'raw_digest' file holds all the resulting fragments
    """
    if test_unicode(args.g):
        fasta = gzip.open(args.g, 'rt')
    else:
        fasta = open(args.g)

    begin, gen_ls, seq = 0, [], ''

    for line in fasta:
        if line.startswith('>') and seq:
            seq_ls = digest_seq(begin, seq, args)
            gen_ls.extend(seq_ls)
            begin += len(seq)
            seq = ''
            chr_name = line.rstrip()[1:].replace(' ', '_').replace(',', '')
        elif line.startswith('>'):
            chr_name = line.rstrip()[1:].replace(' ', '_').replace(',', '')
        else:
            seq += line.rstrip().upper()

    seq_ls = digest_seq(begin, seq, args)
    gen_ls.extend(seq_ls)
    begin += len(seq)
    fasta.close()

    df = pd.DataFrame(gen_ls,
                      columns=['seq', 'start', 'end', 'm1', 'm2', 'internal'])

    return df


def digest_seq(begin, seq, args):
    """
    for every chromosome (seq), find all RE recognition positions
    and preserve max bp ahead as a possible template (fragment)

    each item in seq_ls is [sequence, start, end]
    """
    seq_ls = []
    for motif1 in args.motif_len.keys():
        for idx in re.finditer('(?=' + motif1 + ')', seq):
            start = idx.start()
            fragment = seq[start: start + args.max]
            seq_ls.extend(digest_frag(fragment, motif1, begin+start, args))

    return seq_ls


def digest_frag(fragment, motif1, f_start, args):
    '''
    further search each RE starting point + args.max for more
    RE sites, return list of seq, start, end, m1, m2
    '''
    frag_ls = []

    for motif2 in args.motif_len.keys():
        for i, idx in enumerate(re.finditer('(?=' + motif2 + ')',
                                fragment[1:])):
            end = idx.start()+1
            internals = internal_sites(fragment[1:end+args.motif_len[motif2]-1],
                                       args.motif_len.keys())
            frag_ls.append([fragment[:end+args.motif_len[motif2]],
                            f_start,
                            f_start+end,
                            motif1,
                            motif2,
                            internals])

    return frag_ls


def internal_sites(subseq, all_motifs):
    internals = 0
    for motif3 in all_motifs:
        hits = [m.start() for m in re.finditer('(?=' + motif3 + ')', subseq)]
        internals += len(hits)

    return internals

