#!/usr/bin/env python

import gzip
import pandas as pd
import re

from scripts.gzip_test import test_unicode


def main(motif_dt, args):
    """
    process fasta sequences as restricion enzyme fragments

    input args.genome is fasta

    output 'raw_digest' file holds all the resulting fragments
    """

    if test_unicode(args.genome):
        fasta = gzip.open(args.genome, 'rt')
    else:
        fasta = open(args.genome)

    begin, gen_ls, seq = 0, [], ''

    for line in fasta:
        if line.startswith('>') and seq:
            seq_ls = digest_seq(begin, seq, motif_dt, args.max)
            gen_ls.extend(seq_ls)
            begin += len(seq)
            seq = ''
            chr_name = line.rstrip()[1:].replace(' ', '_').replace(',','')
        elif line.startswith('>'):
            chr_name = line.rstrip()[1:].replace(' ', '_').replace(',','')
        else:
            seq += line.rstrip().upper()

    seq_ls = digest_seq(begin, seq, motif_dt, args.max)
    gen_ls.extend(seq_ls)
    begin += len(seq)
    fasta.close()

    df = pd.DataFrame(gen_ls, columns=['seq', 'start', 'end', 'm1', 'm2', 'internal'])

    return df


def digest_seq(begin, seq, motif_dt, frag_len):
    """
    for every chromosome (seq), find all RE recognition positions
    and preserve frag_len bp ahead as a possible template (fragment)

    each item in seq_ls is [sequence, start, end]
    """
    seq_ls = []
    for motif1 in motif_dt.keys():
        for idx in re.finditer('(?=' + motif1 + ')', seq):
            start = idx.start()
            fragment = seq[start: start + frag_len]
            seq_ls.extend(digest_frag(fragment, motif_dt, motif1, begin+start))

    return seq_ls


def digest_frag(fragment, motif_dt, motif1, f_start):
    '''
    further search each RE starting point + frag_len for more
    RE sites, return list of seq, start, end, m1, m2
    '''
    frag_ls  = []

    for motif2 in motif_dt.keys():
        for i, idx in enumerate(re.finditer('(?=' + motif2 + ')', fragment[1:])):
            end = idx.start()+1
            internals = internal_sites(fragment[1:end+len(motif2)-1], motif_dt.keys())
            frag_ls.append([fragment[:end+len(motif2)],
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



