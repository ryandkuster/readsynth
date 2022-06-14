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
    """
    seq_ls = []

    for motif1 in args.motif_len.keys():
        for idx in re.finditer('(?=' + motif1 + ')', seq):
            start = idx.start()
            end = start + args.motif_len[motif1]
            fragment = seq[start: end]
            seq_ls.append([fragment,
                           begin+start,
                           begin+end,
                           motif1,
                           motif1,
                           0])

    return seq_ls

