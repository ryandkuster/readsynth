#!/usr/bin/env python

import gzip
import pandas as pd
import re

from scripts.gzip_test import test_unicode


def main(args):
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
            seq_ls = digest_seq(begin, seq, args.motif_dt, args.max)
            gen_ls.extend(seq_ls)
            begin += len(seq)
            seq = ''
            chr_name = line.rstrip()[1:].replace(' ', '_').replace(',','')
        elif line.startswith('>'):
            chr_name = line.rstrip()[1:].replace(' ', '_').replace(',','')
        else:
            seq += line.rstrip().upper()

    seq_ls = digest_seq(begin, seq, args.motif_dt, args.max)
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
        mot_len = get_query_len(motif1)
        for idx in re.finditer('(?=' + motif1 + ')', seq):
            start = idx.start()
            end = start + mot_len
            fragment = seq[start: end]
            seq_ls.append([fragment,
                            begin+start,
                            begin+end,
                            motif1,
                            motif1,
                            0])

    return seq_ls


def get_query_len(query):
    mot_len, count = 0, True
    for i in query:
        if i == '[':
            count = False
        elif i == ']':
            count = True
            mot_len += 1
        else:
            if count is True:
                mot_len += 1

    return mot_len




