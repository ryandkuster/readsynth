#!/usr/bin/env python

import argparse
import gzip
import os
import pickle
import re

from scripts.gzip_test import test_unicode


def main(args):
    '''
    process fasta sequences as restricion enzyme fragments

    input args.genome is fasta

    output 'raw_digest' file holds all the resulting fragments
    '''

    args.m1, args.m2 = check_for_enzymes(args)
    args.motif_dt = {}
    args.motif_dt1 = iupac_motifs(args.m1)
    args.motif_dt.update(args.motif_dt1)
    args.motif_dt2 = iupac_motifs(args.m2)
    args.motif_dt.update(args.motif_dt2)

    if test_unicode(args.genome):
        fasta = gzip.open(args.genome, 'rt')
    else:
        fasta = open(args.genome)

    begin, seq = 0, ''
    seq_dt = {k: [] for k in args.motif_dt}

    for line in fasta:
        if line.startswith('>') and seq:
            seq_dt_chrom = digest_seq(begin, seq, args.motif_dt)
            seq_dt = {k: seq_dt[k] + v for k, v in seq_dt_chrom.items()}
            begin += len(seq)
            seq = ''
        elif line.startswith('>'):
            pass
        else:
            seq += line.rstrip().upper()

    seq_dt_chrom = digest_seq(begin, seq, args.motif_dt)
    seq_dt = {k: seq_dt[k] + v for k, v in seq_dt_chrom.items()}

    fasta.close()

    for k, v in seq_dt.items():
        print(f'{k} : {len(v)}')


def parse_user_input():
    parser = argparse.ArgumentParser(description='simulate RAD libary')

    parser.add_argument('-genome', type=str, required=True,
                        help='path to file genome')

    parser.add_argument('-o', type=str, required=True,
                        help='path to store output')

    parser.add_argument('-m1', type=str, required=True, nargs='+',
                        help='space separated list of RE motifs (e.g., AluI = AG/CT, HindIII = A/AGCTT, SmlI = C/TYRAG)')

    parser.add_argument('-m2', type=str, required=True, nargs='+',
                        help='space separated list of RE motifs (e.g., AluI = AG/CT, HindIII = A/AGCTT, SmlI = C/TYRAG)')

    args = parser.parse_args()

    return args


def check_for_enzymes(args):
    with open(os.path.join(os.path.dirname(__file__),
              'resources/type_ii_enzymes.pickle'), 'rb') as type_ii_file:
        re_dt = pickle.load(type_ii_file)

    m1 = [re_dt[i.lower()] if i.lower() in re_dt.keys() else i for i in args.m1]
    m2 = [re_dt[i.lower()] if i.lower() in re_dt.keys() else i for i in args.m2]

    return m1, m2


def iupac_motifs(arg_m):
    '''
    given a list of RE cut motifs, return a dictionary of regex
    compatible iupac redundancy codes as keys and cleaving site
    as the value
    '''
    motif_dt = {}
    iupac_dt = {'/': '',
                'A': 'A',
                'C': 'C',
                'G': 'G',
                'T': 'T',
                'R': '[AG]',
                'Y': '[CT]',
                'S': '[GC]',
                'W': '[AT]',
                'K': '[GT]',
                'M': '[AC]',
                'B': '[CGT]',
                'D': '[AGT]',
                'H': '[ACT]',
                'V': '[ACG]',
                'N': '[ACGT]'}

    for motif in arg_m:
        reg_motif = ''
        for char in motif.upper():
            reg_motif += iupac_dt[char]
        motif_dt[reg_motif] = motif.index('/')

    return motif_dt


def digest_seq(begin, seq, motif_dt):
    '''
    for every chromosome (seq), find all RE recognition positions
    and preserve frag_len bp ahead as a possible template (fragment)

    each item in seq_ls is [sequence, start, end]
    '''
    seq_dt_chrom = {k: [] for k in args.motif_dt}
    for motif1 in motif_dt.keys():
        for idx in re.finditer('(?=' + motif1 + ')', seq):
            start = idx.start()
            seq_dt_chrom[motif1].append(begin+start)

    return seq_dt_chrom


if __name__ == '__main__':
    args = parse_user_input()
    main(args)
