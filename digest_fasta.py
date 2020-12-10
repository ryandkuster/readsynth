#!/usr/bin/env python

import argparse
import matplotlib as plt
import pandas as pd
import re


def parse_user_input():
    parser = argparse.ArgumentParser(description='simulate RE digest on genome')

    parser.add_argument('-genome', type=str, required=True, metavar='',
            help='path to file genome')

    parser.add_argument('-o', type=str, required=True,  metavar='',
            help='path to store output')

    parser.add_argument('-m1', type=str, required=True, nargs='+', metavar='',
            help='space separated list of RE motifs (e.g., AluI = AG/CT, HindIII = A/AGCTT), SmlI = C/TYRAG')

    parser.add_argument('-l', type=int, required=False, metavar='',
            help='max read length after first cut (optional, defaults to 1000bp)')

    parser.add_argument('-complete', type=int, required=False, metavar='',
            help='use \'1\' for complete digestion of fragments (fragments will not contain internal RE sites)')

    args = parser.parse_args()

    return args


def iupac_motifs():
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

    for motif in args.m1:
        reg_motif = ''
        for char in motif.upper():
            reg_motif += iupac_dt[char]
        motif_dt[reg_motif] = motif.index('/')

    return motif_dt


def test_chrom(motif_dt, read_len):
    with open(args.genome) as fasta, open('chopped_dna.fastq', 'w') as outfile:
        for line in fasta:
            if line.startswith('>'):
                try:
                    seq_ls = digest_seq(seq, motif_dt, read_len)
                    second_digest(seq_ls, motif_dt, read_len)
                except UnboundLocalError:
                    pass
                seq = ''
            else:
                seq += line.rstrip().upper()
        seq_ls = digest_seq(seq, motif_dt, read_len)
        second_digest(seq_ls, motif_dt, read_len)


def digest_seq(seq, motif_dt, read_len):
    seq_ls = []
    for motif, offset in motif_dt.items():
        for idx in re.finditer(motif, seq):
            start = idx.start()+offset
            end = start + read_len
            seq_ls.append(seq[start:end])

    seq_ls = [reverse_comp(i) for i in seq_ls]
    return seq_ls
    #TODO keep only seqs that fall within a distribution
    #outfile.write(seq[start:end] + '\n')


def reverse_comp(seq):
    revc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    new = ''
    for base in reversed(seq):
        new += revc[base]
    return new


def second_digest(seq_ls, motif_dt, read_len):
    second_ls = []
    for seq in seq_ls:
        second_ls += digest_seq(seq, motif_dt, read_len)
    if args.complete == 1:
        second_ls = [complete_digest(i, motif_dt) for i in second_ls]
        second_ls = [i for i in second_ls if i]
    print(second_ls)


def complete_digest(seq, motif_dt):
    for motif in motif_dt:
        if motif in seq [1:-1]:
            return False
    return seq


def simulate_profile():
    pass


if __name__ == '__main__':
    args = parse_user_input()
    motif_dt = iupac_motifs()
    read_len = args.l if args.l else 1000
    test_chrom(motif_dt, read_len)

    #add mean and std dev for distribution of sizes
    #check if genome compressed
    #be sensitive to uppercase in genome
    #allow user input of possible barcodes (optional, from anemone)
    #allow user input of T-overhang or other motif
    #check error profile file (optional, from crinoid) use pandas
    #convert genome into readable object by '>' entry
    #search genome and write cut sites to output file
    #be sure that rev_comp sensitivity is considered
    #build in error-profile creation (optional)

