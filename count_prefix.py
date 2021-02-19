#!/usr/bin/env python

import argparse
import pandas as pd
import sys

prefix_dt = {}

def parse_user_input():
    parser = argparse.ArgumentParser(description='count frequency of first n bases of a newline separated file of sequences')

    parser.add_argument('-f', type=str, required=True, metavar='',
            help='path to file of sequences csv')

    parser.add_argument('-rev', type=str, required=False, metavar='',
            help='use 1 to consider only the reverse complement of reads')

    parser.add_argument('-n', type=int, required=True,  metavar='',
            help='number of positions to calculate')

    parser.add_argument('-min', type=int, required=True,  metavar='',
            help='minimum read length to consider')

    parser.add_argument('-max', type=int, required=True,  metavar='',
            help='maximum read length to consider')

    parser.add_argument('-p', type=float, required=False,  metavar='',
            help='percent of reads in defined size range at which to use for motif')

    args = parser.parse_args()

    return args


def reverse_comp(seq):
    revc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    new = ''
    for base in reversed(seq):
        new += revc[base]

    return new


def get_top_hits(df, prefix_dt):
    total = len(df)
    with open(args.f + '_oligos.txt', 'w') as o:
        for seq in reversed(sorted(prefix_dt.items(), key = lambda kv:(kv[1], kv[0]))):
            percent = (seq[1]/total)*100
            if args.p:
                if percent >= args.p:
                    o.write(seq[0] + '\n')
            else:
                o.write(seq[0] + '\n')


if __name__ == '__main__':
    args = parse_user_input()
    df = pd.read_csv(args.f, header=None)
    df = df[(df.iloc[:, 1] >= args.min) & (df.iloc[:, 1] <= args.max)]
    df = df[~df.iloc[:, 0].str.contains('N')]

    for seq in df.iloc[:, 0].values:
        if args.rev:
            seq = reverse_comp(seq)
        if seq[:args.n] in prefix_dt:
            prefix_dt[seq[:args.n]] += 1
        else:
            prefix_dt[seq[:args.n]] = 1

    get_top_hits(df, prefix_dt)
    print(sum(prefix_dt.values()))
    print(len(df))

