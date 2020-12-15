#!/usr/bin/env python

import argparse
import matplotlib as plt
import numpy as np
import pandas as pd
import random
import re

import sys #TODO delete this is not used

def parse_user_input():
    parser = argparse.ArgumentParser(description='simulate RE digest on genome')

    parser.add_argument('-genome', type=str, required=True, metavar='',
            help='path to file genome')

    parser.add_argument('-o', type=str, required=True,  metavar='',
            help='path to store output')

    parser.add_argument('-m1', type=str, required=True, nargs='+', metavar='',
            help='space separated list of RE motifs (e.g., AluI = AG/CT, HindIII = A/AGCTT), SmlI = C/TYRAG')

    parser.add_argument('-l', type=int, required=False, metavar='',
            help='desired read length of final simulated reads (defaults to 250 or given q1/q2 profiles)')

    parser.add_argument('-n', type=int, required=True, metavar='',
            help='number of reads to simulate (currently per chromosome)')

    parser.add_argument('-a1s', type=int, required=True, metavar='',
            help='index (beginning with 0) of read SBS start site within R1 adapter')

    parser.add_argument('-a2s', type=int, required=True, metavar='',
            help='index (beginning with 0) of read SBS start site within R2 adapter')

    parser.add_argument('-f', type=int, required=False, metavar='',
            help='max fragment length after first cut (optional, defaults to 1000bp)')

    parser.add_argument('-complete', type=int, required=False, metavar='',
            help='use \'1\' for complete digestion of fragments (fragments will not contain internal RE sites)')

    parser.add_argument('-a1', type=str, required=True, metavar='',
            help='file containing tab/space-separated adapters and identifiers that attach 5\' to read')

    parser.add_argument('-a2', type=str, required=True, metavar='',
            help='file containing tab/space-separated adapters and identifiers that attach 3\' to read')

    parser.add_argument('-q1', type=str, required=True, metavar='',
            help='file containing R1 q scores in csv format (see ngsComposer tool crinoid)')

    parser.add_argument('-q2', type=str, required=True, metavar='',
            help='file containing R2 q scores in csv format (see ngsComposer tool crinoid)')

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


def get_adapters(arg):
    adapters_dt = {}
    with open(arg) as f:
        for line in f:
            adapter, id = line.rstrip().split()
            adapters_dt[adapter] = id

    return adapters_dt


def get_qscores(arg):
    print('getting q score profile')
    scores = list('!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJK')
    nums = [i for i in range(43)]
    df = pd.read_csv(arg, names=scores)
    df = df.loc[:, df.any()]
    scores_ls = df.columns.tolist()
    scores_dt = {v: 10**(-k/10) for k,v in zip(nums, scores) if v in scores_ls}

    prob_mx = []
    for idx, row in df.iterrows():
        q_ls = [i/sum(row) for i in row]
        prob_mx.append(q_ls)
    #TODO remove extra length in profile or autodetect length
    return prob_mx, scores_dt, scores_ls


def test_chrom(motif_dt, frag_len):
    print('converting fasta')
    with open(args.genome) as fasta,\
         open('simulated_R1.fastq', 'w') as r1,\
         open('simulated_R2.fastq', 'w') as r2:

        for line in fasta:
            if line.startswith('>'):
                chr_name = line.rstrip()[1:].replace(' ', '_')
                try:
                    seq_ls = digest_seq(seq, motif_dt, frag_len)
                    seq_ls = second_digest(seq_ls, motif_dt, frag_len)
                    #TODO add steps up here
                except UnboundLocalError:
                    pass
                seq = ''
            else:
                seq += line.rstrip().upper()
        seq_ls = digest_seq(seq, motif_dt, frag_len)
        seq_ls = second_digest(seq_ls, motif_dt, frag_len)
        seq_ls, len_ls = simulate_adapters(seq_ls)
        seq_ls = simulate_length(seq_ls, len_ls) #TODO digest_seq, second_digest, simulate_length, etc. need organization
        read_writer(seq_ls, r1, r2, chr_name)


def digest_seq(seq, motif_dt, frag_len):
    seq_ls = []
    for motif, offset in motif_dt.items():
        for idx in re.finditer(motif, seq):
            start = idx.start()+offset
            end = start + frag_len
            seq_ls.append(seq[start:end])

    seq_ls = [reverse_comp(i) for i in seq_ls]

    return seq_ls


def reverse_comp(seq):
    revc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    new = ''
    for base in reversed(seq):
        new += revc[base]

    return new


def second_digest(seq_ls, motif_dt, frag_len):
    print('performing second digest')
    second_ls = []
    for seq in seq_ls:
        second_ls += digest_seq(seq, motif_dt, frag_len)
    if args.complete == 1:
        print('including only complete digests')
        second_ls = [complete_digest(i, motif_dt) for i in second_ls]
        second_ls = [i for i in second_ls if i]

    return second_ls


def complete_digest(seq, motif_dt):
    for motif in motif_dt:
        if motif in seq [1:-1]:
            return False

    return seq


def simulate_adapters(seq_ls):
    print('adding random adapters')
    adapt_ls, len_ls = [], []
    for seq in seq_ls:
        r1_dapt, r1_id = random.choice(list(args.a1.items()))
        r2_dapt, r2_id = random.choice(list(args.a2.items()))
        ligated_seq = r1_dapt + seq + r2_dapt
        adapt_ls.append([ligated_seq, r1_id, r2_id])
        len_ls.append(len(seq))

    return adapt_ls, len_ls


def simulate_length(seq_ls, len_ls):
    df = pd.DataFrame(seq_ls, columns = ['sequence', 'r1_id', 'r2_id'])
    df['length'] = [len(i[0]) for i in seq_ls]
    df['fragment_length'] = len_ls
    df.sort_values(['length'], ascending=[True], inplace=True)
    df.reset_index(inplace=True, drop=True)

    print(df)

    df.to_csv('raw_digest_df.csv', index=None) #TODO add option to write raw digests to file
    read_count = len(df) #TODO keep this in case you want to multiply to get coverage
    mean = 400 #TODO add argument MEAN
    sd = 100 #TODO add argument SD
    total_reads = args.n #TODO add argument COVERAGE

    draw_ls = np.random.normal(loc=mean,scale=sd,size=total_reads)
    draw_ls = [round(i) for i in draw_ls]
    draw_dt = {}

    for i in range(max(min(df['length']), min(draw_ls)), min(max(df['length']), max(draw_ls))+1):
        draw_dt[i] = draw_ls.count(i)

    sampled_df = pd.DataFrame(columns=['sequence', 'r1_id', 'r2_id', 'length', 'fragment_length'])
    counts = []

    for length, draws in draw_dt.items():
        tmp_df =  df.loc[df['length'] == length]
        if len(tmp_df) == 0:
            continue
        indices = [i for i in range(len(tmp_df))]
        sampled_idx = random.choices(indices, k=draws)
        counts += [sampled_idx.count(idx) for idx in indices]
        sampled_df = pd.concat([sampled_df, tmp_df])

    sampled_df['counts'] = counts

    index_names = sampled_df[sampled_df['counts'] == 0].index
    sampled_df.drop(index_names, inplace=True)

    sampled_df.reset_index(inplace=True, drop=True)

    sampled_df.to_csv('sampled_df.csv') #TODO add option to write raw digests to file

    print(f'the length of sample_seqs is {len(sampled_df)}')
    print(sampled_df['counts'].sum())

    return sampled_df


def read_writer(sampled_df, r1, r2, chr_name):
    '''
    args.a1s is where to begin in the R1 adapter (33 for ours, I think)
    args.a2s is where to begin in the R2 adapter (34 for ours, I think)
    '''
    id = 0
    for idx, row in sampled_df.iterrows():
        r1_seq = row[0][args.a1s:args.a1s+args.l].ljust(args.l, 'G')
        r2_seq = reverse_comp(row[0])[args.a2s:args.a2s+args.l].ljust(args.l, 'G')
        for count in range(row[5]):
            r1_mut, r1_score = read_mutator(r1_seq, args.prob_mx1, args.q_dt1, args.q_ls1)
            r2_mut, r2_score = read_mutator(r2_seq, args.prob_mx2, args.q_dt2, args.q_ls2)
            r1.write(f'@{id}:{idx}:{row[3]}:{row[4]}:{row[1]}:{chr_name} 1\n{r1_mut}\n+\n{r1_score}\n')
            r2.write(f'@{id}:{idx}:{row[3]}:{row[4]}:{row[2]}:{chr_name} 2\n{r2_mut}\n+\n{r2_score}\n')
            id += 1


def read_mutator(seq, prob_mx, scores_dt, q_ls):
#    #TODO OPTIMIZATION NEEDED!
#    score = ''
#    mut_seq = ''
#    for idx, base in enumerate(seq):
#        q = np.random.choice(q_ls, 1, p=prob_mx[idx])[0]
#        score += q
#        if base == 'N':
#            mut_seq += 'N'
#            continue
#        p = scores_dt[q]
#        base_ls = ['A', 'C', 'G', 'T']
#        base_ls.remove(base)
#        base_ls.append(base)
#        p_ls = [p/3, p/3, p/3, 1-p]
#        mut_seq += np.random.choice(base_ls, 1, p=p_ls)[0]
#    return mut_seq, score
    #TODO OPTIMIZATION NEEDED!
    mut_seq = ''
    score = ''.join([np.random.choice(q_ls, 1, p=i)[0] for i in prob_mx])
    for base, q in zip(seq, score):
        p = scores_dt[q]
        base_ls = ['A', 'C', 'G', 'T', 'N']
        base_ls.remove(base)
        base_ls.append(base)
        p_ls = [p/4, p/4, p/4, p/4, 1-p]
        mut_seq += np.random.choice(base_ls, 1, p=p_ls)[0]
    return mut_seq, score



if __name__ == '__main__':
    args = parse_user_input()

    args.l = args.l if args.l else 250

    if args.a1:
        args.a1 = get_adapters(args.a1)

    if args.a2:
        args.a2 = get_adapters(args.a2)
        args.a2 = {reverse_comp(adapter): id for adapter, id in args.a2.items()}

    if args.q1:
        args.prob_mx1, args.q_dt1, args.q_ls1 = get_qscores(args.q1)
        args.l = len(args.prob_mx1)

    if args.q2:
        args.prob_mx2, args.q_dt2, args.q_ls2 = get_qscores(args.q2)
        args.l = len(args.prob_mx2)

    motif_dt = iupac_motifs()
    frag_len = args.f if args.f else 1000
    test_chrom(motif_dt, frag_len)

    #TODO check if genome compressed

