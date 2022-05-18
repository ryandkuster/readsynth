#!/usr/bin/env python

import numpy as np
import pandas as pd
import random


def main(df, r1, r2, gen_name, args):
    '''
    using the sampled read distribution, write to file as fastq format
    if adapters, randomly add sample and concatenate adapters
    if sample fastq data provided, mutate samples for error simulation
    '''

    if not args.a1s:
        args.a1s, args.a2s = 0, 0

    if args.r1 and args.r2:
        read_writer_samples(df, r1, r2, gen_name, args)
    else:
        read_writer_basic(df, r1, r2, gen_name, args)


def read_writer_samples(df, r1, r2, gen_name, args):
    '''
    args.a1s is where to begin in the R1 adapter
    args.a2s is where to begin in the R2 adapter
    '''

    sampled_q1 = pd.read_csv(args.r1 + '_sampled_scores.csv',
                             names=['score'], sep='\t')
    sampled_q2 = pd.read_csv(args.r2 + '_sampled_scores.csv',
                             names=['score'], sep='\t')
    sampled_q1 = list(sampled_q1['score'].values)
    sampled_q2 = list(sampled_q2['score'].values)

    df = np.array(df)
    r_no = 0
    a1s = args.a1s
    a1e = args.a1s + args.l
    a2s = args.a2s
    a2e = args.a2s + args.l

    for idx, i in enumerate(df):
        for j in range(i[3]):
            a1, a2, r1_seq, r2_seq, full = \
                create_seq(i, a1s, a1e, a2s, a2e, args)
            r1_seq, r1_score = even_score_seq(r1_seq, sampled_q1)
            r2_seq, r2_score = even_score_seq(r2_seq, sampled_q2)
            header = f'@{r_no}:{idx}:{full}:{i[2]}:{a1[2]}:{a2[2]}:{gen_name}'
            r1.write(f'{header} 1\n{r1_seq}\n+\n{r1_score}\n')
            r2.write(f'{header} 2\n{r2_seq}\n+\n{r2_score}\n')
            r_no += 1


def create_seq(i, a1s, a1e, a2s, a2e, args):
    a1 = random.choice(args.a1)
    a2 = random.choice(args.a2)
    r1_seq = a1[0] + i[0] + a2[1]
    r2_seq = a2[0] + i[1] + a1[1]
    full = len(r1_seq)
    r1_seq = r1_seq[a1s:a1e].ljust(args.l, 'G')
    r2_seq = r2_seq[a2s:a2e].ljust(args.l, 'G')

    return a1, a2, r1_seq, r2_seq, full


def even_score_seq(seq, sampled_q):
    score = random.choice(sampled_q)
    adj_l = len(min(seq, score, key=len))

    return seq[:adj_l], score[:adj_l]


def read_writer_basic(df, r1, r2, gen_name, args):
    '''
    create a fastq-formatted output with no attempt at error profiling

    args.a1s is where to begin in the R1 adapter
    args.a2s is where to begin in the R2 adapter
    '''
    df = np.array(df)
    r_no = 0
    score = 'I' * args.l
    a1s = args.a1s
    a1e = args.a1s + args.l
    a2s = args.a2s
    a2e = args.a2s + args.l

    for idx, i in enumerate(df):
        for j in range(i[3]):
            a1, a2, r1_seq, r2_seq, full = \
                create_seq(i, a1s, a1e, a2s, a2e, args)
            header = f'@{r_no}:{idx}:{full}:{i[2]}:{a1[2]}:{a2[2]}:{gen_name}'
            r1.write(f'{header} 1\n{r1_seq}\n+\n{score}\n')
            r2.write(f'{header} 2\n{r2_seq}\n+\n{score}\n')
            r_no += 1
