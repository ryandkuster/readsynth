import numpy as np
import os
import pandas as pd
import random

def main(df, r1, r2, gen_name, args):
    """
    using the sampled read distribution, write to file as fastq format
    if adapters, randomly add sample and concatenate adapters
    if sample fastq data provided, mutate samples for error simulation
    """

    if not args.a1s:
        args.a1s, args.a2s = 0, 0

    if args.r1 and args.r2:
        read_writer_samples(df, r1, r2, gen_name, args)
    else:
        read_writer_basic(df, r1, r2, gen_name, args)


def read_writer_samples(df, r1, r2, gen_name, args):
    """
    args.a1s is where to begin in the R1 adapter
    args.a2s is where to begin in the R2 adapter
    """

    scores = list('!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJK')
    scores_dt = {j: 10**(-i/10) for i, j in enumerate(scores)}

    sampled_q1 = pd.read_csv(args.r1 + '_sampled_scores.csv', names=['score'], sep = '\t')
    sampled_q2 = pd.read_csv(args.r2 + '_sampled_scores.csv', names=['score'], sep = '\t')
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
            a1 = random.choice(args.a1)
            a2 = random.choice(args.a2)
            r1_seq = a1[0] + i[0] + a2[1]
            r2_seq = a2[0] + i[1] + a1[1]
            r1_seq = r1_seq[a1s:a1e].ljust(args.l, 'G')
            r2_seq = r2_seq[a2s:a2e].ljust(args.l, 'G')
            r1_score = random.choice(sampled_q1)
            r2_score = random.choice(sampled_q2)
            full = len(r1_seq)
            header = f'@{r_no}:{idx}:{full}:{i[2]}:{a1[2]}:{a2[2]}:{gen_name}'
            r1.write(f'{header} 1\n{r1_seq}\n+\n{r1_score}\n')
            r2.write(f'{header} 2\n{r2_seq}\n+\n{r2_score}\n')
            r_no += 1


def read_writer_basic(df, r1, r2, gen_name, args):
    """
    create a fastq-formatted output with no attempt at error profiling

    args.a1s is where to begin in the R1 adapter
    args.a2s is where to begin in the R2 adapter
    """
    df = np.array(df)
    r_no = 0
    score = 'I' * args.l
    a1s = args.a1s
    a1e = args.a1s + args.l
    a2s = args.a2s
    a2e = args.a2s + args.l

    for idx, i in enumerate(df):
        for j in range(i[3]):
            a1 = random.choice(args.a1)
            a2 = random.choice(args.a2)
            r1_seq = a1[0] + i[0] + a2[1]
            r2_seq = a2[0] + i[1] + a1[1]
            r1_seq = r1_seq[a1s:a1e].ljust(args.l, 'G')
            r2_seq = r2_seq[a2s:a2e].ljust(args.l, 'G')
            full = len(r1_seq)
            header = f'@{r_no}:{idx}:{full}:{i[2]}:{a1[2]}:{a2[2]}:{gen_name}'
            r1.write(f'{header} 1\n{r1_seq}\n+\n{score}\n')
            r2.write(f'{header} 2\n{r2_seq}\n+\n{score}\n')
