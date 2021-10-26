#!/usr/bin/env python

import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import random
import re
import sys
import time

import n_copies
import prob_n_copies #TODO
import digest_genomes
import size_selection


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

    parser.add_argument('-l', type=int, required=False,
            help='desired read length of final simulated reads (defaults to 250 or given q1/q2 profiles)')

    parser.add_argument('-test', dest='test', action='store_true',
            help='test mode: create newline-separated file of RE digested sequences only')

    parser.add_argument('-t', type=int, required=False,
            help='number of subprocesses to run while simulating copy number')

    parser.add_argument('-n', type=int, required=True,
            help='genome copy (depth per locus)')

    parser.add_argument('-mean', type=int, required=True,
            help='mean (in bp) of read lengths after size selection')

    parser.add_argument('-sd', type=int, required=True,
            help='standard deviation (in bp) of read lengths after size selection')

    parser.add_argument('-min', type=int, required=False,
            help='min distance between cuts (optional, defaults to 6bp)')

    parser.add_argument('-max', type=int, required=False,
            help='max fragment length after first cut (optional, defaults to mean + 6 stdevs)')

    parser.add_argument('-cut_prob', type=float, required=True,
            help='percent probability of per-site cut; use \'1\' for complete digestion of fragments (fragments will not contain internal RE sites)')

    parser.add_argument('-a1', type=str, required=False,
            help='file containing tab/space-separated adapters and barcode that attach 5\' to read')

    parser.add_argument('-a2', type=str, required=False,
            help='file containing tab/space-separated adapters and barcode that attach 3\' to read')

    parser.add_argument('-a1s', type=int, required=False,
            help='manually provide bp length of adapter a1 before SBS begins')

    parser.add_argument('-a2s', type=int, required=False,
            help='manually provide bp length of adapter a1 before SBS begins')

    parser.add_argument('-q1', type=str, required=False,
            help='file containing R1 q scores in csv format (see ngsComposer tool crinoid)')

    parser.add_argument('-q2', type=str, required=False,
            help='file containing R2 q scores in csv format (see ngsComposer tool crinoid)')

    parser.add_argument('-r1', type=str, required=False,
            help='R1 fastq file to sample Q scores')

    parser.add_argument('-r2', type=str, required=False,
            help='R2 fastq file to sample Q scores')

    parser.add_argument('-p', type=int, required=False,
            help='if using r1/r2 for profile, percent of reads to sample')

    args = parser.parse_args()

    return args


def iupac_motifs(arg_m):
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


def get_adapters(arg):
    """
    open adapters file and store adapters and barcodes in dt
    """
    adapters_dt = {}
    with open(arg) as f:
        for line in f:
            adapter, id = line.rstrip().split()
            adapters_dt[adapter] = id

    return adapters_dt


def get_sbs_start(adapter_ls):
    """
    if a1s/a2s not provided, an automated attempt to find the SBS start
    site is performed by finding the first position where adapters
    deviate in sequence
    """
    if len(adapter_ls) == 1:
        return len(adapter_ls[0])

    for idx, pos in enumerate(adapter_ls[0]):
        for seq in adapter_ls:
            if adapter_ls[0][:idx] not in seq:
                return idx-1


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

    return prob_mx, scores_dt, scores_ls


def open_fastq(fastq):
    print(f'sampling {args.p} percent of scores from {fastq}')
    perc_keep = [int(args.p)/100, 1-(int(args.p)/100)]
    i = 0
    with open(fastq) as f, open(fastq + '_sampled_scores.csv', 'w') as out:
        for line in f:
            i += 1
            if i == 4:
                keep = np.random.choice([True, False], 1, p=perc_keep)
                if keep:
                    out.write(line)
                i = 0


def reverse_comp(seq):
    '''
    return the reverse complement of an input sequence
    '''
    revc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    new = ''
    for base in reversed(seq):
        new += revc[base]

    return new


def process_df(df, digest_file):
    print(df.head())
    df['length'] = df['seq'].str.len()
    df['forward'] = np.where((df['m1'].isin(motif_dt1.keys()) & \
                             df['m2'].isin(motif_dt2.keys())), 1, 0)
    df['reverse'] = np.where((df['m1'].isin(motif_dt2.keys()) & \
                             df['m2'].isin(motif_dt1.keys())), 1, 0)

    # remove unviable combos after getting site_ls
    df.drop(df[(df['forward'] == 0) & (df['reverse'] == 0)].index, inplace = True)
    df = df.reset_index(drop=True)

    # create a column of reverse complement sequences
    df['revc'] = [reverse_comp(i) for i in df['seq'].to_list()]

    # duplicate all the reads that work both ways (if RE in m1 and m2)
    tmp_df = df[(df['forward'] == 1) & (df['reverse'] == 1)]
    tmp_df = tmp_df.reset_index(drop=True)

    # recategorize the seqs in df as being forward
    df.loc[(df['forward'] == 1) & (df['reverse'] == 1), 'reverse'] = 0

    # convert remaining reverse strand sequences to the reverse complement
    df.loc[df['reverse'] == 1, 'seq'] = df['revc']

    # make all tmp_df reads reverse complements and recategorize as reverse
    tmp_df['seq'] = tmp_df['revc'].values
    tmp_df.loc[:,'forward'] = 0
    tmp_df.loc[:,'reverse'] = 1

    df = pd.concat([df, tmp_df])
    df = df.sort_values(by=['length'])
    df = df.reset_index(drop=True)
    df.drop('revc', axis=1, inplace=True)
    df.drop('forward', axis=1, inplace=True)

    # add a quick step that removes appropriate over/underhang
    for mot, front in motif_dt.items():
        back = len(mot) - front
        df.loc[(df['m1'] == mot) & (df['reverse'] == 0), 'seq'] = \
                df['seq'].str[front:]
        df.loc[(df['m1'] == mot) & (df['reverse'] == 1), 'seq'] = \
                df['seq'].str[:-back]
        df.loc[(df['m2'] == mot) & (df['reverse'] == 0), 'seq'] = \
                df['seq'].str[:-back]
        df.loc[(df['m2'] == mot) & (df['reverse'] == 1), 'seq'] = \
                df['seq'].str[front:]

    #TODO readjust length?
    df.to_csv(digest_file, index=None)

    return digest_file


def save_hist(proj, read_file, title, leglab):
    df = pd.read_csv(read_file)

    length_ls = []
    if 'counts'  in [col for col in df]:
        for idx, row in df.iterrows():
            length_ls.extend([row['full_length'] for i in range(row['counts'])])
    elif 'copies' in [col for col in df]:
        for idx, row in df.iterrows():
            length_ls.extend([row['length'] for i in range(row['copies'])])
    else:
        length_ls = df['length'].to_list()

    plt.hist(length_ls, bins=100, range=[0, args.max], label=leglab, alpha=0.75)
    plt.xlabel('Fragment Length')
    plt.ylabel('Count')
    plt.title('Distribution of ' + title)
    plt.legend()
    plt.savefig(os.path.join(proj, 'hist_' + os.path.basename(read_file)[:-4] \
                + '.png'), facecolor='white', transparent=False)


def write_reads(proj, sampled_file):
    """
    using the sampled read distribution, write to file as fastq format
    if adapters, randomly add sample and concatenate adapters
    if sample fastq data provided, mutate samples for error simulation
    """
    df = pd.read_csv(sampled_file)

    #TODO add adapters on the fly after size selection
    if not args.a1s:
        args.a1s, args.a2s = 0, 0

    with open(os.path.join(proj, os.path.basename(args.genome) + '_R1.fastq'), 'w') as r1,\
         open(os.path.join(proj, os.path.basename(args.genome) + '_R2.fastq'), 'w') as r2:

        if args.r1:
            read_writer_samples(df, r1, r2)
        elif args.q1:
            read_writer(df, r1, r2)
        else:
            read_writer_basic(df, r1, r2)


def simulate_adapters(dup_file):
    """
    simulate ligation of adapters to sequences

    adapters must contain RE motif if sticky ends are used
    the associated barcode is saved alongside the read as r1/2_id
    """
    df = pd.read_csv(dup_file)
    lig_seq_ls = []
    r1_ls = []
    r2_ls = []
    lig_len_ls = []

    for seq in df['seq']:
        r1_dapt, r1_id = random.choice(list(args.a1.items()))
        r2_dapt, r2_id = random.choice(list(args.a2.items()))
        ligated_seq = r1_dapt + seq + r2_dapt
        lig_seq_ls.append(ligated_seq)
        r1_ls.append(r1_id)
        r2_ls.append(r2_id)
        lig_len_ls.append(len(ligated_seq))

    df['seq'] = lig_seq_ls
    df['full_length'] = lig_len_ls
    df['r1_id'] = r1_ls
    df['r2_id'] = r2_ls

    adapt_file = os.path.join(proj, 'adapt_' + os.path.basename(args.genome) + '.csv')

    df = df[['seq','start','end','strand','length','full_length','r1_id','r2_id']]
    df.to_csv(adapt_file, index=None)

    return adapt_file


def read_writer(df, r1, r2):
    """
    args.a1s is where to begin in the R1 adapter
    args.a2s is where to begin in the R2 adapter
    """
    gen_name = os.path.basename(args.genome)
    id = 0
    for idx, row in df.iterrows():
        r1_seq = row['seq'][args.a1s:args.a1s+args.l].ljust(args.l, 'G')
        r2_seq = reverse_comp(row['seq'])[args.a2s:args.a2s+args.l].ljust(args.l, 'G')
        for count in range(row['counts']):
            r1_mut, r1_score = read_mutator(r1_seq, args.prob_mx1, args.q_dt1, args.q_ls1)
            r2_mut, r2_score = read_mutator(r2_seq, args.prob_mx2, args.q_dt2, args.q_ls2)
            r1.write(f'@{id}:{idx}:{row[5]}:{row[4]}:{row[6]}:{row[7]}:' \
                     f'{gen_name} 1\n{r1_mut}\n+\n{r1_score}\n')
            r2.write(f'@{id}:{idx}:{row[5]}:{row[4]}:{row[6]}:{row[7]}:' \
                     f'{gen_name} 2\n{r2_mut}\n+\n{r2_score}\n')
            id += 1


def read_mutator(seq, prob_mx, scores_dt, q_ls):
    """
    using pure, per-base probability profile, mutate bases
    """
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


def read_writer_samples(df, r1, r2):
    """
    args.a1s is where to begin in the R1 adapter
    args.a2s is where to begin in the R2 adapter
    """
    sampled_q1 = pd.read_csv(args.r1 + '_sampled_scores.csv', names=['score'], sep = '\t')
    sampled_q2 = pd.read_csv(args.r2 + '_sampled_scores.csv', names=['score'], sep = '\t')
    sampled_q1 = list(sampled_q1['score'].values)
    sampled_q2 = list(sampled_q2['score'].values)
    args.l = len(max(sampled_q1, key=len))
    gen_name = os.path.basename(args.genome)
    id = 0
    for idx, row in df.iterrows():
        r1_seq = row['seq'][args.a1s:args.a1s+args.l].ljust(args.l, 'G')
        r2_seq = reverse_comp(row['seq'])[args.a2s:args.a2s+args.l].ljust(args.l, 'G')
        for count in range(row['counts']):
            r1_mut, r1_score = read_mutator_samples(r1_seq, args.q_dt1, sampled_q1)
            r2_mut, r2_score = read_mutator_samples(r2_seq, args.q_dt2, sampled_q2)
            r1.write(f'@{id}:{idx}:{row[5]}:{row[4]}:{row[6]}:{row[7]}:' \
                     f'{gen_name} 1\n{r1_mut}\n+\n{r1_score}\n')
            r2.write(f'@{id}:{idx}:{row[5]}:{row[4]}:{row[6]}:{row[7]}:' \
                     f'{gen_name} 2\n{r2_mut}\n+\n{r2_score}\n')
            id += 1


def read_mutator_samples(seq, scores_dt, sampled_q):
    """
    using actual sampled per-read q scores, mutate bases
    """
    mut_seq = ''
    score = random.choice(sampled_q)
    for base, q in zip(seq, score):
        p = scores_dt[q]
        base_ls = ['A', 'C', 'G', 'T', 'N']
        base_ls.remove(base)
        base_ls.append(base)
        p_ls = [p/4, p/4, p/4, p/4, 1-p]
        mut_seq += np.random.choice(base_ls, 1, p=p_ls)[0]
    return mut_seq, score


def read_writer_basic(df, r1, r2):
    """
    create a fastq-formatted output with no attempt at error profiling

    args.a1s is where to begin in the R1 adapter
    args.a2s is where to begin in the R2 adapter
    """
    df = np.array(df)
    gen_name = os.path.basename(args.genome)
    id = 0
    score = 'I' * args.l
    args.a1 = list(args.a1.items())
    args.a2 = list(args.a2.items())
    a1s = args.a1s
    a1e = args.a1s + args.l
    a2s = args.a2s
    a2e = args.a2s + args.l
    # seq   start   end m1  m2  length  reverse copies  full    counts
    # 0     1       2   3   4   5       6       7       8       9
    for idx, i in enumerate(df):
        seq = i[0]
        for j in range(i[9]):
            a1 = random.choice(args.a1)
            a2 = random.choice(args.a2)
            seq = a1[0] + seq + a2[0]
            r1_seq = seq[a1s:a1e].ljust(args.l, 'G')
            r2_seq = reverse_comp(seq)[a2s:a2e].ljust(args.l, 'G')
            r1.write(f'@{id}:{idx}:{i[8]}:{i[5]}:{a1[1]}:{a2[1]}:' \
                     f'{gen_name} 1\n{r1_seq}\n+\n{score}\n')
            r2.write(f'@{id}:{idx}:{i[8]}:{i[5]}:{a1[1]}:{a2[1]}:' \
                     f'{gen_name} 2\n{r2_seq}\n+\n{score}\n')
            id += 1


if __name__ == '__main__':
    args = parse_user_input()

    if args.o is None:
        proj = os.path.dirname(os.path.abspath(__file__))
    elif os.path.exists(args.o) is True:
        proj = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))

    motif_dt = {}
    motif_dt1 = iupac_motifs(args.m1)
    motif_dt.update(motif_dt1)
    motif_dt2 = iupac_motifs(args.m2)
    motif_dt.update(motif_dt2)

    args.t = args.t if args.t else 1
    args.min = args.min if args.min else 6
    args.max = args.max if args.max else (args.mean + (6*args.sd))
    frag_len = args.max

    if args.r1 and not args.q1:
        sys.exit('please provide q scores profile for R1')

    args.l = args.l if args.l else 250

    if args.a1:
        args.a1 = get_adapters(args.a1)
        if not args.a1s:
            args.a1s = get_sbs_start(list(args.a1.keys()))

    if args.a2:
        args.a2 = get_adapters(args.a2)
        args.a2 = {reverse_comp(adapter): id for adapter, id in args.a2.items()}
        if not args.a2s:
            args.a2s = get_sbs_start(list(args.a2.keys()))

    if args.q1:
        args.prob_mx1, args.q_dt1, args.q_ls1 = get_qscores(args.q1)
        args.l = len(args.prob_mx1)

    if args.q2:
        args.prob_mx2, args.q_dt2, args.q_ls2 = get_qscores(args.q2)
        args.l = len(args.prob_mx2)

    if args.r1:
        open_fastq(args.r1)

    if args.r2:
        open_fastq(args.r2)

    digest_file = os.path.join(proj, 'raw_digest_' + os.path.basename(args.genome) + '.csv')

    #TODO add option to input existing digest_file and skip digest_genome
    #TODO consider chromosome by chromosome approach
        # this would only temporarily store sequences in the raw digest
        # raw digests (per chromosome) will be digested and stored as copies

    print('\nsimulating restriction digest\n')
    df = digest_genomes.main(motif_dt, frag_len, args)
    digest_file  = process_df(df, digest_file)
    save_hist(proj, digest_file, 'Possible Raw Fragments', 'possible')

    print('\nsimulating genome copy number\n')
    #dup_file = n_copies.main(proj, digest_file, args) #TODO
    dup_file = prob_n_copies.main(proj, digest_file, args) #TODO
    save_hist(proj, dup_file, f'Fragments of {args.n}X Copy Number', \
              f'{args.n} copies')

    print('\nsimulating size selection\n')
    sampled_file = size_selection.main(dup_file, proj, args)
    save_hist(proj, sampled_file, 'Size Selected Fragments', 'size selected')

    if args.test:
        sys.exit()

    print('\nsimulating fastq formatted sequence reads\n')
    write_reads(proj, sampled_file)

