#!/usr/bin/env python

import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import random
import re
import seaborn as sns
import sys
import time

import prob_n_copies
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

    parser.add_argument('-up_bound', type=int, required=True,
            help='the upper end of a range (in bp) of read lengths to size select')

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


def check_genomes(genome_file):
    df = pd.read_csv(genome_file, names=['genome', 'abundance'])
    df['abundance'] = df['abundance'] / df['abundance'].sum()

    for idx in range(df.shape[0]):
        genome = df.iloc[idx]['genome']
        assert os.path.exists(genome), f'path to {genome} not found'

    return df


def process_genomes(args, genomes_df):
    digest_ls, prob_ls = [], []
    total_freqs = pd.DataFrame(columns=['length', 'sum_prob', 'name', 'counts_file'])

    for idx in range(genomes_df.shape[0]):
        args.genome = genomes_df.iloc[idx]['genome']
        args.comp = genomes_df.iloc[idx]['abundance']
        digest_file = os.path.join(args.o, 'raw_digest_' +
                os.path.basename(args.genome) + '.csv')

        print(f'processing genome {os.path.basename(args.genome)}')
        df = digest_genomes.main(args)
        digest_file  = process_df(df, digest_file, args)
        prob_file, len_freqs = prob_n_copies.main(digest_file, args)

        save_individual_hist(prob_file, args)

        digest_ls.append(digest_file)
        prob_ls.append(prob_file)
        tmp_df = pd.DataFrame(len_freqs.items(), columns=['length', 'sum_prob'])
        tmp_df['name'] = os.path.basename(args.genome)
        tmp_df['counts_file'] = prob_file
        total_freqs = pd.concat([total_freqs, tmp_df], axis=0)

    genomes_df['digest_file'] = digest_ls
    genomes_df['prob_file'] = prob_ls

    return genomes_df, total_freqs


def reverse_comp(seq):
    '''
    return the reverse complement of an input sequence
    '''
    revc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
            'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
            'N': 'N'}
    new = ''
    for base in reversed(seq):
        new += revc[base]

    return new


def process_df(df, digest_file, args):
    df['length'] = df['seq'].str.len()
    df['forward'] = np.where((df['m1'].isin(args.motif_dt1.keys()) & \
                             df['m2'].isin(args.motif_dt2.keys())), 1, 0)
    df['reverse'] = np.where((df['m1'].isin(args.motif_dt2.keys()) & \
                             df['m2'].isin(args.motif_dt1.keys())), 1, 0)

    # remove unviable combos after getting site_ls
    df.drop(df[(df['forward'] == 0) & (df['reverse'] == 0)].index,
            inplace = True)
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
    for mot, front in args.motif_dt.items():
        back = len(mot) - front
        df.loc[(df['m1'] == mot) & (df['reverse'] == 0), 'seq'] = \
                df['seq'].str[front:]
        df.loc[(df['m1'] == mot) & (df['reverse'] == 1), 'seq'] = \
                df['seq'].str[:-back]
        df.loc[(df['m2'] == mot) & (df['reverse'] == 0), 'seq'] = \
                df['seq'].str[:-back]
        df.loc[(df['m2'] == mot) & (df['reverse'] == 1), 'seq'] = \
                df['seq'].str[front:]

    df.to_csv(digest_file, index=None)

    return digest_file


def save_individual_hist(prob_file, args):
    df = pd.read_csv(prob_file)
    if df.shape[0] == 0:
        print(f'no fragments found in {read_file}')
        return

    sns.histplot(data=df, x=df['length'], binwidth=6,
                 alpha=0.75, color='red')
    sns.histplot(data=df, x=df['length'], weights=df['probability'], binwidth=6,
                 alpha=0.75, color='gold')
    sns.histplot(data=df, x=df['length'], weights=df['adj_prob'], binwidth=6,
                 alpha=0.75, color='blue')
    plt.savefig(os.path.join(args.o, 'hist_' + os.path.basename(prob_file)[:-4] +
                '.png'),
                 facecolor='white', transparent=False)
    plt.close()


def save_combined_hist(total_freqs, image_name, weights, args):
    ax = sns.histplot(data=total_freqs, x='length', hue='name',
                      weights=total_freqs[weights], multiple="stack",
                      binwidth=6, element="step")
    old_legend = ax.legend_
    handles = old_legend.legendHandles
    labels = [t.get_text() for t in old_legend.get_texts()]
    title = old_legend.get_title().get_text()
    ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left',
              borderaxespad=0)
    plt.savefig(os.path.join(args.o, f'_{image_name}.pdf'), bbox_inches='tight')


def write_genomes(comb_file, fragment_comps, adjustment):
    comb_df = pd.read_csv(comb_file)
    count_files_ls  = list(set(comb_df['counts_file'].to_list()))

    args.a1 = list(args.a1.items())
    args.a2 = list(args.a2.items())
    total = 0 #TODO #TODO #TODO #TODO
    with open(os.path.join(args.o, 'sim_metagenome_R1.fastq'), 'w') as r1,\
         open(os.path.join(args.o, 'sim_metagenome_R2.fastq'), 'w') as r2:

        for count_file in count_files_ls:
            gen_name = os.path.basename(count_file)[7:-4]
            df = pd.read_csv(count_file)
            df['counts'] = df['length'].map(fragment_comps)
            df['counts'] = round(df['counts'] * df['adj_prob'] \
                                  * adjustment)
            df.dropna(subset=['counts'], inplace=True)
            df = df[['seq', 'length', 'counts']]
            df['counts'] = df['counts'].astype(int)
            total += df['counts'].sum()
            write_reads(df, r1, r2, gen_name)
    print(total)


def write_reads(df, r1, r2, gen_name):
    """
    using the sampled read distribution, write to file as fastq format
    if adapters, randomly add sample and concatenate adapters
    if sample fastq data provided, mutate samples for error simulation
    """

    if not args.a1s:
        args.a1s, args.a2s = 0, 0

    if args.r1:
        read_writer_samples(df, r1, r2)
    elif args.q1:
        read_writer(df, r1, r2)
    else:
        read_writer_basic(df, r1, r2, gen_name)


#def simulate_adapters(dup_file):
#    """
#    simulate ligation of adapters to sequences
#
#    adapters must contain RE motif if sticky ends are used
#    the associated barcode is saved alongside the read as r1/2_id
#    """
#    df = pd.read_csv(dup_file)
#    lig_seq_ls = []
#    r1_ls = []
#    r2_ls = []
#    lig_len_ls = []
#
#    for seq in df['seq']:
#        r1_dapt, r1_id = random.choice(list(args.a1.items()))
#        r2_dapt, r2_id = random.choice(list(args.a2.items()))
#        ligated_seq = r1_dapt + seq + r2_dapt
#        lig_seq_ls.append(ligated_seq)
#        r1_ls.append(r1_id)
#        r2_ls.append(r2_id)
#        lig_len_ls.append(len(ligated_seq))
#
#    df['seq'] = lig_seq_ls
#    df['full_length'] = lig_len_ls
#    df['r1_id'] = r1_ls
#    df['r2_id'] = r2_ls
#
#    adapt_file = os.path.join(args.o, 'adapt_' + os.path.basename(args.genome) + '.csv')
#
#    df = df[['seq','start','end','strand','length','full_length','r1_id','r2_id']]
#    df.to_csv(adapt_file, index=None)
#
#    return adapt_file


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


def read_writer_basic(df, r1, r2, gen_name):
    """
    create a fastq-formatted output with no attempt at error profiling

    args.a1s is where to begin in the R1 adapter
    args.a2s is where to begin in the R2 adapter
    """
    df = np.array(df)
    id = 0
    score = 'I' * args.l
    a1s = args.a1s
    a1e = args.a1s + args.l
    a2s = args.a2s
    a2e = args.a2s + args.l

    for idx, i in enumerate(df):
        seq = i[0]
        for j in range(i[2]):
            a1 = random.choice(args.a1)
            a2 = random.choice(args.a2)
            seq = a1[0] + seq + a2[0]
            full = len(seq)
            r1_seq = seq[a1s:a1e].ljust(args.l, 'G')
            r2_seq = reverse_comp(seq)[a2s:a2e].ljust(args.l, 'G')
            header = f'@{id}:{idx}:{full}:{i[1]}:{a1[1]}:{a2[1]}:{gen_name}'
            r1.write(f'{header} 1\n{r1_seq}\n+\n{score}\n')
            r2.write(f'{header} 2\n{r2_seq}\n+\n{score}\n')
            id += 1


if __name__ == '__main__':
    args = parse_user_input()

    if args.o is None:
        args.o = os.path.dirname(os.path.abspath(__file__))
    elif os.path.exists(args.o) is True:
        args.o = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))

    args.m1, args.m2 = check_for_enzymes(args)
    args.motif_dt = {}
    args.motif_dt1 = iupac_motifs(args.m1)
    args.motif_dt.update(args.motif_dt1)
    args.motif_dt2 = iupac_motifs(args.m2)
    args.motif_dt.update(args.motif_dt2)

    args.t = args.t if args.t else 1
    args.sd = int(round(0.08*args.mean, 0)) # using Sage Science CV of 8%
    args.sd = max(args.sd, int(round((args.up_bound - args.mean)/2, 0)))
    args.min = args.min if args.min else 6
    args.max = args.max if args.max else (args.mean + (6*args.sd))

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

    '''
    begin processing pipeline
    '''
    genomes_df = check_genomes(args.genome)
    genomes_df, total_freqs = process_genomes(args, genomes_df)

    save_combined_hist(total_freqs, 'fragment_distributions', 'sum_prob', args)
    print('simulating size selection')
    fragment_comps, adjustment = size_selection.main(total_freqs, args)
    total_freqs['counts'] = total_freqs['length'].map(fragment_comps)
    total_freqs['counts'] = round(total_freqs['counts'] * total_freqs['sum_prob'] \
                          * adjustment)
    print(total_freqs['counts'].sum())
    comb_file = os.path.join(args.o, 'combined.csv')
    total_freqs.to_csv(comb_file)
    save_combined_hist(total_freqs, 'read_distributions', 'counts', args)

    if args.test:
        sys.exit()

    print('simulating fastq formatted sequence reads')
    write_genomes(comb_file, fragment_comps, adjustment)

    #TODO add option to input existing digest_file and skip digest_genome
    #TODO consider chromosome by chromosome approach
        # this would only temporarily store sequences in the raw digest
        # raw digests (per chromosome) will be digested and stored as copies
