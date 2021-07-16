#!/usr/bin/env python

import argparse
import gzip
import math
import numpy as np
import os
import pandas as pd
import random
import re
import sys

from scripts.gzip_test import test_unicode


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

    parser.add_argument('-t', type=str, required=False,
            help='test mode: create newline-separated file of RE digested sequences only')

    parser.add_argument('-n', type=int, required=True,
            help='genome copy (depth per locus)')

    parser.add_argument('-mean', type=int, required=True,
            help='mean (in bp) of read lengths after size selection')

    parser.add_argument('-sd', type=int, required=True,
            help='standard deviation (in bp) of read lengths after size selection')

    parser.add_argument('-f', type=int, required=False,
            help='max fragment length after first cut (optional, defaults to mean + 6 stdevs)')

    parser.add_argument('-complete', type=int, required=False,
            help='use \'1\' for complete digestion of fragments (fragments will not contain internal RE sites)')

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


def digest_genome(motif_dt, frag_len, proj):
    """
    process fasta sequences as restricion enzyme fragments

    input args.genome is fasta

    output 'raw_digest' file holds all the resulting fragments
    """
    digest_file = os.path.join(proj, 'raw_digest_' + os.path.basename(args.genome) + '.csv')
    dfile = open(digest_file, 'w')
    dfile.write('seq,start,end,strand,length\n')

    if test_unicode(args.genome):
        fasta = gzip.open(args.genome, 'rt')
    else:
        fasta = open(args.genome)

    begin = 0
    gen_ls = []
    seq = ''
    for line in fasta:
        if line.startswith('>') and seq:
            seq_ls = digest_seq(begin, seq, motif_dt, frag_len)
            seq_ls = second_digest(seq_ls, motif_dt, frag_len)
            gen_ls.extend([i+[len(i[0])] for i in seq_ls])
            begin += len(seq)
            seq = ''
            chr_name = line.rstrip()[1:].replace(' ', '_').replace(',','')
        elif line.startswith('>'):
            chr_name = line.rstrip()[1:].replace(' ', '_').replace(',','')
        else:
            seq += line.rstrip().upper()

    seq_ls = digest_seq(begin, seq, motif_dt, frag_len)
    seq_ls = second_digest(seq_ls, motif_dt, frag_len)
    gen_ls.extend([i+[len(i[0])] for i in seq_ls])
    begin += len(seq)

    if len(gen_ls) > 0:
        digest_stats(gen_ls, begin)
        write_digested_reads(gen_ls, dfile)

    dfile.close()
    fasta.close()

    return digest_file


def digest_seq(begin, seq, motif_dt, frag_len):
    """
    for every chromosome (seq), find all RE recognition positions
    and preserve frag_len bp ahead as a possible template (fragment)

    each item in seq_ls is [sequence, start, end]
    """
    seq_ls = []
    for motif in motif_dt.keys():
        for idx in re.finditer('(?=' + motif + ')', seq):
            start = idx.start()
            end = start + frag_len
            fragment = seq[start:end]
            if fragment not in motif_dt:
                seq_ls.append([fragment, begin+start, begin+(len(fragment)+start)])

    return [[reverse_comp(i[0]), i[1], i[2]] for i in seq_ls]


def reverse_comp(seq):
    revc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    new = ''
    for base in reversed(seq):
        new += revc[base]

    return new


def second_digest(seq_ls, motif_dt, frag_len):
    second_ls = []
    for i in seq_ls:
        subseqs = digest_seq(0, i[0], motif_dt, frag_len)
        subseqs = [[j[0], i[1], i[2]-j[1], '+'] for j in subseqs]
        second_ls += subseqs

    second_ls += [[reverse_comp(i[0]), i[1], i[2], '-'] for i in second_ls]
    second_ls = [orientation_test(*i) for i in second_ls]
    second_ls = [i for i in second_ls if i[0]]

    if args.complete == 1:
        second_ls = [[complete_digest(i[0], motif_dt)] + i[1:] for i in second_ls]
        second_ls = [i for i in second_ls if i[0]]

    if len(second_ls) == 0:
        return second_ls

    return second_ls


def orientation_test(seq, start, end, strand):
    """
    reads passed here contain intact RE motifs
    (i.e., AG/CT remains AGCT, not CT)

    this is necessary if sticky ends will be used in ligating
    fragments to adapters
    """
    for motif1, offset1 in motif_dt1.items():
        if re.search('^'+motif1, seq):
            for motif2, offset2 in motif_dt2.items():
                if re.search(motif2+'$', seq):
                    start += offset1
                    end -= len(motif2)-offset2
                    if offset2 == 0:
                        return [seq[offset1:], start, end, strand]
                    elif len(motif2) == offset2:
                        return [seq[offset1:], start, end, strand]
                    else:
                        return [seq[offset1:-(len(motif2)-offset2)], start, end, strand]

    return [None, start, end, strand]


def complete_digest(seq, motif_dt):
    """
    remove any fragments containing an internal cut site for any of the
    restriction enzymes used
    """
    for motif in motif_dt:
        if motif in seq [1:-1]:
            return False

    return seq


def digest_stats(gen_ls, begin):
    """
    prints a histogram to screen for every input sequence
    """
    seq_ls = [seq[0] for seq in gen_ls]
    stats_dt = {}
    len_ls = [len(seq) for seq in seq_ls]
    for num in range(args.sd, max(max(len_ls), args.sd)+args.sd, args.sd):
        stats_dt[int(math.ceil(num / 10.0)) * 10] = 0

    print(f'input genome length: {begin} bp')
    print(f'average fragment length: {round(sum(map(len, seq_ls)) / len(seq_ls), 2)} bp')
    print(f'number of fragments: {len(seq_ls)}')

    col, row = os.get_terminal_size()
    for leng in len_ls:
        for bin in stats_dt.keys():
            if leng <= bin:
                stats_dt[bin] += 1
                break

    highest = stats_dt[max(stats_dt, key=lambda i: stats_dt[i])]
    denom = highest/(col-10)
    for k, v in stats_dt.items():
        print(str(k) + ':' + ' ' * (len(str(max(stats_dt)))-len(str(k))) + '|' * round(v/denom))


def write_digested_reads(gen_ls, dfile):
    for i in gen_ls:
        dfile.write(f'{i[0]},{str(i[1])},{str(i[2])},{i[3]},{str(i[4])}\n')


def add_position_weights(digest_file):
    """
    using the start and end position for all fragments, create a
    per-fragment weight that accounts for multiple fragments mapping to
    the same locus (can occur due to incomplete digest or nesting RE
    motifs or fragments that occur in the template and reverse template
    """
    df = pd.read_csv(digest_file)
    pos_dt = count_pos(df)
    weight_ls = add_pos_weight(pos_dt, df)
    print(f'positions covered: {len(pos_dt)}')
    df['weight'] = weight_ls
    results = sum([i*j for (i,j) in zip(df['weight'],df['length'])])
    print(f'sum of bp weights: {len(pos_dt)}')
    df.to_csv(digest_file, index=None)


def count_pos(df):
    """
    create a dictionary where keys are genomic loci and values are
    counts of fragments mapping to those loci
    """
    pos_dt = {}
    for start, end in zip(df['start'], df['end']):
        for i in range(start, end):
            if i in pos_dt:
                pos_dt[i] += 1
            else:
                pos_dt[i] = 1
    return pos_dt


def add_pos_weight(pos_dt, df):
    """
    iterate through the sequence start and end sites and sum the bases
    where multiple reads are present; divide fragment length by this
    sum to create a position-specific weight
    """
    weight_ls = []
    for idx, row in df.iterrows():
        frag_total = 0
        for i in range(row['start'], row['end']):
            frag_total += pos_dt[i]
        weight_ls.append(row['length']/frag_total)

    return weight_ls


def save_histogram(proj, digest_file):
    df = pd.read_csv(digest_file)
    ax = df['length'].hist(bins=100, range=[0, args.f])
    fig = ax.get_figure()
    fig.savefig(os.path.join(proj, 'hist_' + os.path.basename(digest_file)[:-4] + '.png'))


def size_selection(proj, digest_file):
    """
    using the previously created digest file containing raw fragments
    of genomic dna of varying length, optionally add adapters, then
    simulate a random sampling of reads in a desired size range
    """
    if args.a1:
        digest_file = simulate_adapters(digest_file)
    else:
        digest_file = no_adapters(digest_file)

    sampled_df = simulate_length(digest_file, proj)

    if args.t:
        sys.exit()

    with open(os.path.join(proj, os.path.basename(args.genome) + '_R1.fastq'), 'w') as r1,\
         open(os.path.join(proj, os.path.basename(args.genome) + '_R2.fastq'), 'w') as r2:

        if args.r1:
            read_writer_samples(sampled_df, r1, r2)
        elif args.q1:
            read_writer(sampled_df, r1, r2)
        else:
            read_writer_basic(sampled_df, r1, r2)


def simulate_adapters(digest_file):
    """
    simulate ligation of adapters to sequences

    adapters must contain RE motif if sticky ends are used
    the associated barcode is saved alongside the read as r1/2_id
    """
    df = pd.read_csv(digest_file)
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

    df = df[['seq','start','end','strand','length','full_length','r1_id','r2_id','weight']]
    df.to_csv(adapt_file, index=None)

    return adapt_file


def no_adapters(digest_file):
    """
    adjust df for lack of adapters
    """
    df = pd.read_csv(digest_file)
    df['full_length'] = df['length']
    df['r1_id'] = ['na' for i in df['seq']]
    df['r2_id'] = ['na' for i in df['seq']]

    df = df[['seq','start','end','strand','length','full_length','r1_id','r2_id','weight']]
    df.to_csv(digest_file, index=None)

    return digest_file


def simulate_length(digest_file, proj):
    """
    sampling is performed by first assessing the count of fragments for
    each possible read length and then generating a normal distribution
    that includes enough counts in the mean+2sd range to include a 1X
    coverage of the genome

    this random sampling distribution is then multiplied by args.n to
    create nX coverage that roughly reflects size selection in molecular
    library prep
    """
    df = pd.read_csv(digest_file)
    col_names = [col for col in df.columns]
    df.sort_values(['full_length'], ascending=[True], inplace=True)
    df.reset_index(inplace=True, drop=True)

    total_reads = 10
    keep_going = True
    len_dt = {}

    # create a len_dt, storing the weight-adjusted count for each length
    for i in range(args.mean, args.mean + 2*args.sd):
        len_count = df[df.full_length == i]['weight'].sum()
        len_dt[i] = len_count

    # produce a normal distribution that includes mean + 2sd counts
    while keep_going is True:
        keep_going = False
        draw_ls = np.random.normal(loc=args.mean,scale=args.sd,size=total_reads)
        draw_ls = [round(i) for i in draw_ls]
        for i, len_count in len_dt.items():
            if len_count > draw_ls.count(i):
                total_reads = round(total_reads*1.1)
                keep_going = True
                break

    # create a dictionary of draw numbers
    draw_dt = {}

    for i in range(min(draw_ls), max(draw_ls)+1):
        draw_counts = draw_ls.count(i) * args.n
        data_counts = round(df[df.full_length == i]['weight'].sum() * args.n)
        draw_dt[i] = min(draw_counts, data_counts)

    # for each fragment length, randomly draw reads
    sampled_df = pd.DataFrame(columns=col_names)
    counts = []

    for length, draws in draw_dt.items():
        tmp_df = df.loc[df['full_length'] == length]
        if len(tmp_df) == 0:
            continue
        indices = [i for i in range(len(tmp_df))]
        sampled_idx = random.choices(indices, k=draws)
        counts += [sampled_idx.count(idx) for idx in indices]
        sampled_df = pd.concat([sampled_df, tmp_df])

    sampled_df['counts'] = counts
    index_names = sampled_df[sampled_df['counts'] == 0].index
    sampled_df.drop(index_names, inplace=True)
    samp_no = sampled_df['counts'].sum()
    print(f'fragments sampled around mean of {args.mean}bp : {samp_no}')
    sampled_df.reset_index(inplace=True, drop=True)
    sampled_file = os.path.join(proj, 'sampled_' + os.path.basename(args.genome) + '.csv')
    sampled_df.to_csv(sampled_file, index=None)

    #TODO start testing visual
    histogram_seqs = pd.DataFrame(columns=['full_length'])
    length_ls = []

    for idx, row in sampled_df.iterrows():
        length_ls.extend([row['full_length'] for i in range(row['counts'])])

    histogram_seqs['full_length'] = length_ls
    ax = histogram_seqs['full_length'].hist(bins=100, range=[0, args.f])
    fig = ax.get_figure()
    fig.savefig(os.path.join(proj, 'hist_' + os.path.basename(sampled_file)[:-4] + '.png'))

    #TODO end testing visual

    return sampled_df


def read_writer(sampled_df, r1, r2):
    """
    args.a1s is where to begin in the R1 adapter
    args.a2s is where to begin in the R2 adapter
    """
    gen_name = os.path.basename(args.genome)
    id = 0
    for idx, row in sampled_df.iterrows():
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


def read_writer_samples(sampled_df, r1, r2):
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
    for idx, row in sampled_df.iterrows():
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


def read_writer_basic(sampled_df, r1, r2):
    """
    create a fastq-formatted output with no attempt at error profiling

    args.a1s is where to begin in the R1 adapter
    args.a2s is where to begin in the R2 adapter
    """
    gen_name = os.path.basename(args.genome)
    id = 0
    for idx, row in sampled_df.iterrows():
        r1_seq = row['seq'][args.a1s:args.a1s+args.l].ljust(args.l, 'G')
        r2_seq = reverse_comp(row['seq'])[args.a2s:args.a2s+args.l].ljust(args.l, 'G')
        score = 'I' * args.l
        for count in range(row['counts']):
            r1.write(f'@{id}:{idx}:{row[5]}:{row[4]}:{row[6]}:{row[7]}:' \
                    f'{gen_name} 1\n{r1_seq}\n+\n{score}\n')
            r2.write(f'@{id}:{idx}:{row[5]}:{row[4]}:{row[6]}:{row[7]}:' \
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

    args.f = args.f if args.f else (args.mean + (6*args.sd))
    frag_len = args.f

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

    print('\nsimulating restriction digest\n')
    digest_file = digest_genome(motif_dt, frag_len, proj)
    save_histogram(proj, digest_file)
    print('\ncalculating per base weights\n')
    add_position_weights(digest_file)
    print('\nsimulating size selection\n')
    size_selection(proj, digest_file)

