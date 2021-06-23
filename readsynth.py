#!/usr/bin/env python

import argparse
import math
import matplotlib as plt
import numpy as np
import os
import pandas as pd
import random
import re
import sys

#TODO consider removing l argument, as q scores are necessary
#TODO automate detection of a1s/a2s? (find common adapter? assumes N composition of pos1...)
#TODO check if genome compressed
#TODO check if adapters contain RE motif offsets

def parse_user_input():
    parser = argparse.ArgumentParser(description='simulate RAD libary')

    parser.add_argument('-genome', type=str, required=True,
            help='path to file genome')

    parser.add_argument('-o', type=str, required=True,
            help='path to store output')

    parser.add_argument('-m1', type=str, required=True, nargs='+',
            help='space separated list of RE motifs (e.g., AluI = AG/CT, HindIII = A/AGCTT, SmlI = C/TYRAG)')

    parser.add_argument('-m2', type=str, required=False, nargs='+',
            help='space separated list of RE motifs (e.g., AluI = AG/CT, HindIII = A/AGCTT, SmlI = C/TYRAG)')

    parser.add_argument('-l', type=int, required=False,
            help='desired read length of final simulated reads (defaults to 250 or given q1/q2 profiles)')

    parser.add_argument('-t', type=str, required=False,
            help='test mode: create newline-separated file of RE digested sequences only')

    parser.add_argument('-n', type=int, required=True,
            help='number of reads to simulate (currently per chromosome)')

    parser.add_argument('-mean', type=int, required=True,
            help='mean (in bp) of read lengths after size selection')

    parser.add_argument('-sd', type=int, required=True,
            help='standard deviation (in bp) of read lengths after size selection')

    parser.add_argument('-f', type=int, required=False,
            help='max fragment length after first cut (optional, defaults to mean + 6 stdevs)')

    parser.add_argument('-complete', type=int, required=False,
            help='use \'1\' for complete digestion of fragments (fragments will not contain internal RE sites)')

    parser.add_argument('-a1', type=str, required=True,
            help='file containing tab/space-separated adapters and barcode that attach 5\' to read')

    parser.add_argument('-a2', type=str, required=True,
            help='file containing tab/space-separated adapters and barcode that attach 3\' to read')

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
    '''
    open adapters file and store adapters and barcodes in dt
    '''
    adapters_dt = {}
    with open(arg) as f:
        for line in f:
            adapter, id = line.rstrip().split()
            adapters_dt[adapter] = id

    sbs_start = get_sbs_start(list(adapters_dt.keys()))

    return adapters_dt, sbs_start


def get_sbs_start(adapter_ls):
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


'''
begin main simulation
'''


def readsynth_main(motif_dt, frag_len, proj):
    '''
    process fastq sequences as adapter-ligated fragments
    '''
    with open(args.genome) as fasta,\
         open(os.path.join(proj, 'simulated_R1.fastq'), 'w') as r1,\
         open(os.path.join(proj, 'simulated_R2.fastq'), 'w') as r2:
        seq = ''
        for line in fasta:
            if line.startswith('>') and seq:
                seq_ls = digest_seq(seq, motif_dt, frag_len)
                seq_ls = second_digest(seq_ls, motif_dt, frag_len)
                #TODO add simulate_adapters/simulate_lenght at the end
                #seq_ls, len_ls = simulate_adapters(seq_ls)
                #seq_ls = simulate_length(seq_ls, len_ls, chr_name, proj)

                #if args.r1:
                #    read_writer_samples(seq_ls, r1, r2, chr_name)
                #else:
                #    read_writer(seq_ls, r1, r2, chr_name)

                seq = ''
                chr_name = line.rstrip()[1:].replace(' ', '_')
                print(f'processing {chr_name}')
            elif line.startswith('>'):
                chr_name = line.rstrip()[1:].replace(' ', '_')
                print(f'processing {chr_name}')
            else:
                seq += line.rstrip().upper()

        seq_ls = digest_seq(seq, motif_dt, frag_len)
        seq_ls = second_digest(seq_ls, motif_dt, frag_len)
        #TODO add simulate_adapters/simulate_lenght at the end
        #seq_ls, len_ls = simulate_adapters(seq_ls)
        #seq_ls = simulate_length(seq_ls, len_ls, chr_name, proj)

        #if args.r1:
        #    read_writer_samples(seq_ls, r1, r2, chr_name)
        #else:
        #    read_writer(seq_ls, r1, r2, chr_name)


def digest_seq(seq, motif_dt, frag_len):
    '''
    for every chromosome (seq), find all RE recognition positions
    and preserve frag_len bp ahead as a possible template (fragment)

    if args.m2, preserves motif cut site
    '''
    seq_ls = []
    for motif, offset in motif_dt.items():
        offset = 0 if args.m2 else offset
        for idx in re.finditer(motif, seq):
            start = idx.start()+offset
            end = start + frag_len
            fragment = seq[start:end]
            if fragment not in motif_dt:
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
    #TODO add test for empty lists
    second_ls = []
    for seq in seq_ls:
        second_ls += digest_seq(seq, motif_dt, frag_len)
    second_ls += [reverse_comp(seq) for seq in second_ls]

    if args.m2:
        second_ls = [orientation_test(i) for i in second_ls]
        second_ls = [i for i in second_ls if i]

    if args.complete == 1:
        second_ls = [complete_digest(i, motif_dt) for i in second_ls]
        second_ls = [i for i in second_ls if i]

    if len(second_ls) == 0:
        return second_ls

    return second_ls


def orientation_test(seq):
    '''
    run only if m2 defined
    reads passed here contain intact RE motifs
    (i.e., AG/CT remains AGCT, not CT)

    this is necessary if sticky ends will be used in ligating
    fragments to adapters
    '''
    for motif1, offset1 in motif_dt1.items():
        if re.search('^'+motif1, seq):
            for motif2, offset2 in motif_dt2.items():
                if re.search(motif2+'$', seq):
                    if offset2 == 0:
                        return seq[offset1:]
                    else:
                        return seq[offset1:-offset2]

    return None


def complete_digest(seq, motif_dt):
    '''
    remove any fragments containing an internal cut site for any of the
    restriction enzymes used
    '''
    for motif in motif_dt:
        if motif in seq [1:-1]:
            return False

    return seq


def digest_stats(gen_ls):
    '''
    prints a histogram to screen for every input sequence
    '''
    seq_ls = [seq[0] for seq in gen_ls]
    stats_dt = {}
    len_ls = [len(seq) for seq in seq_ls]
    for num in range(args.sd, max(max(len_ls), args.sd)+args.sd, args.sd):
        stats_dt[int(math.ceil(num / 10.0)) * 10] = 0

    print(f'average fragment length: {round(sum(map(len, seq_ls)) / len(seq_ls), 2)}')
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


def simulate_adapters(seq_ls):
    '''
    simulate ligation of adapters to sequences
    adapters must contain RE motif if sticky ends are used
    the associated barcode is saved alongside the read as r1/2_id
    '''
    print('adding random adapters')
    adapt_ls, len_ls = [], []
    for seq in seq_ls:
        r1_dapt, r1_id = random.choice(list(args.a1.items()))
        r2_dapt, r2_id = random.choice(list(args.a2.items()))
        ligated_seq = r1_dapt + seq + r2_dapt
        adapt_ls.append([ligated_seq, r1_id, r2_id])
        len_ls.append(len(seq))

    return adapt_ls, len_ls


def simulate_length(seq_ls, len_ls, chr_name, proj):
    df = pd.DataFrame(seq_ls, columns = ['sequence', 'r1_id', 'r2_id'])
    df['length'] = [len(i[0]) for i in seq_ls]
    df['fragment_length'] = len_ls
    df.sort_values(['length'], ascending=[True], inplace=True)
    df.reset_index(inplace=True, drop=True)

    df.to_csv(os.path.join(proj, 'raw_digest_' + chr_name + '.csv'), index=None, header=False) #TODO add option to write raw digests to file per file

    read_count = len(df) #TODO keep this in case you want to multiply to get coverage
    total_reads = args.n #TODO add argument COVERAGE

    draw_ls = np.random.normal(loc=args.mean,scale=args.sd,size=total_reads)
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

    sampled_df.to_csv(os.path.join(proj, 'sampled_df.csv')) #TODO add option to write raw digests to file

    print(f'the length of sample_seqs is {len(sampled_df)}')
    reads_to_write = sampled_df['counts'].sum()
    print(f'writing {reads_to_write} simulated reads')

    return sampled_df


def read_writer(sampled_df, r1, r2, chr_name):
    '''
    args.a1s is where to begin in the R1 adapter
    args.a2s is where to begin in the R2 adapter
    '''
    id = 0
    for idx, row in sampled_df.iterrows():
        r1_seq = row[0][args.a1s:args.a1s+args.l].ljust(args.l, 'G')
        r2_seq = reverse_comp(row[0])[args.a2s:args.a2s+args.l].ljust(args.l, 'G')
        for count in range(row[5]):
            r1_mut, r1_score = read_mutator(r1_seq, args.prob_mx1, args.q_dt1, args.q_ls1)
            r2_mut, r2_score = read_mutator(r2_seq, args.prob_mx2, args.q_dt2, args.q_ls2)
            r1.write(f'@{id}:{idx}:{row[3]}:{row[4]}:{row[1]}:{row[2]}:{chr_name} 1\n{r1_mut}\n+\n{r1_score}\n')
            r2.write(f'@{id}:{idx}:{row[3]}:{row[4]}:{row[1]}:{row[2]}:{chr_name} 2\n{r2_mut}\n+\n{r2_score}\n')
            id += 1


def read_mutator(seq, prob_mx, scores_dt, q_ls):
    '''
    using pure, per-base probability profile, mutate bases
    '''
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


def read_writer_samples(sampled_df, r1, r2, chr_name):
    '''
    args.a1s is where to begin in the R1 adapter (33 for ours, I think)
    args.a2s is where to begin in the R2 adapter (34 for ours, I think)
    '''
    print(f'writing simulated reads')
    sampled_q1 = list(pd.read_csv(args.r1 + '_sampled_scores.csv', names=['score'])['score'].values)
    sampled_q2 = list(pd.read_csv(args.r2 + '_sampled_scores.csv', names=['score'])['score'].values)
    args.l = len(max(sampled_q1, key=len))
    id = 0
    for idx, row in sampled_df.iterrows():
        r1_seq = row[0][args.a1s:args.a1s+args.l].ljust(args.l, 'G')
        r2_seq = reverse_comp(row[0])[args.a2s:args.a2s+args.l].ljust(args.l, 'G')
        for count in range(row[5]):
            r1_mut, r1_score = read_mutator_samples(r1_seq, args.q_dt1, sampled_q1)
            r2_mut, r2_score = read_mutator_samples(r2_seq, args.q_dt2, sampled_q2)
            r1.write(f'@{id}:{idx}:{row[3]}:{row[4]}:{row[1]}:{row[2]}:{chr_name} 1\n{r1_mut}\n+\n{r1_score}\n')
            r2.write(f'@{id}:{idx}:{row[3]}:{row[4]}:{row[1]}:{row[2]}:{chr_name} 2\n{r2_mut}\n+\n{r2_score}\n')
            id += 1


def read_mutator_samples(seq, scores_dt, sampled_q):
    '''
    using actual sampled per-read q scores, mutate bases
    '''
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


'''
begin simulation of restriction digests only
'''


def readsynth_digest_only(motif_dt, frag_len, proj):
    '''
    process fastq sequences for digestion simulation only
    '''
    #TODO make fasta fasta_file pulled from directory walk
    digest_file = os.path.join(proj, 'raw_digest_' + os.path.basename(args.genome) + '.csv')
    with open(args.genome) as fasta, open(digest_file, 'w') as dfile:
        gen_ls = []
        seq = ''
        for line in fasta:
            if line.startswith('>') and seq:
                seq_ls = digest_seq(seq, motif_dt, frag_len)
                seq_ls = second_digest(seq_ls, motif_dt, frag_len)
                gen_ls.extend([(i, len(i), chr_name) for i in seq_ls])
                seq = ''
                chr_name = line.rstrip()[1:].replace(' ', '_').replace(',','')
            elif line.startswith('>'):
                chr_name = line.rstrip()[1:].replace(' ', '_').replace(',','')
            else:
                seq += line.rstrip().upper()

        seq_ls = digest_seq(seq, motif_dt, frag_len)
        seq_ls = second_digest(seq_ls, motif_dt, frag_len)
        gen_ls.extend([(i, len(i), chr_name) for i in seq_ls])

        if len(gen_ls) > 0:
            digest_stats(gen_ls)
            write_digested_reads(gen_ls, dfile)

    save_histogram(proj, digest_file)


    #TODO make this a separate function
    #with open(os.path.join(proj, 'digested_R1.txt'), 'w') as r1,\
    #     open(os.path.join(proj, 'digested_R2.txt'), 'w') as r2:
    #        seq_ls = basic_simulate_length(digest_file, proj)
    #        basic_writer(seq_ls, r1, r2)


def write_digested_reads(gen_ls, dfile):
    for seq in gen_ls:
        dfile.write(f'{seq[0]},{str(seq[1])},{seq[2]}\n')


def save_histogram(proj, digest_file):
    df = pd.read_csv(digest_file, names=['sequence', 'length', 'chr'])
    ax = df['length'].hist(bins=100)
    fig = ax.get_figure()
    fig.savefig(os.path.join(proj, 'hist_' + os.path.basename(digest_file)[:-4] + '.pdf'))


def basic_simulate_length(dfile, proj):
    df = pd.read_csv(dfile, names=['sequence', 'length', 'chr'])
    df.sort_values(['length'], ascending=[True], inplace=True)
    df.reset_index(inplace=True, drop=True)
    print(df)

    read_count = len(df) #TODO keep this in case you want to multiply to get coverage
    total_reads = args.n #TODO add argument COVERAGE

    draw_ls = np.random.normal(loc=args.mean,scale=args.sd,size=total_reads)
    draw_ls = [round(i) for i in draw_ls]
    draw_dt = {}

    for i in range(max(min(df['length']), min(draw_ls)), min(max(df['length']), max(draw_ls))+1):
        draw_dt[i] = draw_ls.count(i)

    sampled_df = pd.DataFrame(columns=['sequence', 'length'])
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
    sampled_df.to_csv(os.path.join(proj, 'sampled_df.csv'))
    print(sampled_df)
    return sampled_df


def basic_writer(seq_ls, r1, r2):
    for r1_seq in seq_ls['sequence'].values:
        r2_seq = reverse_comp(r1_seq)
        r1.write(r1_seq + '\n')
        r2.write(r2_seq + '\n')


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

    if args.m2:
        motif_dt2 = iupac_motifs(args.m2)
        motif_dt.update(motif_dt2)

    frag_len = args.f if args.f else (args.mean + (6*args.sd))

    if args.t:
        readsynth_digest_only(motif_dt, frag_len, proj)
        sys.exit()

    if args.r1 and not args.q1:
        sys.exit('please provide q scores profile for R1')

    args.l = args.l if args.l else 250

    if args.a1:
        args.a1, args.a1s = get_adapters(args.a1)

    if args.a2:
        args.a2, args.a2s = get_adapters(args.a2)
        args.a2 = {reverse_comp(adapter): id for adapter, id in args.a2.items()}

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

    readsynth_main(motif_dt, frag_len, proj)


