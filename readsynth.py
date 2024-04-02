#!/usr/bin/env python
"""
File: readsynth.py
Author: Ryan Kuster
Date: March 4, 2024
Description:
  Main script in readsynth. Imports necessary helper scripts for:
    - simulated digestion (ddradseq and iso approaches)
    - per-fragment probability estimation (ddradseq and iso approaches)
    - size selection simulation
License: Apache-2.0 license
"""

import argparse
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import re
import seaborn as sns
import subprocess
import sys

import scripts.digest_genomes as digest_genomes
import scripts.prob_n_copies as prob_n_copies
import scripts.size_selection as size_selection
import scripts.write_reads as write_reads

import scripts.digest_genomes_iso as digest_genomes_iso
import scripts.prob_n_copies_iso as prob_n_copies_iso


def parse_user_input():
    parser = argparse.ArgumentParser(description='simulate RAD libary')

    parser.add_argument('-g', type=str, required=True,
                        help='path to file genome')

    parser.add_argument('-o', type=str, required=True,
                        help='path to store output')

    parser.add_argument('-m1', type=str, required='-iso' not in sys.argv, nargs='+',
                        help='space separated list of search motifs (e.g., RE motifs AluI or AG/CT, or 16S primer /CCTACGGGNGGCWGCAG)')

    parser.add_argument('-m2', type=str, required='-iso' not in sys.argv, nargs='+',
                        help='space separated list of search motifs (e.g., RE motifs SmlI or C/TYRAG, or 16S primer /GACTACHVGGGTATCTAANCC)')

    parser.add_argument('-iso', type=str, required=False,
                        help='optional type IIB RE motif (e.g., BcgI or NN/NNNNNNNNNNCGANNNNNNTGCNNNNNNNNNNNN/)')

    parser.add_argument('-l', '-l1', type=int, required=True,
                        help='desired R1 read length of final simulated reads')

    parser.add_argument('-n', type=int, required=True,
                        help='total read number')

    parser.add_argument('-u', type=int, required='-sd' in sys.argv,
                        help='mean (in bp) of read lengths after size selection')

    parser.add_argument('-sd', type=int, required='-u' in sys.argv,
                        help='standard deviation (in bp) of read lengths after size selection')

    parser.add_argument('-x', type=int, required=False,
                        help='fragment length where fragment distribution intersects size distribution')

    parser.add_argument('-d', type=str, required=False,
                        help='json dictionary of fragment length:count for all expected bp fragments range')

    parser.add_argument('-lp', type=int, required=False,
                        help='low-pass mode: defines maximum expected fragment size, distribution free')

    parser.add_argument('-c', type=float, required=False,
                        help='optional: percent probability of per-site cut; default 1 for complete digestion of fragments (fragments will not contain internal RE sites)')

    parser.add_argument('-a1', type=str, required=False,
                        help='optional: file containing tab/space-separated adapters and barcode that attach 5\' to read')

    parser.add_argument('-a2', type=str, required=False,
                        help='optional: file containing tab/space-separated adapters and barcode that attach 3\' to read')

    parser.add_argument('-a1s', type=int, required=False,
                        help='optional: manually provide bp length of adapter a1 before SBS begins')

    parser.add_argument('-a2s', type=int, required=False,
                        help='optional: manually provide bp length of adapter a1 before SBS begins')

    parser.add_argument('-q1', type=str, required=False,
                        help='optional: file containing newline-separated R1 Q scores >= length -l')

    parser.add_argument('-q2', type=str, required=False,
                        help='optional: file containing newline-separated R2 Q scores >= length -l')

    parser.add_argument('-e', type=str, required=False,
                        help='optional: filler base to use if full adapter contaminaton occurs')

    parser.add_argument('-l2', type=int, required=False,
                        help='optional: desired R2 read length of final simulated reads')

    parser.add_argument('-test', dest='test', action='store_true',
                        help='test mode: skip writing simulated fastq files')

    args = parser.parse_args()

    return args


def open_enzyme_file(args):
    """
    open pickle file of common, named RE motifs
    return dictionary of enzyme names : motifs
    """
    with open(os.path.join(os.path.dirname(__file__),
              args.enzyme_file), 'rb') as type_iip_file:
        re_dt = pickle.load(type_iip_file)

    return re_dt


def check_for_enzymes(args, re_dt):
    """
    check if user input is named restriction enzyme
    return either dictionary value or unmodified user input (m1, m2)
    """
    m1 = [re_dt[i.lower()] if i.lower() in re_dt.keys() else i for i in args.m1]
    m2 = [re_dt[i.lower()] if i.lower() in re_dt.keys() else i for i in args.m2]

    return m1, m2


def check_for_enzymes_iso(args, re_dt):
    """
    check if user input is named restriction enzyme
    return either dictionary value or unmodified user input (m1)
    """
    if args.iso.lower() in re_dt.keys():
        m1 = re_dt[args.iso.lower()]
    else:
        m1 = args.iso

    return m1


def check_for_palindrome(arg_m):
    """
    most RE motifs are palindromic, but for those that are not, apply
    a cut position in the reverse complement at the corresponding site
    as the input and return a list of sequence(s) with cut index
    per motif:
      - find the cut position ("/") index
      - get the reverse complement of the input sequence
      - check if reverse complement = input (palindromic)
      - if not palindromic, 
    """
    tmp_ls = []
    for motif in arg_m:
        cut_pos = motif.index('/')
        orig = motif.replace('/', '')
        revc_orig = reverse_comp(orig)
        if orig == revc_orig:
            pass
        else:
            revc_orig = revc_orig[:cut_pos] + '/' + revc_orig[cut_pos:]
            tmp_ls.append(revc_orig)
    arg_m += tmp_ls

    return arg_m


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


def iupac_motifs_iso(arg_m):
    '''
    given a single iso-length cut motif, return a dictionary of regex
    compatible iupac redundancy codes as keys and a list of the two
    cleaving sites as the value

    the reverse complement is considered as these are not palindromic
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

    arg_m.extend([reverse_comp(i) for i in arg_m])

    for motif in arg_m:
        reg_motif = ''
        for char in motif.upper():
            reg_motif += iupac_dt[char]
        motif_dt[reg_motif] = [m.start() for m in re.finditer('/', motif)]
        motif_dt[reg_motif][-1] = motif_dt[reg_motif][-1] - 1

    return motif_dt


def get_motif_regex_len(args):
    """
    returns a dictionary (motif_len) of RE motifs and their final
    length (adjusting for redundant IUPAC bases)
    """
    motif_len = {}
    for motif in args.motif_dt.keys():
        mot_len, count = 0, True
        for j in motif:
            if j == '[':
                count = False
            elif j == ']':
                count = True
                mot_len += 1
            else:
                if count is True:
                    mot_len += 1
        motif_len[motif] = mot_len

    return  motif_len


def check_custom_distribution(args):
    """
    if user provides a custom json distribution of fragment lengths,
    load file and create dictionary (confirm data format works)
    return the maximum fragment size for use in simulations
    """
    with open(args.d) as f_o:
        tmp_dt = json.load(f_o)
    tmp_dt = {int(k): int(v) for k, v in tmp_dt.items()}

    return max(list(tmp_dt.keys()))


def get_adapters(adapter_file):
    """
    open adapters file and store adapters and barcodes in list of tuples
    """
    adapters_ls = []
    with open(adapter_file) as f:
        for line in f:
            a_top, a_bot, a_id = line.rstrip().split()
            adapters_ls.append((a_top, a_bot, a_id))

    return adapters_ls


def create_adapters(args):
    """
    default adapters used for simulating reads where fragment less than
    input read length
    returns list of [adapter, adapter name] where sequence is modified
    with the expected overhang based on the RE motifs
    """
    a1 = ['AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG',
          'CTGTCTCTTATACACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
          'rs1']
    a2 = ['CAAGCAGAAGACGGCATACGAGATGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG',
          'CTGTCTCTTATACACATCTCCGAGCCCACGAGACATCTCGTATGCCGTCTTCTGCTTG',
          'rs2']
    m1 = list(args.motif_dt1.keys())[0]
    m2 = list(args.motif_dt2.keys())[0]
    a1[0] = a1[0] + args.m1[0][:args.motif_dt1[m1]]
    a1[1] = args.m1[0][args.motif_dt1[m1]+1:] + a1[1]
    a2[0] = a2[0] + args.m2[0][:args.motif_dt2[m2]]
    a2[1] = args.m2[0][args.motif_dt2[m2]+1:] + a2[1]

    return [a1], [a2]


def create_adapters_iso(args):
    """
    default adapters used for simulating reads where fragment less than
    input read length
    returns list of [adapter, adapter name] where sequence is modified
    with the expected overhang based on the iso-length RE motifs

    note: the overhang behavior differs from standard RE motifs and is a
    fill of "N" for the distance to the cut site
    """
    a1 = ['AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG',
          'CTGTCTCTTATACACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
          'rs1']
    a2 = ['CAAGCAGAAGACGGCATACGAGATGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG',
          'CTGTCTCTTATACACATCTCCGAGCCCACGAGACATCTCGTATGCCGTCTTCTGCTTG',
          'rs2']
    m1 = list(args.motif_dt.keys())[0]
    m2 = list(args.motif_dt.keys())[1]
    a1[0] = a1[0] + 'N' * args.motif_dt[m1][0]
    a1[1] = 'N' * args.motif_dt[m2][0] + a1[1]
    a2[0] = a2[0] + 'N' * args.motif_dt[m1][0]
    a2[1] = 'N' * args.motif_dt[m2][0] + a2[1]

    return [a1], [a2]


def check_genomes(genome_file):
    '''
    open 'genome_file' abundance profile
    assert each fasta file exists as listed
    '''
    df = pd.read_csv(genome_file, names=['genome', 'abundance'])
    df['abundance'] = df['abundance'] / df['abundance'].sum()

    for idx in range(df.shape[0]):
        genome = df.iloc[idx]['genome']
        assert os.path.exists(genome), f'path to {genome} not found'

    return df


def process_genomes(args, genomes_df):
    '''
    len_freqs is a dictionary where each key is a fragment length from
    a digested genome, the value is the sum of all fragment probabilities
    for that length after adjusting for composition and size selection

    total_freqs collects the frequency for fragment lengths from all the
    genomes to be processed
    '''
    digest_ls, prob_ls = [], []
    total_freqs = pd.DataFrame({'length': pd.Series(dtype='int64'),
                                'sum_prob': pd.Series(dtype='float64'),
                                'name': pd.Series(dtype='object'),
                                'counts_file': pd.Series(dtype='object')})

    for idx in range(genomes_df.shape[0]):
        args.g = genomes_df.iloc[idx]['genome']
        print(args.g)
        args.comp = genomes_df.iloc[idx]['abundance']
        digest_file = os.path.join(args.o, 'raw_digests', 'raw_digest_' +
                                   os.path.basename(args.g) + '.csv')

        df = digest_genomes.main(args)

        if df.shape[0] == 0:
            digest_ls.append(None)
            prob_ls.append(None)
            continue

        digest_file = process_df(df, digest_file, args)

        if digest_file is None:
            digest_ls.append(None)
            prob_ls.append(None)
            continue

        prob_file, len_freqs = prob_n_copies.main(digest_file, args)
        save_individual_hist(prob_file, args)
        digest_ls.append(digest_file)
        prob_ls.append(prob_file)
        tmp_df = pd.DataFrame(len_freqs.items(), columns=['length', 'sum_prob'])
        tmp_df['name'] = os.path.basename(args.g)
        tmp_df['counts_file'] = prob_file
        if tmp_df.empty is False:
            total_freqs = pd.concat([total_freqs, tmp_df], axis=0)

    total_freqs = total_freqs.reset_index(drop=True)
    genomes_df['digest_file'] = digest_ls
    genomes_df['prob_file'] = prob_ls

    return genomes_df, total_freqs


def process_genomes_iso(args, genomes_df):
    '''
    len_freqs is a dictionary where each key is a fragment length from
    a digested genome, the value is the sum of all fragment probabilities
    for that length after adjusting for composition and size selection

    total_freqs collects the frequency for fragment lengths from all the
    genomes to be processed
    '''
    digest_ls, prob_ls = [], []
    total_freqs = pd.DataFrame({'length': pd.Series(dtype='int64'),
                                'sum_prob': pd.Series(dtype='float64'),
                                'name': pd.Series(dtype='object'),
                                'counts_file': pd.Series(dtype='object')})

    for idx in range(genomes_df.shape[0]):
        args.g = genomes_df.iloc[idx]['genome']
        print(args.g)
        args.comp = genomes_df.iloc[idx]['abundance']
        digest_file = os.path.join(args.o, 'raw_digests', 'raw_digest_' +
                                   os.path.basename(args.g) + '.csv')

        df = digest_genomes_iso.main(args)

        if df.shape[0] == 0:
            digest_ls.append(None)
            prob_ls.append(None)
            continue

        digest_file = process_df_iso(df, digest_file, args)
        prob_file, len_freqs = prob_n_copies_iso.main(digest_file, args)
        digest_ls.append(digest_file)
        prob_ls.append(prob_file)
        tmp_df = pd.DataFrame(len_freqs.items(), columns=['length', 'sum_prob'])
        tmp_df['name'] = os.path.basename(args.g)
        tmp_df['counts_file'] = prob_file
        if tmp_df.empty is False:
            total_freqs = pd.concat([total_freqs, tmp_df], axis=0)

    total_freqs = total_freqs.reset_index(drop=True)
    genomes_df['digest_file'] = digest_ls
    genomes_df['prob_file'] = prob_ls

    return genomes_df, total_freqs


def reverse_comp(seq):
    '''
    return the reverse complement of an input sequence
    '''
    revc = {'/': '/',
            'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
            'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
            'N': 'N'}
    new = ''
    for base in reversed(seq):
        new += revc[base]

    return new


def process_df(df, digest_file, args):
    """
    upon completion of simulated digest:
      - filter the pandas dataframe keeping only viable combinations of
        REs that will produce a library
      - replace redundant bases with N
      - account for viable reverse complements
      - trim over/underhang from the produced fragments
      - save as dataframe
    """

    df['forward'] = np.where((df['m1'].isin(args.motif_dt1.keys()) &
                             df['m2'].isin(args.motif_dt2.keys())), 1, 0)
    df['reverse'] = np.where((df['m1'].isin(args.motif_dt2.keys()) &
                             df['m2'].isin(args.motif_dt1.keys())), 1, 0)

    """
    remove unviable combos after getting site_ls
    """
    df.drop(df[(df['forward'] == 0) & (df['reverse'] == 0)].index,
            inplace=True)
    df = df.reset_index(drop=True)

    """
    remove fragments where RE motif sites overlap
    """
    df['min_len_m1'] = df['m1'].map(args.motif_len)
    df['min_len_m2'] = df['m2'].map(args.motif_len)
    df['min_len'] = df['min_len_m1'] + df['min_len_m2']
    df = df[df['seq'].str.len() > df['min_len']]
    df.drop('min_len_m1', axis=1, inplace=True)
    df.drop('min_len_m2', axis=1, inplace=True)
    df.drop('min_len', axis=1, inplace=True)

    """
    convert all redundant IUPAC codes to 'N'
    """
    df['seq'] = df['seq'].str.replace('[RYSWKMBDHV]', 'N', regex=True)

    """
    create a column of reverse complement sequences
    """
    df['revc'] = [reverse_comp(i) for i in df['seq'].to_list()]

    """
    use tmp_df to temporarilty store bidirectional reads
    duplicate all the reads that work both ways (if RE in m1 and m2)
    """
    tmp_df = df[(df['forward'] == 1) & (df['reverse'] == 1)]
    tmp_df = tmp_df.reset_index(drop=True)

    """
    recategorize the bidirectional seqs in df as being forward
    """
    df.loc[(df['forward'] == 1) & (df['reverse'] == 1), 'reverse'] = 0

    """
    convert unidirectional reverse strand sequences to the reverse complement
    """
    df.loc[df['reverse'] == 1, 'tmp_seq'] = df['revc']
    df.loc[df['reverse'] == 1, 'revc'] = df['seq']
    df.loc[df['reverse'] == 1, 'seq'] = df['tmp_seq']
    df.drop('tmp_seq', axis=1, inplace=True)

    """
    make all tmp_df reads reverse complements and recategorize as reverse
    """
    tmp_df['tmp_seq'] = tmp_df['revc'].values
    tmp_df['revc'] = tmp_df['seq'].values
    tmp_df['seq'] = tmp_df['tmp_seq'].values
    tmp_df.drop('tmp_seq', axis=1, inplace=True)
    tmp_df.loc[:, 'forward'] = 0
    tmp_df.loc[:, 'reverse'] = 1

    df = pd.concat([df, tmp_df])
    df.drop('forward', axis=1, inplace=True)
    df = df.reset_index(drop=True)

    """
    add a quick step that removes appropriate over/underhang
    """

    for mot, front in args.motif_dt.items():
        back = args.motif_len[mot] - front
        df.loc[(df['m1'] == mot) & (df['reverse'] == 0), 'seq'] = \
            df['seq'].str[front:]
        if back != 0:
            df.loc[(df['m1'] == mot) & (df['reverse'] == 0), 'revc'] = \
                df['revc'].str[:-back]

        df.loc[(df['m1'] == mot) & (df['reverse'] == 1), 'seq'] = \
            df['seq'].str[:-back]
        df.loc[(df['m1'] == mot) & (df['reverse'] == 1), 'revc'] = \
            df['revc'].str[front:]

        df.loc[(df['m2'] == mot) & (df['reverse'] == 0), 'seq'] = \
            df['seq'].str[:-back]
        df.loc[(df['m2'] == mot) & (df['reverse'] == 0), 'revc'] = \
            df['revc'].str[front:]

        df.loc[(df['m2'] == mot) & (df['reverse'] == 1), 'seq'] = \
            df['seq'].str[front:]
        if back != 0:
            df.loc[(df['m2'] == mot) & (df['reverse'] == 1), 'revc'] = \
                df['revc'].str[:-back]


    df['length'] = df['seq'].str.len()
    #df = df[(df['seq'].str.len() > 0) & (df['revc'].str.len() > 0)]

    df = df.sort_values(by=['length'])
    df = df.reset_index(drop=True)
    df.to_csv(digest_file, index=None)

    if df.shape[0] == 0:
        digest_file = None

    return digest_file


def process_df_iso(df, digest_file, args):
    """
    upon completion of simulated digest:
      - replace redundant bases with N
      - account for viable reverse complements
      - trim over/underhang from the produced fragments
      - save as dataframe
    """

    df['forward'] = np.where(df['m1'] == list(args.motif_dt.keys())[0], 1, 0)
    df['reverse'] = np.where(df['m1'] == list(args.motif_dt.keys())[1], 1, 0)

    """
    convert all redundant IUPAC codes to 'N'
    """
    df['seq'] = df['seq'].str.replace('[RYSWKMBDHV]', 'N', regex=True)

    """
    create a column of reverse complement sequences
    """
    df['revc'] = [reverse_comp(i) for i in df['seq'].to_list()]

    """
    swap seq and revc for fragments on the reverse direction
    """
    tmp_df = df['reverse'] == 1
    df.loc[tmp_df, ['seq', 'revc']] = (df.loc[tmp_df, ['revc', 'seq']].values)

    df.drop('forward', axis=1, inplace=True)
    df = df.reset_index(drop=True)

    """
    add a quick step that removes appropriate over/underhang
    """
    m1 = list(args.motif_dt.keys())[0]
    m1_f = args.motif_dt[m1][0]
    m1_b = args.motif_dt[m1][1]

    df['seq'] = df['seq'].str[m1_f:m1_b]
    df['revc'] = df['revc'].str[m1_f:m1_b]

    df['length'] = df['seq'].str.len()
    df = df.sort_values(by=['start'])
    df = df.reset_index(drop=True)
    df.to_csv(digest_file, index=None)

    return digest_file


def save_individual_hist(prob_file, args):
    """
    per genome in the sample, create a distribution of frament lengths:
      - red: all fragments (raw counts)
      - gold: fragments after applying frequency/size distributin
      - blue: frequency/size adjusted fragments after relative abundance
    """
    df = pd.read_csv(prob_file)
    if df.shape[0] == 0:
        print(f'no fragments found in {prob_file}')
        return

    try:
        sns.histplot(data=df,
                     x=df['length'],
                     binwidth=6,
                     alpha=0.75,
                     color='red')
        sns.histplot(data=df,
                     x=df['length'],
                     weights=df['probability'],
                     binwidth=6,
                     alpha=0.75,
                     color='gold')
        sns.histplot(data=df,
                     x=df['length'],
                     weights=df['adj_prob'],
                     binwidth=6,
                     alpha=0.75,
                     color='blue')
        plt.savefig(os.path.join(args.o, 'individual_histograms', 'hist_' +
                    os.path.basename(prob_file)[:-4] + '.png'),
                    facecolor='white', transparent=False)
        plt.close()

    except ValueError:
        #print(f'too few bins to produce histogram for {prob_file}') #LOG
        return


def save_combined_hist(total_freqs, image_name, weights, args):
    """
    stacked histogram of all members of community
    """
    try:
        ax = sns.histplot(data=total_freqs, x='length', hue='name',
                          weights=total_freqs[weights], multiple="stack",
                          binwidth=6, element="step")
    except IndexError:
        print('singular read lengths, cannot produce histogram')
        return

    old_legend = ax.legend_
    handles = old_legend.legend_handles
    labels = [t.get_text() for t in old_legend.get_texts()]
    ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left',
              borderaxespad=0)
    plt.savefig(os.path.join(args.o, f'{image_name}.pdf'),
                bbox_inches='tight')
    plt.close()


def prob_to_counts(comb_file, fragment_comps, adjustment, genomes_df):
    """
    after calculated per-fragment adjusted probabilities, use the pre-
    determined "adjustment" constant to update the "adj_prob" column in
    each digested genome to rounded counts of reads based on the desired
    library size

    save the dataframes with "counts" column
    """
    comb_df = pd.read_csv(comb_file)
    count_files_ls = list(set(comb_df['counts_file'].to_list()))
    total = 0

    basenames = [os.path.basename(i) for i in genomes_df['genome']]
    genomes_df['genome'] = basenames
    genomes_df['reads'] = np.nan
    genomes_df['sites'] = np.nan

    for count_file in count_files_ls:
        df = pd.read_csv(count_file)
        df['counts'] = df['length'].map(fragment_comps)
        df['counts'] = round(df['counts'] * df['adj_prob'] * adjustment)
        df.dropna(subset=['counts'], inplace=True)
        df['counts'] = df['counts'].astype(int)
        df = df[df['counts'] > 0]
        total += df['counts'].sum()
        df = df.reset_index(drop=True)
        df.to_csv(count_file)

        b_name = os.path.basename(count_file)[7:-4]
        genomes_df.loc[genomes_df.genome == b_name, 'reads'] = df['counts'].sum()
        genomes_df.loc[genomes_df.genome == b_name, 'sites'] = df.shape[0]

    return genomes_df


def write_final_file(args, genomes_df):
    """
    calculate simple summary statistics for genomes_df and write to file
    """
    genomes_df['avg_depth'] = genomes_df['reads'] / genomes_df['sites']
    genomes_df['depth_abundance'] = genomes_df['avg_depth'] / \
        genomes_df['avg_depth'].sum()
    genomes_df['read_abundance'] = genomes_df['reads'] / \
        genomes_df['reads'].sum()
    genomes_df.to_csv(os.path.join(args.o, 'metagenome_summary.csv'))


def write_genomes(comb_file, fragment_comps, adjustment):
    """
    write probability and size-selection adjusted fragment counts as
    adapter-ligated, simulated fastq paired-end read files (sim1/sim2)

    each count file (as listed in the comb_file) is opened as a pandas
    dataframe and passed to the write_reads script
    """
    comb_df = pd.read_csv(comb_file)
    count_files_ls = list(set(comb_df['counts_file'].to_list()))

    sim1 = os.path.join(args.o, 'sim_metagenome_R1.fastq')
    sim2 = os.path.join(args.o, 'sim_metagenome_R2.fastq')
    error1 = os.path.join(args.o, 'error_sim_metagenome_R1.fastq')
    error2 = os.path.join(args.o, 'error_sim_metagenome_R2.fastq')

    with open(sim1, 'w') as r1, open(sim2, 'w') as r2:
        for count_file in count_files_ls:
            gen_name = os.path.basename(count_file)[7:-4]
            df = pd.read_csv(count_file)
            df = df[['seq', 'revc', 'length', 'counts']]
            write_reads.main(df, r1, r2, gen_name, args)

    if args.q1 and args.q2:
        print('applying error profile')
        command = os.path.join(os.path.dirname(__file__), "src", "apply_error")
        simulate_error(command, sim1, error1)
        os.remove(sim1)
        simulate_error(command, sim2, error2)
        os.remove(sim2)


def simulate_error(command, sim_in, error_out):
    """
    if selected, calls apply_error.cpp (from src directory) to directly
    mutate single bases based on the input phred likelihood of miscall
    stored in the q scores
    """
    try:
        process = subprocess.Popen([command, sim_in, error_out], shell=False)
        out, err = process.communicate()
        errcode = process.returncode
        process.kill()
        process.terminate()
    except FileNotFoundError:
        sys.exit('please run \'make apply_error\' in the src directory of readsynth')


if __name__ == '__main__':
    args = parse_user_input()

    if args.o is None:
        args.o = os.path.dirname(os.path.abspath(__file__))
    elif os.path.exists(args.o) is True:
        args.o = os.path.abspath(args.o)
        os.mkdir(os.path.join(args.o, "raw_digests"))
        os.mkdir(os.path.join(args.o, "individual_histograms"))
        os.mkdir(os.path.join(args.o, "individual_counts"))
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))

    if args.iso:
        args.enzyme_file = 'resources/type_iib_enzymes.pickle'
        re_dt = open_enzyme_file(args)
        args.iso = check_for_enzymes_iso(args, re_dt)
        args.motif_dt = iupac_motifs_iso([args.iso])
    else:
        args.enzyme_file = 'resources/type_iip_enzymes.pickle'
        args.motif_dt = {}
        re_dt = open_enzyme_file(args)
        args.m1, args.m2 = check_for_enzymes(args, re_dt)
        args.m1 = check_for_palindrome(args.m1)
        args.m2 = check_for_palindrome(args.m2)
        args.motif_dt1 = iupac_motifs(args.m1)
        args.motif_dt.update(args.motif_dt1)
        args.motif_dt2 = iupac_motifs(args.m2)
        args.motif_dt.update(args.motif_dt2)

    args.motif_len = get_motif_regex_len(args)

    if args.d:
        args.max = check_custom_distribution(args)
    elif args.lp:
        args.max = args.lp
    elif args.u and args.sd:
        args.max = args.u + (6*args.sd)
    else:
        sys.exit('must define -d, -lp, or -u/-sd')

    if not args.x:
        args.x = args.u

    if args.a1 and args.a2:
        args.a1 = get_adapters(args.a1)
        if not args.a1s:
            args.a1s = len(args.a1[0])
        args.a2 = get_adapters(args.a2)
        if not args.a2s:
            args.a2s = len(args.a2[0])
    elif args.m1:
        args.a1, args.a2 = create_adapters(args)
        args.a1s, args.a2s = 62, 58
    else:
        args.a1, args.a2 = create_adapters_iso(args)
        args.a1s, args.a2s = 62, 58

    args.l1 = args.l

    if not args.l2:
        args.l2 = args.l1

    if args.q1 or args.q2:
        q_src = os.path.join(os.path.dirname(__file__), "src", "apply_error")
        if not os.path.exists(q_src):
            sys.exit('please run \'make apply_error\' in the src directory of'\
                     'readsynth')
        if not args.q1 or not args.q2:
            sys.exit('arguments q1 and q2 required')

    if not args.c:
        args.c = 1

    if args.e:
        args.e = args.e[0]
    else:
        args.e = 'G'

    '''
    1.
    digest genomes one by one, producing raw digest files
    '''
    genomes_df = check_genomes(args.g)

    print('\n1. finding target sequences\n')
    if args.iso:
        genomes_df, total_freqs = process_genomes_iso(args, genomes_df)
    else:
        genomes_df, total_freqs = process_genomes(args, genomes_df)

    if total_freqs['sum_prob'].sum() == 0:
        sys.exit('no fragments produced, exiting')

    save_combined_hist(total_freqs, 'fragment_distributions', 'sum_prob', args)

    '''
    2.
    combine relative probabilities of fragments from all input genome digests
    and perform simulation of pooled size selection
    '''
    print('\n2. simulating size selection\n')
    if args.iso or args.lp:
        fragment_comps = \
            total_freqs.groupby('length')['sum_prob'].apply(list).to_dict()
        fragment_comps = {k: sum(v) for k, v in fragment_comps.items()}
        adjustment = args.n / sum(fragment_comps.values())
        fragment_comps = \
            {k: 1 if v > 0 else 0 for k, v in fragment_comps.items()}
    else:
        fragment_comps, adjustment = size_selection.main(total_freqs, args)

    total_freqs['counts'] = total_freqs['length'].map(fragment_comps)
    total_freqs['counts'] = \
        round(total_freqs['counts'] * total_freqs['sum_prob'] * adjustment)

    comb_file = os.path.join(args.o, 'combined.csv')
    total_freqs.to_csv(comb_file)
    save_combined_hist(total_freqs, 'read_distributions', 'counts', args)
    genomes_df = prob_to_counts(
        comb_file, fragment_comps, adjustment, genomes_df)
    write_final_file(args, genomes_df)

    if args.test:
        sys.exit()

    '''
    3.
    write fragments to fastq-formatted file with adapters concatenated
    '''
    print('\n3. simulating fastq formatted sequence reads\n')
    write_genomes(comb_file, fragment_comps, adjustment)

