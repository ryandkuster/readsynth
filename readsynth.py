#!/usr/bin/env python

import argparse
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

    parser.add_argument('-genome', type=str, required=True,
                        help='path to file genome')

    parser.add_argument('-o', type=str, required=True,
                        help='path to store output')

    parser.add_argument('-m1', type=str, required='-iso' not in sys.argv, nargs='+',
                        help='space separated list of RE motifs (e.g., AluI or AG/CT, HindIII or A/AGCTT, SmlI or C/TYRAG)')

    parser.add_argument('-m2', type=str, required='-iso' not in sys.argv, nargs='+',
                        help='space separated list of RE motifs (e.g., AluI or AG/CT, HindIII or A/AGCTT, SmlI or C/TYRAG)')

    parser.add_argument('-iso', type=str, required=False,
                        help='optional type IIB RE motif (e.g., NN/NNNNNNNNNNCGANNNNNNTGCNNNNNNNNNNNN/)')

    parser.add_argument('-l', type=int, required=True,
                        help='desired read length of final simulated reads')

    parser.add_argument('-test', dest='test', action='store_true',
                        help='test mode: create newline-separated file of RE digested sequences only')

    parser.add_argument('-n', type=int, required=True,
                        help='total read number')

    parser.add_argument('-mean', type=int, required=True,
                        help='mean (in bp) of read lengths after size selection')

    parser.add_argument('-up_bound', type=int, required=True,
                        help='the upper end of a range (in bp) of read lengths to size select')

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
                        help='file containing newline-separated R1 Q scores >= length -l')

    parser.add_argument('-q2', type=str, required=False,
                        help='file containing newline-separated R2 Q scores >= length -l')

    args = parser.parse_args()

    return args


def check_for_enzymes(args):
    with open(os.path.join(os.path.dirname(__file__),
              'resources/type_iip_enzymes.pickle'), 'rb') as type_iip_file:
        re_dt = pickle.load(type_iip_file)

    m1 = [re_dt[i.lower()] if i.lower() in re_dt.keys() else i for i in args.m1]
    m2 = [re_dt[i.lower()] if i.lower() in re_dt.keys() else i for i in args.m2]

    return m1, m2


def check_for_enzymes_iso(args):
    with open(os.path.join(os.path.dirname(__file__),
              'resources/type_iib_enzymes.pickle'), 'rb') as type_iip_file:
        re_dt = pickle.load(type_iip_file)

    if args.iso.lower() in re_dt.keys():
        m1 = re_dt[args.iso.lower()]
    else:
        m1 = args.iso

    return m1


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
    total_freqs = pd.DataFrame(columns=['length',
                                        'sum_prob',
                                        'name',
                                        'counts_file'])

    for idx in range(genomes_df.shape[0]):
        args.genome = genomes_df.iloc[idx]['genome']
        args.comp = genomes_df.iloc[idx]['abundance']
        digest_file = os.path.join(args.o, 'raw_digest_' +
                                   os.path.basename(args.genome) + '.csv')

        print(f'processing genome {os.path.basename(args.genome)}')
        df = digest_genomes.main(args)

        if df.shape[0] == 0:
            digest_ls.append(None)
            prob_ls.append(None)
            print(f'no fragments found in {args.genome}\n')
            continue

        digest_file = process_df(df, digest_file, args)
        prob_file, len_freqs = prob_n_copies.main(digest_file, args)
        save_individual_hist(prob_file, args)
        digest_ls.append(digest_file)
        prob_ls.append(prob_file)
        tmp_df = pd.DataFrame(len_freqs.items(), columns=['length', 'sum_prob'])
        tmp_df['name'] = os.path.basename(args.genome)
        tmp_df['counts_file'] = prob_file
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
    total_freqs = pd.DataFrame(columns=['length', 'sum_prob', 'name', 'counts_file'])

    for idx in range(genomes_df.shape[0]):
        args.genome = genomes_df.iloc[idx]['genome']
        args.comp = genomes_df.iloc[idx]['abundance']
        digest_file = os.path.join(args.o, 'raw_digest_' +
                                   os.path.basename(args.genome) + '.csv')

        print(f'processing genome {os.path.basename(args.genome)}')
        df = digest_genomes_iso.main(args)

        if df.shape[0] == 0:
            digest_ls.append(None)
            prob_ls.append(None)
            print(f'no fragments found in {args.genome}\n')
            continue

        digest_file = process_df_iso(df, digest_file, args)
        prob_file, len_freqs = prob_n_copies_iso.main(digest_file, args)
        digest_ls.append(digest_file)
        prob_ls.append(prob_file)
        tmp_df = pd.DataFrame(len_freqs.items(), columns=['length', 'sum_prob'])
        tmp_df['name'] = os.path.basename(args.genome)
        tmp_df['counts_file'] = prob_file
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
    df['forward'] = np.where((df['m1'].isin(args.motif_dt1.keys()) &
                             df['m2'].isin(args.motif_dt2.keys())), 1, 0)
    df['reverse'] = np.where((df['m1'].isin(args.motif_dt2.keys()) &
                             df['m2'].isin(args.motif_dt1.keys())), 1, 0)

    # remove unviable combos after getting site_ls
    df.drop(df[(df['forward'] == 0) & (df['reverse'] == 0)].index,
            inplace=True)
    df = df.reset_index(drop=True)

    # convert all redundant IUPAC codes to 'N'
    df['seq'] = df['seq'].str.replace('[RYSWKMBDHV]', 'N')

    # create a column of reverse complement sequences
    df['revc'] = [reverse_comp(i) for i in df['seq'].to_list()]

    # use tmp_df to temporarilty store bidirectional reads
    # duplicate all the reads that work both ways (if RE in m1 and m2)
    tmp_df = df[(df['forward'] == 1) & (df['reverse'] == 1)]
    tmp_df = tmp_df.reset_index(drop=True)

    # recategorize the bidirectional seqs in df as being forward
    df.loc[(df['forward'] == 1) & (df['reverse'] == 1), 'reverse'] = 0

    # convert unidirectional reverse strand sequences to the reverse complement
    df.loc[df['reverse'] == 1, 'tmp_seq'] = df['revc']
    df.loc[df['reverse'] == 1, 'revc'] = df['seq']
    df.loc[df['reverse'] == 1, 'seq'] = df['tmp_seq']
    df.drop('tmp_seq', axis=1, inplace=True)

    # make all tmp_df reads reverse complements and recategorize as reverse
    tmp_df['tmp_seq'] = tmp_df['revc'].values
    tmp_df['revc'] = tmp_df['seq'].values
    tmp_df['seq'] = tmp_df['tmp_seq'].values
    tmp_df.drop('tmp_seq', axis=1, inplace=True)
    tmp_df.loc[:, 'forward'] = 0
    tmp_df.loc[:, 'reverse'] = 1

    df = pd.concat([df, tmp_df])
    df.drop('forward', axis=1, inplace=True)
    df = df.reset_index(drop=True)

    # add a quick step that removes appropriate over/underhang
    for mot, front in args.motif_dt.items():
        back = len(mot) - front
        df.loc[(df['m1'] == mot) & (df['reverse'] == 0), 'seq'] = \
            df['seq'].str[front:]
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
        df.loc[(df['m2'] == mot) & (df['reverse'] == 1), 'revc'] = \
            df['revc'].str[:-back]

    df['length'] = df['seq'].str.len()
    df = df.sort_values(by=['length'])
    df = df.reset_index(drop=True)
    df.to_csv(digest_file, index=None)

    return digest_file


def process_df_iso(df, digest_file, args):
    df['forward'] = np.where(df['m1'] == list(args.motif_dt.keys())[0], 1, 0)
    df['reverse'] = np.where(df['m1'] == list(args.motif_dt.keys())[1], 1, 0)

    # convert all redundant IUPAC codes to 'N'
    df['seq'] = df['seq'].str.replace('[RYSWKMBDHV]', 'N')

    # create a column of reverse complement sequences
    df['revc'] = [reverse_comp(i) for i in df['seq'].to_list()]

    # swap seq and revc for fragments on the reverse direction
    tmp_df = df['reverse'] == 1
    df.loc[tmp_df, ['seq', 'revc']] = (df.loc[tmp_df, ['revc', 'seq']].values)

    df.drop('forward', axis=1, inplace=True)
    df = df.reset_index(drop=True)

    # add a quick step that removes appropriate over/underhang
    m1 = list(args.motif_dt.keys())[0]
    m1_f = args.motif_dt[m1][0]
    m1_b = args.motif_dt[m1][1]
    m2 = list(args.motif_dt.keys())[1]
    m2_f = args.motif_dt[m2][0]
    m2_b = args.motif_dt[m2][1]

    df['seq'] = df['seq'].str[m1_f:m1_b]
    df['revc'] = df['revc'].str[m2_f:m2_b]

    df['length'] = df['seq'].str.len()
    df = df.sort_values(by=['start'])
    df = df.reset_index(drop=True)
    df.to_csv(digest_file, index=None)

    return digest_file


def save_individual_hist(prob_file, args):
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
        plt.savefig(os.path.join(args.o, 'hist_' +
                    os.path.basename(prob_file)[:-4] + '.png'),
                    facecolor='white', transparent=False)
        plt.close()

    except ValueError:
        print(f'too few bins to produce histogram for {prob_file}')
        return


def save_combined_hist(total_freqs, image_name, weights, args):
    try:
        ax = sns.histplot(data=total_freqs, x='length', hue='name',
                          weights=total_freqs[weights], multiple="stack",
                          binwidth=6, element="step")
    except IndexError:
        print('singular read lengths, cannot produce histogram')
        return

    old_legend = ax.legend_
    handles = old_legend.legendHandles
    labels = [t.get_text() for t in old_legend.get_texts()]
    ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left',
              borderaxespad=0)
    plt.savefig(os.path.join(args.o, f'_{image_name}.pdf'),
                bbox_inches='tight')


def prob_to_counts(comb_file, fragment_comps, adjustment, genomes_df):
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
    genomes_df['avg_depth'] = genomes_df['reads'] / genomes_df['sites']
    genomes_df['depth_abundance'] = genomes_df['avg_depth'] / \
        genomes_df['avg_depth'].sum()
    genomes_df['read_abundance'] = genomes_df['reads'] / \
        genomes_df['reads'].sum()
    genomes_df.to_csv(os.path.join(args.o, 'metagenome_summary.csv'))


def write_genomes(comb_file, fragment_comps, adjustment):
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
        simulate_error(command, sim2, error2)


def simulate_error(command, sim_in, error_out):
    process = subprocess.Popen([command, sim_in, error_out], shell=False)
    out, err = process.communicate()
    errcode = process.returncode
    process.kill()
    process.terminate()


if __name__ == '__main__':
    args = parse_user_input()

    if args.o is None:
        args.o = os.path.dirname(os.path.abspath(__file__))
    elif os.path.exists(args.o) is True:
        args.o = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))

    if args.iso:
        args.iso = check_for_enzymes_iso(args)
        args.motif_dt = iupac_motifs_iso([args.iso])
    else:
        args.motif_dt = {}
        args.m1, args.m2 = check_for_enzymes(args)
        args.motif_dt1 = iupac_motifs(args.m1)
        args.motif_dt.update(args.motif_dt1)
        args.motif_dt2 = iupac_motifs(args.m2)
        args.motif_dt.update(args.motif_dt2)

    args.sd = int(round(0.08*args.mean, 0)) # using Sage Science CV of 8%
    args.sd = max(args.sd, int(round((args.up_bound - args.mean)/2, 0)))
    args.max = args.max if args.max else (args.mean + (6*args.sd))

    if args.a1:
        args.a1 = get_adapters(args.a1)
        if not args.a1s:
            args.a1s = len(args.a1[0])

    if args.a2:
        args.a2 = get_adapters(args.a2)
        if not args.a2s:
            args.a2s = len(args.a2[0])

    if args.q1 or args.q2:
        if not args.q1 or not args.q2:
            sys.exit('arguments q1 and q2 required')

    '''
    1.
    digest genomes one by one, producing raw digest files
    '''
    genomes_df = check_genomes(args.genome)

    print('\n1. simulating enzyme digests\n')
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
    if args.iso:
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

    read_no = "{:,}".format(int(total_freqs['counts'].sum()))
    print(f'\n{read_no} reads simulated\n')

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
    print('\nsimulating fastq formatted sequence reads')
    write_genomes(comb_file, fragment_comps, adjustment)

    #TODO consider chromosome by chromosome approach
        # this would only temporarily store sequences in the raw digest
        # raw digests (per chromosome) will be digested and stored as copies
