#!/usr/bin/env python

import numpy as np
import os
import pandas as pd
import random
import sys


def main(dup_file, proj, args):
    """
    sampling is performed by first assessing the count of fragments for
    each possible read length and then generating a normal distribution
    that includes enough counts in the mean+2sd range to include a 1X
    coverage of the genome
    """
    df = pd.read_csv(dup_file)
    col_names = [col for col in df.columns]

    if len(df) == 0:
        sys.exit('no fragments produced with current settings')

    df = modify_length(df, args)
    len_dt, total_reads = length_dict(df, args)
    print('len dict')
    draw_ls = normal_distribution(args, len_dt, total_reads)
    print('norm dist')
    draw_dt = draw_dict(df, draw_ls)
    print('draw dict')
    sampled_df = draw_reads(df, col_names, draw_dt)

    #TODO the sampled file is now complete
    print('reads drawn')
    index_names = sampled_df[sampled_df['counts'] == 0].index
    sampled_df.drop(index_names, inplace=True)
    samp_no = sampled_df['counts'].sum()
    print(f'fragments sampled around mean of {args.mean}bp : {samp_no}')
    sampled_df.reset_index(inplace=True, drop=True)
    sampled_file = os.path.join(proj, 'sampled_' + os.path.basename(args.genome)     + '.csv')
    sampled_df.to_csv(sampled_file, index=None)

    return sampled_file


def modify_length(df, args):
    '''
    if adapters are present, calculate average length for each
    sum averages to create a modifier for the fragment length
    after adapter ligation

    sort and return df with full_length column
    '''
    if args.a1 and args.a2:
        a1_avg = sum([len(i) for i in args.a1.keys()]) / len(args.a1.keys())
        a2_avg = sum([len(i) for i in args.a2.keys()]) / len(args.a2.keys())
        modifier = round(a1_avg + a2_avg)
    elif args.a1:
        a1_avg = sum([len(i) for i in args.a1.keys()]) / len(args.a1.keys())
        modifier = round(2 * a1_avg)
    else:
        modifier = 0

    df['full_length'] = df['length'] + modifier
    df.sort_values(['full_length'], ascending=[True], inplace=True)
    df.reset_index(inplace=True, drop=True)

    return df


def length_dict(df, args):
    '''
    create a len_dt, storing the count for each length
    '''
    len_dt = {}
    for i in range(args.mean, args.mean + 2*args.sd):
        len_dt[i] = df[df.full_length == i]['copies'].sum()
    total_reads = sum(len_dt.values()) * 2

    return len_dt, total_reads


def normal_distribution(args, len_dt, total_reads):
    '''
    produce a normal distribution that includes mean + 2sd counts
    '''
    keep_going = True

    while keep_going is True:
        keep_going = False
        print(total_reads)
        draw_ls = np.random.normal(loc=args.mean,scale=args.sd,size=total_reads)
        draw_ls = [round(i) for i in draw_ls]
        for i, len_count in len_dt.items():
            if len_count > draw_ls.count(i):
                total_reads = round(total_reads*1.1)
                keep_going = True
                break

    return draw_ls


def draw_dict(df, draw_ls):
    '''
    create a dictionary of draw numbers
    '''
    draw_dt = {}

    for i in range(min(draw_ls), max(draw_ls)+1):
        draw_counts = draw_ls.count(i)
        data_counts = df[df.full_length == i]['copies'].sum()
        draw_dt[i] = min(draw_counts, data_counts)

    return draw_dt


def draw_reads(df, col_names, draw_dt):
    '''
    for each fragment length, randomly draw reads
    '''
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

    return sampled_df


if __name__ == '__main__': #TODO add args parse
    dup_file = sys.argv[1]
    proj = sys.argv[2]
    main(dup_file, proj, args)
