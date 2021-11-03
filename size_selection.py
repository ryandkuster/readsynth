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
    len_dt = length_dict(df, args.mean, args.sd)
    a = [len_dt[i] for i in range(args.mean, (args.up_bound+(args.up_bound-args.mean))+1)]
    avg_upper = sum(a)/len(a)
    scale_by = avg_upper/gauss_pdf(args.mean, args.sd, args.up_bound)
    draw_dt = get_draw_dict(args.mean, args.sd, len_dt, scale_by)
    sampled_df = draw_reads(df, col_names, draw_dt)

    index_names = sampled_df[sampled_df['counts'] == 0].index
    sampled_df.drop(index_names, inplace=True)
    samp_no = sampled_df['counts'].sum()
    print(f'fragments sampled around mean of {args.mean}bp : {samp_no}')
    sampled_df.reset_index(inplace=True, drop=True)
    sampled_file = os.path.join(proj, 'sampled_' + os.path.basename(args.genome) + '.csv')
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


def gauss_pdf(mean, sd, x):
    '''
    return the probability density of x from a normal distribution
    note: the sum of pdf for all points x = 1
    '''
    pdf = (1/sd*np.sqrt(2*np.pi))*np.exp((-1/2)*((x-mean)/sd)**2)
    return pdf


def length_dict(df, mean, sd):
    '''
    create a len_dt, storing the count for each length
    '''
    len_dt = {}
    for i in range(0, mean + sd*6 + 1):
        len_dt[i] = df[df.full_length == i]['copies'].sum()

    return len_dt


def get_draw_dict(mean, sd, len_dt, scale_by):
    '''
    create a dictionary of draw numbers
    '''
    draw_dt = {}

    for x in range(0, mean + 6*sd + 1):
        draw_counts = int(round(gauss_pdf(mean, sd, x)*scale_by, 0))
        data_counts = len_dt[x]
        draw_dt[x] = int(round(min(draw_counts, data_counts)*0.65, 0)) #TODO average recovery

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
    sampled_df['counts'] = sampled_df['counts'].astype(int)
    sampled_df['full_length'] = sampled_df['full_length'].astype(int)

    return sampled_df


if __name__ == '__main__': #TODO add args parse
    dup_file = sys.argv[1]
    proj = sys.argv[2]
    main(dup_file, proj, args)
