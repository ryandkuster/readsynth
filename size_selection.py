#!/usr/bin/env python

import numpy as np
import os
import pandas as pd
import random
import sys


def main(df, args):
    """
    sampling is performed by first assessing the count of fragments for
    each possible read length and then generating a normal distribution
    that includes enough counts in the mean+2sd range to include a 1X
    coverage of the genome

    fragment_comps:
    the composition, or relative rate of occurence, per fragment length,
    across the entirety of input genome digests

    adjustment:
    based on the sum of all possible fragment occurences, adjusts the
    value of -n (total read count) to accound for fragment totals that
    are > or < than n (e.g. for 1,000,000 total reads, 100,000 size-selected
    fragments will adjust the weight of each read to 10, while 10,000,000
    fragments will need to be adjusted to 0.1)
    """

    if len(df) == 0:
        sys.exit('no fragments produced with current settings')

    modify_length(df, args)

    len_dt = length_dict(df, args.mean, args.sd)

    a = [len_dt[i] for i in range(max(0, args.mean), (args.up_bound+(args.up_bound-args.mean))+1)]
    avg_upper = sum(a)/len(a)
    scale_by = avg_upper/gauss_pdf(args.mean, args.sd, args.up_bound)
    draw_dt = get_draw_dict(args.mean, args.sd, len_dt, scale_by)
    adjustment = args.n / sum(draw_dt.values()) #TODO
    fragment_comps = {}

    for i in range(0, (args.up_bound+(args.up_bound-args.mean))+1):
        fragment_comps[i] = draw_dt[i] / len_dt[i]

    return fragment_comps, adjustment


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

    args.mean -= modifier
    args.up_bound -= modifier


def gauss_pdf(mean, sd, x):
    '''
    return the probability density of x from a normal distribution
    note: the sum of pdf for all points x = 1
    '''
    pdf = (1/sd*np.sqrt(2*np.pi))*np.exp((-1/2)*((x-mean)/sd)**2)
    return pdf


def length_dict(df, mean, sd):
    '''
    create a len_dt, storing the combined probability for each length
    '''
    len_dt = {}
    for i in range(0, mean + sd*6 + 1):
        len_dt[i] = df[df.length == i]['sum_prob'].sum()

    return len_dt


def get_draw_dict(mean, sd, len_dt, scale_by):
    '''
    create a dictionary of draw numbers
    '''
    draw_dt = {}

    for x in range(0, mean + 6*sd + 1):
        draw_counts = gauss_pdf(mean, sd, x)*scale_by
        data_counts = len_dt[x]
        draw_dt[x] = min(draw_counts, data_counts)*0.65 #TODO average recovery

    return draw_dt


if __name__ == '__main__': #TODO add args parse
    pass
