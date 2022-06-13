#!/usr/bin/env python

import json
import math
import numpy as np
import os
import pandas as pd
import sys


def main(df, args):
    """
    sampling is performed by first assessing the count of fragments for
    each possible read length and then generating a normal distribution
    that includes enough counts in the mean+2sd range to include a 1X
    coverage of the genome

    scale_by is the density of all fragments at the point x that is the
    upper_bound defined by the user; this value is approximately the
    expected counts of fragments this length over all genomes assuming
    1X copies of each genome followed by adjustment for composition

    fragment_comps:
    the composition, or relative rate of occurence, per fragment length,
    across the entirety of input genome digests (ie for each length, this
    value should fall between 0 and 1 and determines the percent of reads
    this size to be produced from size selection)

    adjustment:
    based on the sum of all possible fragment occurences, adjusts the
    value of -n (total read count) to account for fragment totals that
    are > or < than n (e.g. for 1,000,000 total reads, 100,000 size-selected
    fragments will adjust the weight of each read to 10, while 10,000,000
    fragments will need to be adjusted to 0.1)
    """

    if len(df) == 0:
        sys.exit('no fragments produced with current settings')

    len_dt = length_dict(df, args)

    if args.dist:
        fragments_dt = get_custom_dist(len_dt, args)
        draw_dt = get_custom_draw_dict(fragments_dt, args)
    else:
        scale_by = get_scale_gauss(len_dt, args)
        draw_dt = get_draw_dict(len_dt, scale_by, args)

    if sum(draw_dt.values()) == 0:
        print(f'no fragments produced in the size selection range, quitting')
        sys.exit()

    adjustment = args.n / sum(draw_dt.values())
    fragment_comps = {}

    for i in range(0, args.max + 1):
        if len_dt[i] == 0 or draw_dt[i] == 0:
            fragment_comps[i] = 0
        else:
            fragment_comps[i] = draw_dt[i] / len_dt[i]

        if math.isnan(fragment_comps[i]):
            fragment_comps[i] = 0

    return fragment_comps, adjustment


def length_dict(df, args):
    '''
    create a len_dt, storing the combined probability for each length
    '''
    len_dt = {}
    for i in range(0, args.max + 1):
        len_dt[i] = df[df.length == i]['sum_prob'].sum()

    return len_dt


def get_custom_dist(len_dt, args):
    with open(args.dist) as f_o:
        tmp_dt = json.load(f_o)

    tmp_dt = {int(k): int(v) for k, v in tmp_dt.items()}
    fragments_dt = {}

    for i in len_dt:
        if i not in tmp_dt:
            fragments_dt[i] = 0
        else:
            fragments_dt[i] = tmp_dt[i]

    return fragments_dt


def get_custom_draw_dict(fragments_dt, args):
    draw_dt = {}

    for x in range(0, args.max + 1):
        draw_dt[x] = fragments_dt[x]

    return draw_dt


def get_scale_gauss(len_dt, args):
    x_dist = max(args.x - args.mean, args.x - 0)
    x_range = [len_dt[i] for i in range(args.x-x_dist, args.x+x_dist+1)]

    try:
        avg_x = sum(x_range)/len(x_range)
    except ZeroDivisionError:
        print(f'no fragments produced in the range of mean {args.mean}bp ' +
              f'+/- {args.sd}bp (adjusted for adapter lengths), quitting')
        sys.exit()

    scale_by = avg_x / gauss_pdf(args.x, args)

    return scale_by


def gauss_pdf(x, args):
    '''
    return the probability density of x from a normal distribution
    note: the sum of pdf for all points x = 1
    '''
    pdf = (1 / args.sd * np.sqrt(2 * np.pi)) * np.exp((-1 / 2) * ((x - args.mean) / args.sd)**2)

    return pdf


def get_draw_dict(len_dt, scale_by, args):
    '''
    create a dictionary of draw numbers in the range 0 to the maximum
    fragment length based on the Gaussian probability density function
    at a given point, x

    scale_by changes the density for a point (probabilities sum to 1) to
    the units of adj_prob from the compositional dataset

    at the intersection of the above distribution and the possible
    fragments produced by digestion, take the minimum count for each x
    '''
    draw_dt = {}

    for x in range(0, args.max + 1):
        draw_counts = gauss_pdf(x, args)*scale_by
        data_counts = len_dt[x]
        draw_dt[x] = min(draw_counts, data_counts) * 0.65 #Sage average recovery

    return draw_dt


if __name__ == '__main__': #TODO add args parse
    pass
