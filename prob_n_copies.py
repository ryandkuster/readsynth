import numpy as np
import os
import pandas as pd
import random
import sys

from functools import partial
from multiprocessing import Pool


def main(proj, digest_file, args):

    df = pd.read_csv(digest_file)
    internal_max = df['internal'].max()

    copies_dt = copies_dict(internal_max, args.cut_prob, args.n)

    # get set of all unique restriction site positions
    df = apply_approach(df, copies_dt)

    dup_file = os.path.join(proj, 'copies_' + os.path.basename(args.genome) + '.csv')
    df.drop(df[df['probability'] == 0].index, inplace = True)
    df = df.reset_index(drop=True)
    df.to_csv(dup_file, index=None)

    return dup_file


def copies_dict(internal_max, cut_prob, n):
    '''
    return a dictionary of expected probabilities for a fragment
    given the fragment contains i internal cut sites
    '''
    copies_dt = {}

    for i in range(internal_max+1):
        copies_dt[i] = (cut_prob**2) * ((1-cut_prob)**i)

    return copies_dt


def apply_approach(df, copies_dt):
    df['probability'] = df.apply(
        lambda row: get_copies(copies_dt, internal=row['internal']), axis=1)

    return df


def get_copies(copies_dt, internal):
    return copies_dt[internal]


if __name__ == '__main__':
    print('implement standalone approach')
