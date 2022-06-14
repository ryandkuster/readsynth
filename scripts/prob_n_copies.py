#!/usr/bin/env python

import os
import pandas as pd


def main(digest_file, args):

    df = pd.read_csv(digest_file)
    internal_max = df['internal'].max()
    copies_dt = copies_dict(internal_max, args)

    df = apply_approach(df, copies_dt)
    prob_file = os.path.join(args.o, 'counts_' +
                             os.path.basename(args.g) + '.csv')
    df.drop(df[df['probability'] == 0].index, inplace=True)

    if df.duplicated(subset=['start','end']).any():
        df = bidirectional_weights(df)

    df['adj_prob'] = df['probability'] * args.comp
    df = df.reset_index(drop=True)
    df.to_csv(prob_file, index=None)
    len_freqs = get_len_freqs(df, args.max)

    return prob_file, len_freqs


def copies_dict(internal_max, args):
    '''
    return a dictionary of expected probabilities for a fragment
    given the fragment contains i internal cut sites
    '''
    copies_dt = {}

    for i in range(internal_max+1):
        copies_dt[i] = (args.c**2) * ((1-args.c)**i)

    return copies_dt


def apply_approach(df, copies_dt):
    df['probability'] = df.apply(
        lambda row: get_copies(copies_dt, internal=row['internal']), axis=1)

    return df


def get_copies(copies_dt, internal):
    return copies_dt[internal]


def bidirectional_weights(df):
    df.loc[(df.duplicated(subset=['start','end'], keep=False)),
        'probability'] = df['probability']/2

    return df

def get_len_freqs(df, max_len):
    len_freqs = {}

    for i in range(0, max_len+1):
        len_freqs[i] = df[df['length'] == i]['adj_prob'].sum()

    return len_freqs


if __name__ == '__main__':
    print('implement standalone approach')
