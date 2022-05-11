import numpy as np
import os
import pandas as pd
import random
import sys

from functools import partial
from multiprocessing import Pool


def main(digest_file, args):

    df = pd.read_csv(digest_file)
    df = df.rename(columns={'internal': 'overlaps'})
    df = find_overlaps(df)
    calculate_prob(df)
    #TODO add reduced weight for overlapping fragments

    prob_file = os.path.join(args.o, 'counts_' + os.path.basename(args.genome) + '.csv')
    #df['adj_prob'] = df['probability'] * args.comp
    df = df.reset_index(drop=True)
    df.to_csv(prob_file, index=None)

    sys.exit() #TODO
    return prob_file


def find_overlaps(df):
    np_df = np.array(df)
    overlaps = []

    for i in np_df:
        tmp_df = df[(df['start'] > i[1]) & (df['start'] < i[2]) | \
                    (df['start'] < i[1]) & (df['end'] > i[1])]
        if tmp_df.shape[0] > 0:
            overlaps.append(list(tmp_df.index.values))
        else:
            overlaps.append(None)
    df['overlaps'] = overlaps

    return df


def calculate_prob(df):
    '''
    calculate the probability of each fragment by first isolating linked
    fragments, then finding every fragment combination event
    per-fragment probability is the number of events containing that fragment
    divided by the possible events
    '''
    np_df = np.array(df)
    prev_max = 0
    linked = []

    for i in np_df:
        if i[5]:
            print(i[5])
            linked.append(i[5])
            curr_max = max(i[5])
            if prev_max > curr_max:
                print(linked)
                process_links(linked)
                print('\n')
                linked = []

            prev_max = curr_max


def process_links(linked):
    '''
    using the network of linked fragments that overlap, determine every
    possible outcome
    '''
    first_idx = linked[0][0] - 1
    print(first_idx)

    #TODO LEFT OFF HERE?


if __name__ == '__main__':
    print('implement standalone approach')
