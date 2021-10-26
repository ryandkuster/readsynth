import numpy as np
import os
import pandas as pd
import random
import sys

from functools import partial
from multiprocessing import Pool


def main(proj, digest_file, args):
    df = pd.read_csv(digest_file)

    # get set of all unique restriction site positions
    site_ls = list(set(df['start'].tolist() + df['end'].tolist()))
    site_ls.sort()

    print('defining fragments')
    master_dt = define_fragments(site_ls, args.max, args.cut_prob, args.n)

    print('adding copies')
    df = add_copies(df, master_dt)

    dup_file = os.path.join(proj, 'copies_' + os.path.basename(args.genome) + '.csv')
    df.drop(df[df['copies'] == 0].index, inplace = True)
    df = df.reset_index(drop=True)
    df.to_csv(dup_file, index=None)

    return dup_file


def define_fragments(site_ls, max, cut_prob, n):
    master_dt = {}

    for idx, site in enumerate(site_ls):
        fragments = [site]
        for sub_site in site_ls[idx+1:]:
            if sub_site - site <= max:
                fragments.append(sub_site)
            else:
                break

        if len(fragments) > 1:
            master_dt.update(assign_counts(fragments, cut_prob, n))

    return master_dt


def assign_counts(fragments, cut_prob, n):
    fragment_dt = {}
    start = fragments[0]

    for idx, end in enumerate(fragments[1:]):
        frequency = (cut_prob**2) * ((1-cut_prob)**idx)
        fragment_dt[(start, end)] = round(frequency * n)

    return fragment_dt


def add_copies(df, master_dt):
    '''
    using the cumulative fragment counts for each start/end site key
    in master_dt, create 'copies' column in df and add counts to each
    indexed start/end position

    in the case that m1 and m2 orientation can go on either end, randomly
    draw a forward or reverse fragment
    '''
    copy_dt = {}
    df['copies'] = 0

    for pos, count in master_dt.items():
        hits = (df.index[(df['start'] == pos[0]) & (df['end'] == pos[1])].tolist())
        if hits:
            if len(hits) == 2:
                hits = [random.choice(hits) for i in range(count)]
            else:
                hits = [hits[0] for i in range(count)]

            for i in hits:
                if i in copy_dt:
                    copy_dt[i] += 1
                else:
                    copy_dt[i] = 1

    for idx, count in copy_dt.items():
        df.loc[idx, 'copies'] = count

    return df


if __name__ == '__main__':
    print('implement standalone approach')
