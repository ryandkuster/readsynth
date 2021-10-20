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

    pool = Pool(args.t)
    pool_ls = [i+1 for i in range(args.n)]
    pool_part = partial(incomplete_digest, site_ls, args.max, args.cut_prob, args.min)
    indices = pool.map(pool_part, pool_ls)
    pool.close

    master_dt = combine_dts(indices)

    df = add_copies(df, master_dt)

    dup_file = os.path.join(proj, 'copies_' + os.path.basename(args.genome) + '.csv')
    df.drop(df[df['copies'] == 0].index, inplace = True)
    df = df.reset_index(drop=True) 
    df.to_csv(dup_file, index=None)

    return dup_file


def incomplete_digest(site_ls, frag_len, cut_prob, min_dist, copy_no):
    '''
    multiprocessed function

    for every RE cut site, simulate cut probability and return a dictionary
    where keys are tuples of fragment start end list, values are 1
    '''
    cut_ls = generate_cuts(site_ls, cut_prob)
    if min_dist and len(cut_ls) > 0:
        cut_ls = limit_cut_proximity(cut_ls, min_dist)
    cut_pairs = generate_cut_pairs(cut_ls, frag_len)
    indices = {}

    for start, end in cut_pairs:
        indices[(start, end)] = 1

    return indices


def generate_cuts(site_ls, cut_prob):
    '''
    for every RE cut site, simulate a cut (1) or no cut (0)

    return genome positions as a list with each cut site
    '''
    cut_ls = []

    for idx, site in enumerate(site_ls):
        cut = np.random.choice([1, 0], 1, p=[cut_prob, 1-cut_prob])
        cut_ls.extend(cut)

    cut_ls = [i for (i, v) in zip(site_ls, cut_ls) if v]

    return cut_ls


def generate_cut_pairs(cut_ls, frag_len):
    '''
    return a list of start and end sites for resulting fragments

    discard fragments longer than frag_len
    '''
    cut_pairs = []

    for idx, i in enumerate(cut_ls[:-1]):
        if cut_ls[idx+1] - cut_ls[idx] <= frag_len:
            cut_pairs.append(cut_ls[idx:idx+2])

    return cut_pairs


def limit_cut_proximity(cut_ls, min_dist):
    '''
    randomly remove cut sites from cut_ls if distance is less than
    min_dist and store in new list, limit_cut_ls

    compare every cut site with the last value of limit_cut_ls

    '''
    limit_cut_ls = [cut_ls[0]]

    for idx, i in enumerate(cut_ls[1:]):
        if i - limit_cut_ls[-1] < min_dist:
            limit_cut_ls[-1] = random.choice([i, limit_cut_ls[-1]])
        else:
            limit_cut_ls.append(i)

    return limit_cut_ls


def combine_dts(indices):
    master_dt = {}
    for tup_dt in indices:
        tmp_dt = {i:tup_dt[i] if i not in master_dt.keys() else tup_dt[i] + master_dt[i] for i in tup_dt} 
        master_dt.update(tmp_dt)

    return master_dt


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
    digest_file = sys.argv[1]
    frag_len = int(sys.argv[2])
    cut_prob = float(sys.argv[3])
    min_dist = 6
    df = pd.read_csv(digest_file)
    new_df = incomplete_digest(frag_len, cut_prob, min_dist, 1)
    print(len(indices))
