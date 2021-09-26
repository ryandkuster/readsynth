import numpy as np
import pandas as pd
import random
import sys


def incomplete_digest(df, site_ls, frag_len, cut_prob, min_dist, copy_no):
    cut_pairs = generate_cuts(site_ls, cut_prob)
    #TODO add ability to modify cuts if too close (min_dist)

    indices = []

    for start, end in cut_pairs:
        hits = (df.index[(df['start'] == start) & (df['end'] == end)].tolist())
        if len(hits) == 2:
            hits = np.random.choice(hits, 1)
        if hits:
            indices.extend(hits)

    return indices


def generate_cuts(site_ls, cut_prob):
    ones = round(cut_prob * len(site_ls))
    zeroes = len(site_ls) - ones
    cut_ls = [0] * zeroes + [1] * ones
    random.shuffle(cut_ls)
    cut_ls = [i for (i, v) in zip(site_ls, cut_ls) if v]

    cut_pairs = []

    for idx, i in enumerate(cut_ls[:-1]):
        cut_pairs.append(cut_ls[idx:idx+2])

    return cut_pairs


if __name__ == '__main__':
    digest_file = sys.argv[1]
    frag_len = int(sys.argv[2])
    cut_prob = float(sys.argv[3])
    min_dist = 6
    df = pd.read_csv(digest_file)
    new_df = incomplete_digest(df, frag_len, cut_prob, min_dist, 1)
    print(len(indices))
