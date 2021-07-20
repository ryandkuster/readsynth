import numpy as np
import pandas as pd
import sys


def incomplete_digest(df, frag_len, cut_prob, min_dist, copy_no):
    site_ls = list(set(df['start'].tolist() + df['end'].tolist()))
    site_ls.sort()
    cut_ls = probability_digest(site_ls, cut_prob, min_dist)
    site_ls = [i for (i, v) in zip(site_ls, cut_ls) if v]
    start = site_ls[0]

    indices = []

    for site in site_ls[1:]:
        if site - frag_len > start:
            start = site
        else:
            end = site
            hits = (df.index[(df['start'] == start) & (df['end'] == end)].    tolist())
            if len(hits) == 2:
                hits = np.random.choice(hits, 1)
            if hits:
                indices.extend(hits)
            start = site

    return indices


def probability_digest(site_ls, cut_prob, min_dist):
    last_cut = -min_dist
    cut_ls = []
    for idx, site in enumerate(site_ls):
        if site - last_cut < min_dist:
            cut = 0 #TODO adjust for variable enzyme behavior
            cut_ls.append(cut)
        else:
            cut = np.random.choice([1, 0], 1, p=[cut_prob, 1-cut_prob])
            cut_ls.extend(cut)
        if cut:
            last_cut = site

    return cut_ls


if __name__ == '__main__':
    digest_file = sys.argv[1]
    frag_len = int(sys.argv[2])
    cut_prob = float(sys.argv[3])
    min_dist = 6
    df = pd.read_csv(digest_file)
    new_df = incomplete_digest(df, frag_len, cut_prob, min_dist, 1)
    print(len(indices))
