#!/usr/bin/env python

import numpy as np
import os
import pandas as pd


def main(digest_file, args):

    df = pd.read_csv(digest_file)
    df = df.rename(columns={'internal': 'overlaps'})
    df = find_overlaps(df)
    df = calculate_prob(df)

    prob_file = os.path.join(args.o, 'counts_' +
                             os.path.basename(args.genome) + '.csv')
    df['adj_prob'] = df['probability'] * args.comp
    df = df.reset_index(drop=True)
    df.to_csv(prob_file, index=None)
    len_freqs = get_len_freqs(df, args.max)

    return prob_file, len_freqs


def find_overlaps(df):
    '''
    manually iterate all fragments and store a list of overlapping
    adjacent fragments based on shared start and end sites, if present
    '''
    np_df = np.array(df)
    overlaps = []

    for i in np_df:
        tmp_df = df[(df['start'] > i[1]) & (df['start'] < i[2]) |
                    (df['start'] < i[1]) & (df['end'] > i[1])]
        if tmp_df.shape[0] > 0:
            overlaps.append(list(tmp_df.index.values))
        else:
            overlaps.append(None)
    df['overlaps'] = overlaps

    return df


def calculate_prob(df):
    '''
    isolate linked fragments, then find every fragment combination event
    in that cluster of linked fragments

    fragments are linked when any number overlap with another
    clusters of linked fragments end when the maximum index of a set of
    linked fragments is less than in the previous links
    '''
    np_df = np.array(df)
    prev_max = 0
    linked = []
    all_probs = []

    for i in np_df:
        if i[5]:
            linked.append(i[5])
            curr_max = max(i[5])
            if prev_max > curr_max:
                first, last, events = count_events(linked)
                all_probs.extend(apply_probs(first, last, events))
                linked = []
            prev_max = curr_max
        else:
            all_probs.append(1)

    df['probability'] = all_probs

    return df


def count_events(linked):
    '''
    using the network of linked fragments that overlap, determine every
    possible event, or, non exclusive fragment combinations
    '''
    linked_dt = {}
    events = []
    first = linked[0][0] - 1
    last = first + len(linked) - 1

    for idx1 in range(len(linked)):
        current = first + idx1
        for idx2, link in enumerate(linked):
            if idx1 != idx2 and current not in link:
                if current in linked_dt:
                    linked_dt[current].append(first + idx2)
                else:
                    linked_dt[current] = [first + idx2]

    for k in range(first, last+1):
        event = [k]
        try:
            v = linked_dt[k]
        except KeyError:
            events.append(event.copy())
            continue

        for i in v:
            event.append(i)
            for j in linked_dt[i]:
                if j in v:
                    event.append(j)
            events.append(event.copy())
            event = [k]

    sorted_events = [sorted(i) for i in events]
    events = list(set(tuple(sorted(i)) for i in sorted_events))

    return first, last, events


def apply_probs(first, last, events):
    '''
    per-fragment probability is the number of events containing that fragment
    divided by the possible events
    '''
    denom = len(events)
    cluster_probs = []

    for i in range(first, last+1):
        count = 0
        for j in events:
            count += j.count(i)
        cluster_probs.append(count/denom)

    return cluster_probs


def get_len_freqs(df, max_len):
    len_freqs = {}

    for i in range(0, max_len+1):
        len_freqs[i] = df[df['length'] == i]['adj_prob'].sum()

    return len_freqs


if __name__ == '__main__':
    print('implement standalone approach')
