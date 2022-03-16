#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd

'''
inputs:
1 - sampled_df.csv
2 - simulated_R1.fastq
3 - simulated_R2.fastq
4 - a1s (adapter 1 pos where sequencing occurs)
5 - a2s (adapter 2 pos where sequencing occurs)
6 - length of fastq reads
'''

sampled_df = pd.read_csv(sys.argv[1], index_col=0)
print(sampled_df)
read_len = [0 for i in range(int(sys.argv[6]))]
scores = list('!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJK')
scores_dt1 = {i: (read_len.copy(), read_len.copy()) for i in scores}
scores_dt2 = {i: (read_len.copy(), read_len.copy()) for i in scores}


def open_files(r1, r2):
    a1s = int(sys.argv[4])
    a2s = int(sys.argv[5])
    i = 0
    entry1 = []
    entry2 = []
    for line1, line2 in zip(r1, r2):
        i += 1
        entry1.append(line1.rstrip())
        entry2.append(line2.rstrip())
        if i == 4:
            i = 0
            check_errors(entry1, a1s)
            check_errors(entry2, a2s)
            entry1 = []
            entry2 = []


def check_errors(entry, start):
    idx, read = get_header_parts(entry[0])
    mut_seq = entry[1]
    if read == '1':
        seq = sampled_df.iloc[idx,0][start:]
        seq = seq.ljust(len(mut_seq), 'G')
        count_miscalls(seq, mut_seq, entry[3], scores_dt1)
    else:
        seq = reverse_comp(sampled_df.iloc[idx,0])[start:]
        seq = seq.ljust(len(mut_seq), 'G')
        count_miscalls(seq, mut_seq, entry[3], scores_dt2)


def get_header_parts(header):
    idx = int(header.split(':')[1])
    read = header[-1]
    return idx, read


def reverse_comp(seq):
    revc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    new = ''
    for base in reversed(seq):
        new += revc[base]
    return new


def count_miscalls(seq, mut_seq, scores, scores_dt):
    for idx, (i, j, k) in enumerate(zip(seq, mut_seq, scores)):
        if i != j:
            scores_dt[k][0][idx] = scores_dt[k][0][idx] + 1
        scores_dt[k][1][idx] = scores_dt[k][1][idx] + 1


def write_to_csv(scores_dt, outname):
    output_df = pd.DataFrame()
    output_df['base'] = [i+1 for i in range(int(sys.argv[6]))]
    for idx, (k, v) in enumerate(scores_dt.items()):
        expected = 10**(-idx/10)
        output_df[str(k + '_miss')] = v[0]
        output_df[str(k + '_total')] = v[1]
        if sum(output_df[str(k + '_total')]) > 0:
            output_df[str(k + 'O_rate')] = [i/j if j > 0 else 0 for i, j in zip(v[0],v[1])]
            output_df[str(k + 'E_rate')] = [expected if i > 0 else 0 for i in v[1]]
            output_df[str(k + '_O-E')] = [i-j for i, j in zip(output_df[str(k + 'O_rate')], output_df[str(k + 'E_rate')])]

    output_df = output_df.loc[:, (output_df != 0).any(axis=0)]
    output_df.to_csv(outname, index=None)


with open(sys.argv[2]) as r1, open(sys.argv[3]) as r2:
    open_files(r1, r2)

write_to_csv(scores_dt1, 'R1_simulation_summary.csv')
write_to_csv(scores_dt2, 'R2_simulation_summary.csv')
