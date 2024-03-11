#!/usr/bin/env python

import argparse
import gzip
import numpy as np
import os

from gzip_test import test_unicode



def parse_user_input():
    parser = argparse.ArgumentParser(description='sample q scores to length')

    parser.add_argument('-o', type=str, required=True,
                        help='path to store output')

    parser.add_argument('-r1', type=str, required=True,
                        help='R1 fastq file to sample Q scores')

    parser.add_argument('-r2', type=str, required=True,
                        help='R2 fastq file to sample Q scores')

    parser.add_argument('-p', type=int, required=True,
                        help='if using r1/r2 for profile, percent of reads to sample')

    parser.add_argument('-l', type=int, required=False,
                        help='desired read length of final simulated reads')

    args = parser.parse_args()

    return args


def open_fastq(r1_out, r2_out, args):

    print(f'sampling {args.p} percent of scores from {args.r1} and {args.r2}')
    perc_keep = [int(args.p)/100, 1-(int(args.p)/100)]
    i = 0
    if args.l:
        least = args.l + 1

    if test_unicode(args.r1):
        opened_r1 = gzip.open(args.r1, 'rt')
        opened_r2 = gzip.open(args.r2, 'rt')
    else:
        opened_r1 = open(args.r1)
        opened_r2 = open(args.r2)

    r1_out = open(r1_out, 'w')
    r2_out = open(r2_out, 'w')


    for line1, line2 in zip(opened_r1, opened_r2):
        i += 1
        if i == 4:
            i = 0
            keep = np.random.choice([True, False], 1, p=perc_keep)
            if not keep:
                pass
            elif args.l and len(line1) >= least and len(line2) >= least:
                r1_out.write(line1.rstrip()[:args.l] + '\n')
                r2_out.write(line2.rstrip()[:args.l] + '\n')
            else:
                r1_out.write(line1)
                r2_out.write(line2)

    opened_r1.close()
    opened_r2.close()
    r1_out.close()
    r2_out.close()


if __name__ == '__main__':
    args = parse_user_input()
    if args.r1 or args.r2:
        if not args.p:
            sys.exit('please provide input \"-p\" for percent of fq to sample')
        r1_out = os.path.join(args.o, os.path.basename(args.r1)+'_sampled_scores.csv')
        r2_out = os.path.join(args.o, os.path.basename(args.r2)+'_sampled_scores.csv')
        open_fastq(r1_out, r2_out, args)

