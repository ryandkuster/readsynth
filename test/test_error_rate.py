#
import math
import sys

'''
sys.argv[1] is the original sim file, no errors
sys.argv[2] is the error-induced reads
'''

def count_errors(line1, line2, error_dt):
    for pos, (b1, b2) in enumerate(zip(line1, line2)):
        if b1 != b2:
            if pos in error_dt:
                error_dt[pos] += 1
            else:
                error_dt[pos] = 1

    return error_dt


def evaluate_rate(error_dt, read_no):
    print(read_no)
    print(error_dt)

    for i in range(max(error_dt)):
        try:
            prob = error_dt[i]/read_no
        except KeyError:
            prob = 1/10000
        q = -10 * math.log10(prob)
        print(f'{i} : {q}')


if __name__ == "__main__":
    error_dt = {}
    line_no = 0
    read_no = 0

    with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2:
        for line1, line2 in zip(f1, f2):
            line_no += 1
            if line_no == 2:
                error_dt = count_errors(line1, line2, error_dt)
            if line_no == 4:
                read_no += 1
                line_no = 0

    evaluate_rate(error_dt, read_no)
