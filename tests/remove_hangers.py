import pandas as pd
import sys

motif_dt = {'CATG': 0, 'ATGCAT': 5}

df = pd.read_csv(sys.argv[1])

for mot, front in motif_dt.items():
    back = len(mot) - front
    df.loc[(df['m1'] == mot) & (df['reverse'] == 0), 'seq'] = df['seq'].str[front:]
    df.loc[(df['m1'] == mot) & (df['reverse'] == 1), 'seq'] = df['seq'].str[:-back]
    df.loc[(df['m2'] == mot) & (df['reverse'] == 0), 'seq'] = df['seq'].str[:-back]
    df.loc[(df['m2'] == mot) & (df['reverse'] == 1), 'seq'] = df['seq'].str[front:]
print(df)
