import pandas as pd
import sys


df_part = pd.read_csv(sys.argv[1])
#df_total = pd.read_csv(sys.argv[2])
df_part['percent'] = sys.argv[3]
print(df_part)

#df_new = df_part.groupby(df_part['seq']).aggregate({'copies': 'sum'})
df_new = df_part.groupby(['seq','start','end','length'], as_index=False)['copies'].sum()

print(df_new)
