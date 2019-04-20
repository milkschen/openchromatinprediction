import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# df1 = pd.read_csv('GWAS_pvalues_ENCFF342EGB_ENCFF305QBE.csv', sep='\t', index_col=0)
df1 = pd.read_csv('GWAS_pvalues_ENCFF342EGB_ENCFF853TRI.csv', sep='\t', index_col=0)

min_not_inf = 0
for i, v in enumerate(df1.sort_values('PValue')['PValue']):
    if float(v) > 0:
        min_not_inf = v
        # print (i, v)
        break

for i in range(df1.shape[0]):
    if df1['PValue'][df1.index[i]] == 0:
        df1['PValue'][df1.index[i]] = min_not_inf
        # df1['-10logP'][df1.index[i]] = 0

df1['-10logP'] = -10 * np.log10(df1['PValue'])
df1['weighted'] = (df1.index < 537)
g = df1.groupby('weighted').groups

sns.scatterplot(g[0]-537, y=df1.loc[g[0]]['-10logP'], label='weighted', marker='|')
sns.scatterplot(g[1], y=df1.loc[g[1]]['-10logP'], label='count', marker='_')
plt.legend(loc='right')
plt.ylabel('-10log(pvalue)')
plt.xlabel('motif ID')
plt.show()
plt.clf()