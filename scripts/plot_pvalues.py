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

#: histogram using -10log(pvalue), with INF ceiled at P_CEIL
P_CEIL = 1000
df1['logp'] = -np.log(df1['PValue'])
num_zero_pvalue = 0
for i, v in enumerate(df1['logp']):
    if v > P_CEIL:
        df1['logp'].iloc[i] = P_CEIL
        num_zero_pvalue += 1
x = df1['logp']
sns.distplot(x, kde=False, axlabel='-log(pvalue)', bins=100)
plt.yscale('log')
plt.ylabel('counts')
plt.show()
plt.clf()

#: histogram using raw p values
sns.distplot(df1, kde=False, axlabel='pvalue')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('counts')
plt.show()
plt.clf()

#: cumulative histogram
sns.distplot(df1, kde=False, axlabel='pvalue', hist_kws=dict(cumulative=True))
plt.yscale('log')
plt.ylabel('counts')
plt.show()
plt.clf()

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