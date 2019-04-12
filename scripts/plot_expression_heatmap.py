import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr

df_expr = pd.read_csv('../expression_levels.csv', sep='\t')
h = df_expr[['ENCFF297CNO_TPM',  'ENCFF879WBJ_TPM',  'ENCFF285HUZ_TPM', 'ENCFF853TRI_TPM', 'ENCFF305QBE_TPM']]
h.index = df_expr['# MA']
h = h.T

sns.heatmap(h, cmap='YlGnBu')
plt.show()

hh = h.sort_values(by='ENCFF297CNO_TPM', axis=1, ascending=False)
plt.clf()
sns.heatmap(hh, cmap='YlGnBu')
plt.show()

print (pearsonr(h.iloc[0], h.iloc[1]))
print (pearsonr(h.iloc[0], h.iloc[2]))
print (pearsonr(h.iloc[1], h.iloc[2]))
