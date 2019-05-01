import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

### This file plots the expression levels ###


df = pd.read_csv('expression_levels.csv', sep='\t')

sns.lineplot(x=df.index, y='ENCFF853TRI_TPM', data=df, label='GM12878')
sns.lineplot(x=df.index, y='ENCFF297CNO_TPM', data=df, label='K562')
plt.ylabel('TPM')
plt.xlabel('motif ID')
plt.show()
plt.clf()