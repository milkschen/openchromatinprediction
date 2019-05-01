import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

### This file plots the correlation heatmap between features ###

parser = argparse.ArgumentParser()
parser.add_argument('Input_File' , help="The .tsv file containing the TPM matrices (e.g. ENCFF342EGB_ENCFF853TRI_TPM_matrix.tsv)")
args = vars(parser.parse_args())
fname = args['Input_File']

# read in the files
# df1 = pd.read_csv('ENCFF342EGB_ENCFF305QBE_TPM_matrix.tsv', header=None, sep='\t')
df1 = pd.read_csv(fname, header=None, sep='\t')

# get the counts and weighted counts
df_count = df1.iloc[:,:537]
df_weightedcount = df1.iloc[:,537:]

# remove the zero columns
zero_cols_count = [ col for col, is_zero in ((df_count == 0).sum() == df_count.shape[0]).items() if is_zero ]
df_count.drop(zero_cols_count, axis=1, inplace=True)
zero_cols_weightedcount = [ col for col, is_zero in ((df_weightedcount == 0).sum() == df_weightedcount.shape[0]).items() if is_zero ]
df_weightedcount.drop(zero_cols_weightedcount, axis=1, inplace=True)

# final preprocessing
corr_count = df_count.corr()
corr_weightedcount = df_weightedcount.corr()
corr_weightedcount.index = corr_weightedcount.index - 537
corr_weightedcount.columns = corr_weightedcount.columns - 537

# plot
plt.figure(1)
plt.subplot(1, 2, 1)
clustergrid_count = sns.clustermap(corr_count, cmap='YlGnBu')

plt.title('Motif count only')
plt.xlabel('Motif id')
plt.ylabel('Motif id')
plt.subplot(1, 2, 2)
clustergrid_weightedcount = sns.clustermap(corr_weightedcount, cmap='YlGnBu')

plt.title('Weighted motif count')
plt.xlabel('Motif id')
plt.ylabel('Motif id')
plt.show()