import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import argparse

### This file plots the expression levels ###

# input
parser = argparse.ArgumentParser()
parser.add_argument('Input_File' , help=".csv file containing expression levels (e.g. expression_levels.csv)")
args = vars(parser.parse_args())
fname = args['Input_File']


df = pd.read_csv(fname, sep='\t')

sns.lineplot(x=df.index, y='ENCFF853TRI_TPM', data=df, label='GM12878')
sns.lineplot(x=df.index, y='ENCFF297CNO_TPM', data=df, label='K562')
plt.ylabel('TPM')
plt.xlabel('motif ID')
plt.show()
plt.clf()