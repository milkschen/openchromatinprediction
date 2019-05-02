import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import argparse

# input
parser = argparse.ArgumentParser()
parser.add_argument('Input_File1' , help="The .csv file containing the p-values (e.g. GWAS_pvalues_ENCFF342EGB_ENCFF853TRI.csv)")
parser.add_argument('Input_File2' , help="The .csv file containing the p-values (e.g. GWAS_pvalues_K562_rep1_ENCFF285HUZ.csv)")
parser.add_argument('ID_File' , help=".txt file containing the PWM IDs (e.g. human_pwm_ids_sorted.txt)")
args = vars(parser.parse_args())
fname1 = args['Input_File1']
fname2 = args['Input_File2']
id_name = args['ID_File']

# read in the PWM IDs
df_id = pd.read_csv(id_name, sep=' ', header=None, index_col=0)
df_id.index = range(0, df_id.shape[0])
# get the relevant columns - MA and Gene Name
df_id.columns = ['MA', 'gene_name']


df1 = pd.read_csv(fname1, sep='\t', index_col=0)
df2 = pd.read_csv(fname2, sep='\t', index_col=0)

list_sig_first_1 = []
list_sig_second_1 = []
for i, p in enumerate(df1['PValue']):
    if p == 0:
        if i < 537:
            list_sig_first_1.append(i)
        else:
            list_sig_second_1.append(i - 537)
set_sig_first_1 = set(list_sig_first_1)
set_sig_second_1 = set(list_sig_second_1)
set_sig_union_1 = set_sig_first_1.union(set_sig_second_1)
set_sig_inter_1 = set_sig_first_1.intersection(set_sig_second_1)
print ('Num of very-significant features for {0} :'.format(fname1))
print ('  counts only: {0}; weighted counts: {1}'.format(len(set_sig_first_1), len(set_sig_second_1)))
print ('  union: {0}; intersection: {1}'.format(len(set_sig_union_1), len(set_sig_inter_1)))

list_sig_first_2 = []
list_sig_second_2 = []
for i, p in enumerate(df2['PValue']):
    if p == 0:
        if i < 537:
            list_sig_first_2.append(i)
        else:
            list_sig_second_2.append(i - 537)
set_sig_first_2 = set(list_sig_first_2)
set_sig_second_2 = set(list_sig_second_2)
set_sig_union_2 = set_sig_first_2.union(set_sig_second_2)
set_sig_inter_2 = set_sig_first_2.intersection(set_sig_second_2)
print ('Num of very-significant features for {0} :'.format(fname2))
print ('  counts only: {0}; weighted counts: {1}'.format(len(set_sig_first_2), len(set_sig_second_2)))
print ('  union: {0}; intersection: {1}'.format(len(set_sig_union_2), len(set_sig_inter_2)))


set_inter_all = set_sig_inter_1.intersection(set_sig_inter_2)