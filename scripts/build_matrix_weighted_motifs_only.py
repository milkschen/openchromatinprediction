#! /usr/bin/env python
import pandas as pd 
import numpy as np 
from os import path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('motif_matrix' , help="<cell_all_matrix.tsv> number of significant TF motif matches at each site")
parser.add_argument('TF_expression' , help="<exp.csv> expression level of each TF under different conditions")
parser.add_argument('out_directory', help="output directory")

args = vars(parser.parse_args())

matrix_fname = args['motif_matrix']
exp_fname = args['TF_expression']
out_dir = args['out_directory']
dnase_prefix = path.basename(matrix_fname).rstrip('_all_matrix.tsv')

matrix = np.loadtxt(matrix_fname, delimiter='\t')
exp_table = pd.read_csv(exp_fname, sep='\t', index_col=0)

n_sites = matrix.shape[0]
n_TFs = matrix.shape[1]
exp_matrix = np.zeros((n_sites, n_TFs))

# Matrix contains motif counts weighted by expression level of the corresponding TF

#rnaseq_names = ["ENCFF297CNO_TPM", "ENCFF285HUZ_TPM"]
rnaseq_names = ["ENCFF853TRI_TPM", "ENCFF305QBE_TPM"]

for rnaseq in rnaseq_names:
	
	exp_list = exp_table[rnaseq]
	
	for i in range(n_TFs):
		colID = i
		exp_matrix[:, colID] = exp_list[i] * matrix[:, i]
	
	out_fname = '%s/%s_%s_matrix_weighted_counts_only.tsv' %(out_dir, dnase_prefix, rnaseq)
	np.savetxt(out_fname, exp_matrix, delimiter='\t')
	# break

