import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

gm_coef2_fname= "../results/ENCFF342EGB_ENCFF305QBE_TPM_matrix_weighted_counts_only_lr_ranks.tsv"
gm_coef2 = pd.read_csv(gm_coef2_fname, sep='\t', header=None, index_col='MOTIF_ID', names=["MOTIF_ID", "TF NAME", "COEF"])

k562_coef_fname= "../results/K562_rep1_ENCFF285HUZ_TPM_matrix_weighted_counts_only_lr_ranks.tsv"
k562_coef = pd.read_csv(k562_coef_fname, sep='\t', header=None, index_col='MOTIF_ID', names=["MOTIF_ID", "TF NAME", "COEF"])
k562_coef_sorted = k562_coef.sort_index()
gm_coef_sorted = gm_coef2.sort_index()

plt.figure(figsize=(10, 7.5))
plt.plot(k562_coef_sorted['COEF'].values, gm_coef_sorted['COEF'].values, '.', markersize=7, alpha=0.3)
plt.xlabel('K562', FontSize=15)
plt.ylabel('GM12878', FontSize=15)
plt.title('Weighted only coefficients', FontSize=15)
plt.savefig('../figs/weighted_only_coefficients.png', dpi=200)
# plt.show()

gm_coef3_fname = "../results/ENCFF342EGB_lr_ranks.tsv"
gm_coef3 = pd.read_csv(gm_coef3_fname, sep='\t', header=None, index_col='MOTIF_ID', names=["MOTIF_ID", "TF NAME", "COEF"])
k562_coef_fname2= "../results/K562_rep1_lr_scaled_ranks.tsv"
k562_coef2 = pd.read_csv(k562_coef_fname2, sep='\t', header=None, index_col='MOTIF_ID', names=["MOTIF_ID", "TF NAME", "COEF"])
k562_coef_sorted2 = k562_coef2.sort_index()

gm_coef_sorted3 = gm_coef3.sort_index()
plt.figure(figsize=(10, 7.5))
plt.plot(k562_coef_sorted2['COEF'].values, gm_coef_sorted3['COEF'].values, '.', markersize=7, alpha=0.3)
plt.xlabel('K562', FontSize=15)
plt.ylabel('GM12878', FontSize=15)
plt.title('Counts only coefficients', FontSize=15)
plt.savefig('../figs/counts_only_coefficients.png', dpi=200)
# plt.show()

