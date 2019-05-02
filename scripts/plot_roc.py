import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse

x1 = np.loadtxt("../results/ENCFF342EGB_ENCFF305QBE_TPM_lr_scaled_fpr_100.txt")
y1 = np.loadtxt("../results/ENCFF342EGB_ENCFF305QBE_TPM_lr_scaled_tpr_100.txt")
x2 = np.loadtxt("../results/ENCFF342EGB_ENCFF305QBE_TPM_matrix_weighted_counts_only_lr_fpr_100.txt")
y2 = np.loadtxt("../results/ENCFF342EGB_ENCFF305QBE_TPM_matrix_weighted_counts_only_lr_tpr_100.txt")
x3 = np.loadtxt("../results/ENCFF342EGB_lr_fpr_100.txt")
y3 = np.loadtxt("../results/ENCFF342EGB_lr_tpr_100.txt")

x4 = np.loadtxt("../results/K562_rep1_ENCFF285HUZ_TPM_lr_scaled_fpr_100.txt")
y4 = np.loadtxt("../results/K562_rep1_ENCFF285HUZ_TPM_lr_scaled_tpr_100.txt")
x5 = np.loadtxt("../results/K562_rep1_ENCFF285HUZ_TPM_matrix_weighted_counts_only_lr_fpr_100.txt")
y5 = np.loadtxt("../results/K562_rep1_ENCFF285HUZ_TPM_matrix_weighted_counts_only_lr_tpr_100.txt")
x6 = np.loadtxt("../results/K562_rep1_lr_fpr_100.txt")
y6 = np.loadtxt("../results/K562_rep1_lr_tpr_100.txt")

plt.figure()
plt.plot(x1, y1, color='b', label="GM12878 counts + weighted, AUC = 0.90")
plt.plot(x2, y2, color='g', label="GM12878 weighted only, AUC = 0.90")
plt.plot(x3, y3, color='r', label="GM12878 counts only, AUC = 0.90")

plt.plot(x4, y4, '--', color='b', label="K562 counts + weighted, AUC = 0.85")
plt.plot(x5, y5, '--', color='g', label="K562 weighted only, AUC = 0.84")
plt.plot(x6, y6, '--', color='r', label="K562 counts only, AUC = 0.85")
# plt.plot(x6, y6, '--', color='r', label="K562 no RNA-seq")
plt.xlabel("False positive rate")
plt.ylabel("True positive rate")
plt.legend()
plt.savefig("../figs/ROC.png", dpi=200)
# plt.show()