import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

from statsmodels.discrete.discrete_model import Logit
from statsmodels.stats.multitest import multipletests

import argparse

### This file analyzes the p-value of each of the features, by regressing non-zero columns with the labels and reporting the corrected p-values. ###

# input
parser = argparse.ArgumentParser()
parser.add_argument('Labels' , help="The .csv file containing labels (e.g. K562_rep1_labels.csv")
parser.add_argument('File_1' , help="The first feature .tsv file (e.g. K562_rep1_ENCFF285HUZ_TPM_matrix.tsv)")
parser.add_argument('File_2' , help="The second feature .tsv file (e.g. K562_rep1_ENCFF297CNO_TPM_matrix.tsv)")
args = vars(parser.parse_args())
labelsFile = args['ID_File']
file1 = args['File_1']
file2 = args['File_2']

# paths to the TPM and labels files
# labelsFile = '/Users/arun/Desktop/JHU/spring2019/CompGenomics/FinalProject/GWASSelection/ENCFF342EGB_labels.csv'
# file1 = '/Users/arun/Desktop/JHU/spring2019/CompGenomics/FinalProject/GWASSelection/ENCFF342EGB_ENCFF305QBE_TPM_matrix.tsv'
# file2 = '/Users/arun/Desktop/JHU/spring2019/CompGenomics/FinalProject/GWASSelection/ENCFF342EGB_ENCFF853TRI_TPM_matrix.tsv'

# labelsFile = '/home-4/xfeng17@jhu.edu/openchromatinprediction/K562_rep1_labels.csv'
# file1 = '/home-4/xfeng17@jhu.edu/with_expression/K562_rep1_ENCFF285HUZ_TPM_matrix.tsv'
# file2 = '/home-4/xfeng17@jhu.edu/with_expression/K562_rep1_ENCFF297CNO_TPM_matrix.tsv'

labels_matrix = np.loadtxt(labelsFile, delimiter = ',')
file1_matrix = np.loadtxt(file1, delimiter = '\t')
file2_matrix = np.loadtxt(file2, delimiter = '\t')

print(labels_matrix.shape)
# (109516,)
print(file1_matrix.shape)
# (109516, 1074)
print(file2_matrix.shape)
# (109516, 1074)

n_samples = labels_matrix.shape[0]

# we want to only use the non-zero columns, as otherwise, the p-values are meaningless
# file 1 zero/non-zero columns
nz1 = np.nonzero(np.any(file1_matrix != 0, axis = 0))[0]
z1 = np.nonzero(np.all(file1_matrix == 0, axis = 0))[0]
file1_nz = file1_matrix[:,nz1]
file1_z = file1_matrix[:,z1]

# file 2 zero/non-zero columns
nz2 = np.nonzero(np.any(file2_matrix != 0, axis = 0))[0]
z2 = np.nonzero(np.all(file2_matrix == 0, axis = 0))[0]
file2_nz = file2_matrix[:,nz2]
file2_z = file2_matrix[:,z2]

### Using the statsmodels Logit Model ###

# P-values and intercepts for each file
pVals1b = []
pVals2b = []
intercept1 = np.ones((n_samples,1))
intercept2 = np.ones((n_samples,1))

# file 1 - regress each feature with the labels, and store the p-value
for feature in range(0, len(nz1)):
	f = file1_nz[:,feature]
	np.reshape(f, (n_samples,1))
	X = np.column_stack((f, intercept1))
	Y = labels_matrix
	l = Logit(Y, X)
	r = l.fit()
	pVals1b.append(r.pvalues[0])

# carry out multiple hypothesis correction to obtain corrected value (we use FDR Benjamini Hochberg)
rej1b, pValsCorrected1b, aS1b, aB1b = multipletests(pVals1b, method = 'fdr_bh')


# file 2 - regress each feature with the labels, and store the p-value
for feature in range(0, len(nz2)):
	f = file2_nz[:,feature]
	np.reshape(f, (n_samples,1))
	X = np.column_stack((f, intercept1))
	Y = labels_matrix
	l = Logit(Y, X)
	r = l.fit()
	pVals2b.append(r.pvalues[0])

# carry out multiple hypothesis correction to obtain corrected value (we use FDR Benjamini Hochberg)
rej2b, pValsCorrected2b, aS2b, aB2b = multipletests(pVals2b, method = 'fdr_bh')


### Summary ###
# nz1 holds the non-zero features of file 1, pValsCorrected1b holds their p-values
# nz2 holds the non-zero features of file 2, pValsCorrected2b holds their p-values
nz1 = nz1.astype(int)
nz2 = nz2.astype(int)
file1_data = np.array([nz1, pValsCorrected1b])
file1_data = file1_data.T
file2_data = np.array([nz2, pValsCorrected2b])
file2_data = file2_data.T

# save the p-values to files
np.savetxt("../GWAS_pvalues_K562_rep1_ENCFF285HUZ.csv", file1_data, delimiter='\t', header="FeatureNumber\tPValue", fmt="%s")
np.savetxt("../GWAS_pvalues_K562_rep1_ENCFF297CNO.csv", file2_data, delimiter='\t', header="FeatureNumber\tPValue", fmt="%s")
