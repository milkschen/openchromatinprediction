import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

from statsmodels.discrete.discrete_model import Logit
from statsmodels.stats.multitest import multipletests

# paths to the files
labelsFile = '/Users/arun/Desktop/JHU/spring2019/CompGenomics/FinalProject/GWASSelection/ENCFF342EGB_labels.csv'
file1 = '/Users/arun/Desktop/JHU/spring2019/CompGenomics/FinalProject/GWASSelection/ENCFF342EGB_ENCFF305QBE_TPM_matrix.tsv'
file2 = '/Users/arun/Desktop/JHU/spring2019/CompGenomics/FinalProject/GWASSelection/ENCFF342EGB_ENCFF853TRI_TPM_matrix.tsv'

labels_matrix = np.loadtxt(labelsFile, delimiter = ',')
file1_matrix = np.loadtxt(file1, delimiter = '\t')
file2_matrix = np.loadtxt(file2, delimiter = '\t')

print(labels_matrix.shape)
# (109516,)
print(file1_matrix.shape)
# (109516, 1074)
print(file2_matrix.shape)
# (109516, 1074)

# file 1 zero/non-zero columns
nz1 = np.nonzero(np.any(file1_matrix != 0, axis = 0))[0]
z1 = np.nonzero(np.all(file1_matrix == 0, axis = 0))[0]

file1_nz = file1_matrix[:,nz1]
file1_z = file1_matrix[:,z1]

nz2 = np.nonzero(np.any(file2_matrix != 0, axis = 0))[0]
z2 = np.nonzero(np.all(file2_matrix == 0, axis = 0))[0]

file2_nz = file2_matrix[:,nz2]
file2_z = file2_matrix[:,z2]


#### sm.OLS model ###

# pVals1a = []
# pVals2a = []
# intercept1 = np.ones((109516,1))
# intercept2 = np.ones((109516,1))

# # file 1
# for feature in range(0, len(nz1)):
# 	f = file1_nz[:,feature]
# 	np.reshape(f, (109516,1))
# 	X = np.column_stack((f, intercept1))
# 	Y = labels_matrix
# 	model = sm.OLS(Y,X)
# 	results = model.fit()
# 	pVals1a.append(results.pvalues[0])

# # Getting these warnings
# # /Users/arun/anaconda3/lib/python3.6/site-packages/statsmodels/base/model.py:1100: RuntimeWarning: invalid value encountered in true_divide
# #   return self.params / self.bse
# # /Users/arun/anaconda3/lib/python3.6/site-packages/scipy/stats/_distn_infrastructure.py:877: RuntimeWarning: invalid value encountered in greater
# #   return (self.a < x) & (x < self.b)
# # /Users/arun/anaconda3/lib/python3.6/site-packages/scipy/stats/_distn_infrastructure.py:877: RuntimeWarning: invalid value encountered in less
# #   return (self.a < x) & (x < self.b)
# # /Users/arun/anaconda3/lib/python3.6/site-packages/scipy/stats/_distn_infrastructure.py:1831: RuntimeWarning: invalid value encountered in less_equal
# #   cond2 = cond0 & (x <= self.a)

# rej1a, pValsCorrected1a, aS1a, aB1a = multipletests(pVals1a, method = 'fdr_bh')

# # Get this warning
# # /Users/arun/anaconda3/lib/python3.6/site-packages/statsmodels/stats/multitest.py:320: RuntimeWarning: invalid value encountered in less_equal
# #   reject = pvals_sorted <= ecdffactor*alpha

# # file 2
# for feature in range(0, len(nz2)):
# 	f = file2_nz[:,feature]
# 	np.reshape(f, (109516,1))
# 	X = np.column_stack((f, intercept2))
# 	Y = labels_matrix
# 	model = sm.OLS(Y,X)
# 	results = model.fit()
# 	pVals2a.append(results.pvalues[0])

# # Warnings
# # /Users/arun/anaconda3/lib/python3.6/site-packages/statsmodels/base/model.py:1100: RuntimeWarning: invalid value encountered in true_divide
# #   return self.params / self.bse
# # /Users/arun/anaconda3/lib/python3.6/site-packages/scipy/stats/_distn_infrastructure.py:877: RuntimeWarning: invalid value encountered in greater
# #   return (self.a < x) & (x < self.b)
# # /Users/arun/anaconda3/lib/python3.6/site-packages/scipy/stats/_distn_infrastructure.py:877: RuntimeWarning: invalid value encountered in less
# #   return (self.a < x) & (x < self.b)
# # /Users/arun/anaconda3/lib/python3.6/site-packages/scipy/stats/_distn_infrastructure.py:1831: RuntimeWarning: invalid value encountered in less_equal
# #   cond2 = cond0 & (x <= self.a)

# rej2a, pValsCorrected2a, aS2a, aB2a = multipletests(pVals2a, method = 'fdr_bh')

# # /Users/arun/anaconda3/lib/python3.6/site-packages/statsmodels/stats/multitest.py:320: RuntimeWarning: invalid value encountered in less_equal
# #   reject = pvals_sorted <= ecdffactor*alpha




### Logit Model ###

pVals1b = []
pVals2b = []
intercept1 = np.ones((109516,1))
intercept2 = np.ones((109516,1))


# file 1
for feature in range(0, len(nz1)):
	f = file1_nz[:,feature]
	np.reshape(f, (109516,1))
	X = np.column_stack((f, intercept1))
	Y = labels_matrix
	l = Logit(Y, X)
	r = l.fit()
	pVals1b.append(r.pvalues[0])

# Get this error
# LinAlgError: Singular matrix

rej1b, pValsCorrected1b, aS1b, aB1b = multipletests(pVals1b, method = 'fdr_bh')


# file 2
for feature in range(0, len(nz2)):
	f = file2_nz[:,feature]
	np.reshape(f, (109516,1))
	X = np.column_stack((f, intercept1))
	Y = labels_matrix
	l = Logit(Y, X)
	r = l.fit()
	pVals2b.append(r.pvalues[0])

# Get this error
# LinAlgError: Singular matrix

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

np.savetxt("../GWAS_pvalues_ENCFF342EGB_ENCFF305QBE.csv", file1_data, delimiter='\t', header="FeatureNumber\tPValue", fmt="%s")
np.savetxt("../GWAS_pvalues_ENCFF342EGB_ENCFF853TRI.csv", file2_data, delimiter='\t', header="FeatureNumber\tPValue", fmt="%s")

### Manhattan plot ###
# TODO
