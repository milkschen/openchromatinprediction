#!/usr/bin/env python

"""
Trains a open/closed chromatin classifier on the given matrix of PWM matches to open chromation regions
and assesses its performance.

Xinyu Feng, April 4 2018

May 1, 2019: update to scale features to unit variance and zero mean.

Example usage: ./11-classifier.py -r 10-test/00-test_all_matrix.tsv 50 human_pwm_ids_sorted.txt lr 11-test/result
Example 2: ./11-classifier.py -r -p 4 training_data/A549/noNNN/ENCFF045PYX_rep5_all_matrix.tsv 72082 human_pwm_ids_sorted.txt rf 11-test/A549_rep5_RF
RNA-seq example: ./11-classifier2.py -r ENCFF342EGB_ENCFF297CNO_TPM_matrix.tsv 56659 ../human_pwm_ids_sorted.txt lr ENCFF342EGB_ENCFF297CNO_TPM_lr
"""

import argparse

# Parse arguments

parser = argparse.ArgumentParser(description='Train a chromatin accessibility classifier.')
parser.add_argument('matrix' , help="<mat.tsv> matrix of training data")
parser.add_argument('n_true', type=int, help="<int> number of true labels in the matrix")
parser.add_argument('motif_list', help="<pwm_ids.txt> names of pwms in the same order as the matrix")
parser.add_argument('classifier', help='< lr | rf > train a logistic regression or random forest classifier')
parser.add_argument('out_prefix', help="prefix for output")

parser.add_argument('-r', '--roc', action='store_true', help='generate ROC curve')
parser.add_argument('-p', '--threads', type=int, help='[int=4] number of threads to use.')

args = vars(parser.parse_args())

from sklearn.preprocessing import normalize
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression 
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

from sklearn.preprocessing import StandardScaler

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

"""
Evaluate performance of clf, save results to log file and as figures.

input  clf: trained LogisticRegressionCV or RandomeForest instance
"""
def performance(clf, X_test, y_test, out_prefix, will_plot, log):

	test_score = clf.score(X_test, y_test)
	y_score = clf.predict_proba(X_test)
	auc = roc_auc_score(y_test, y_score[:, 1])
	fpr, tpr, thresholds = roc_curve(y_test, y_score[:, 1])

	np.savetxt(out_prefix + '_fpr.txt', fpr)
	np.savetxt(out_prefix + '_tpr.txt', tpr)

	if will_plot:

		plt.figure()
		plt.plot(fpr, tpr)
		plt.xlabel("False positive rate")
		plt.ylabel("True positive rate")
		plt.savefig(out_prefix + '_roc.png')
		plt.close()

	if (type(clf) == LogisticRegressionCV):
		log.write("\nLogistic regression classifier:\n\n")
		log.write("optimized C = %.2f\ntest set score = %.2f\nauc = %.2f\n" %(clf.C_, test_score, auc))

	elif (type(clf) == RandomForestClassifier):
		log.write("\nRandom forest classifier:\n\n")
		log.write("n_estimators = %d\ntest set score = %.2f\nauc = %.2f\n" %(len(clf.estimators_), test_score, auc))

	return

""" 
Save pwms ranked by clf parameters.

input ranks: 1-D array
"""
def save_ranks(clf, ranks, motif_list, out_prefix):

	sort_ids = np.flip(np.argsort(ranks), 0)
	n_motifs = motif_list.shape[0]
	constant_row = pd.DataFrame([["const", "const", "const"]], columns=motif_list.columns)
	motif_list_with_const = constant_row.append(motif_list)

	# DNase-seq only or RNA-seq only
	if len(ranks) == n_motifs+1:
		motif_list_with_const['rank'] = ranks.T

	# DNase-seq and RNA-seq:
	else:
		motif_list_with_const = motif_list_with_const.append(motif_list)
		for i in range(n_motifs+1, 2*n_motifs+1):
			motif_list_with_const.iloc[i, 1] = motif_list_with_const.iloc[i, 1] + "*TPM"
		motif_list_with_const['rank'] = ranks.T

	to_save = motif_list_with_const.iloc[sort_ids, [1, 2, 3]]
	to_save.to_csv(out_prefix + '_ranks.tsv', sep='\t',header=False, index=False)

	return

matrix_filename = args['matrix']
n_true = args['n_true']
out_prefix = args['out_prefix']
will_plot = args['roc']
motif_list = pd.read_csv(args['motif_list'], sep=' ', header=None)

clf_type = args['classifier']
if args['threads']:
	p = args['threads']
else:
	p = 4

log = open(out_prefix + '.log', 'w')
log.write(time.asctime(time.localtime()) + '\n')
log.write('Matrix file: %s\n' %matrix_filename)
log.write('Number of true labels: %d\n' %n_true)
if clf_type == 'rf':
	clf_type_print = 'random forest'
elif clf_type == 'lr':
	clf_type_print = 'logistic regression'
log.write('Training a %s classifier...\n' %clf_type_print)

# Loading data 

X = np.loadtxt(matrix_filename, delimiter='\t')
# X = normalize(X)
X_with_const = np.ones((X.shape[0], X.shape[1]+1))
X_with_const[:, 1:] = X

y = np.zeros((X.shape[0],), dtype=np.uint8)
for i in range(n_true):
	y[i] = 1
X_train, X_test, y_train, y_test = train_test_split(X_with_const, y, test_size=0.33)

scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# Train classifier

if clf_type == 'lr':

	clf = LogisticRegressionCV(cv=5)
	clf.fit(X_train, y_train)
	ranks = np.reshape(clf.coef_, len(clf.coef_[0]))

elif clf_type == 'rf':

	clf = RandomForestClassifier(n_estimators=20, n_jobs=p)
	clf.fit(X_train, y_train)
	ranks = clf.feature_importances_

performance(clf, X_test, y_test, out_prefix, will_plot, log)
save_ranks(clf, ranks, motif_list, out_prefix)

log.close()
