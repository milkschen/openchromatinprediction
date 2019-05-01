import pandas as pd
import numpy as np
import mygene
import ast
import argparse

### This file extracts the expression levels from different .tsv files, so that we have multiple replicates per cell line. ###

parser = argparse.ArgumentParser()
parser.add_argument('ConvertedIDs' , help="The .csv file containing the converted IDs (e.g. converted_id.csv)")
parser.add_argument('EM1' , help="The .tsv file containing the first expression matrix (ENCFF297CNO.tsv)")
parser.add_argument('EM2' , help="The .tsv file containing the second expression matrix (ENCFF879WBJ.tsv)")
parser.add_argument('EM3' , help="The .tsv file containing the third expression matrix (ENCFF285HUZ.tsv)")
parser.add_argument('EM4' , help="The .tsv file containing the fourth expression matrix (ENCFF853TRI.tsv)")
parser.add_argument('EM5' , help="The .tsv file containing the fifth expression matrix (ENCFF305QBE.tsv)")
args = vars(parser.parse_args())
idFile = args['ConvertedIDs']
em1File = args['EM1']
em2File = args['EM2']
em3File = args['EM3']
em4File = args['EM4']
em5File = args['EM5']

# get data from different TSV files
df_converted_id = pd.read_csv(idFile, sep='\t', header=0, index_col=0)
df_expression_matrix_1 = pd.read_csv(em1File, sep='\t', header = 0, index_col = False)
df_expression_matrix_2 = pd.read_csv(em2File, sep='\t', header = 0, index_col = False)
df_expression_matrix_3 = pd.read_csv(em3File, sep='\t', header = 0, index_col = False)
df_expression_matrix_4 = pd.read_csv(em4File, sep='\t', header = 0, index_col = False)
df_expression_matrix_5 = pd.read_csv(em5File, sep='\t', header = 0, index_col = False)

# extract into numpy arrays
converted_id_matrix = df_converted_id.to_numpy()
expression_matrix_1 = df_expression_matrix_1.to_numpy()
expression_matrix_2 = df_expression_matrix_2.to_numpy()
expression_matrix_3 = df_expression_matrix_3.to_numpy()
expression_matrix_4 = df_expression_matrix_4.to_numpy()
expression_matrix_5 = df_expression_matrix_5.to_numpy()

# get the names
converted_id_headers = list(df_converted_id)
expression_matrix_1_headers = list(df_expression_matrix_1)
expression_matrix_2_headers = list(df_expression_matrix_2)
expression_matrix_3_headers = list(df_expression_matrix_3)
expression_matrix_4_headers = list(df_expression_matrix_4)
expression_matrix_5_headers = list(df_expression_matrix_5)

TF_names = converted_id_matrix[:,0]
gene_names = converted_id_matrix[:,2]
ensembl_names = converted_id_matrix[:,3]
transcript_names = converted_id_matrix[:,4]

# tpm matrices
tpm_1 = np.zeros((TF_names.shape[0]))
tpm_2 = np.zeros((TF_names.shape[0]))
tpm_3 = np.zeros((TF_names.shape[0]))
tpm_4 = np.zeros((TF_names.shape[0]))
tpm_5 = np.zeros((TF_names.shape[0]))

##### ENCFF297CNO.tsv #####
for n in range(0, len(ensembl_names)):
	# print(n)
	for r in range(0, len(expression_matrix_1)):
		if pd.isna(ensembl_names[n]):
			continue
		if expression_matrix_1[r][1].startswith(ensembl_names[n]):
			# check if only one transcript is present
			if len(transcript_names[n]) == 15:
				if expression_matrix_1[r][0].startswith(transcript_names[n]):
					tpm_1[n] += expression_matrix_1[r][5]
			else:
				# get a list of transcripts that are present
				transcripts = ast.literal_eval(transcript_names[n])
				for t in transcripts:
					if expression_matrix_1[r][0].startswith(t):
						tpm_1[n] += expression_matrix_1[r][5]

print(tpm_1)

##### ENCFF879WBJ.tsv #####
e1 = 0
e2 = 0
e3 = 0
for n in range(0, len(ensembl_names)):
	# print(n)
	for r in range(0, len(expression_matrix_2)):
		if pd.isna(ensembl_names[n]):
			continue
		# gene_id is first column, not second
		if expression_matrix_2[r][0].startswith(ensembl_names[n]):
			# only one transcript
			if len(transcript_names[n]) == 15:
				if transcript_names[n] in expression_matrix_2[r][1]:

					# check to make sure no additional transcripts in file
					if len(expression_matrix_2[r][1].split(',')) > 1:
						# print('Error1')
						tpm_2[n] += 0.0
						e1 += 1

					# match up exactly
					else:
						tpm_2[n] += expression_matrix_2[r][5]
				
				else:
					tpm_2[n] += 0.0
			else:
				# get a list of transcripts that are present
				transcripts = ast.literal_eval(transcript_names[n])
				# compare them to what was in the file

				# if one of the transcripts is not in the file
				test1 = 0
				for t in transcripts:
					if t not in expression_matrix_2[r][1]:
						test1 += 1
				if test1 > 0:
					# print('Error2')
					e2 += 1
					tpm_2[n] += 0.0
				else:

					# check if the file contains transcripts we don't have
					if len(expression_matrix_2[r][1].split(',')) > len(transcripts):
						# print('Error3')
						e3 += 1
						tpm_2[n] += 0.0

					else:
						tpm_2[n] += expression_matrix_2[r][5]
print(tpm_2)
print(e1,e2,e3)

##### ENCFF285HUZ.tsv #####
for n in range(0, len(ensembl_names)):
	# print(n)
	for r in range(0, len(expression_matrix_3)):
		if pd.isna(ensembl_names[n]):
			continue
		if expression_matrix_3[r][1].startswith(ensembl_names[n]):
			# check if only one transcript is present
			if len(transcript_names[n]) == 15:
				if expression_matrix_3[r][0].startswith(transcript_names[n]):
					tpm_3[n] += expression_matrix_3[r][5]
			else:
				# get a list of transcripts that are present
				transcripts = ast.literal_eval(transcript_names[n])
				for t in transcripts:
					if expression_matrix_3[r][0].startswith(t):
						tpm_3[n] += expression_matrix_3[r][5]

print(tpm_3)

##### ENCFF853TRI.tsv #####
for n in range(0, len(ensembl_names)):
	# print(n)
	for r in range(0, len(expression_matrix_4)):
		if pd.isna(ensembl_names[n]):
			continue
		if expression_matrix_4[r][1].startswith(ensembl_names[n]):
			# check if only one transcript is present
			if len(transcript_names[n]) == 15:
				if expression_matrix_4[r][0].startswith(transcript_names[n]):
					tpm_4[n] += expression_matrix_4[r][5]
			else:
				# get a list of transcripts that are present
				transcripts = ast.literal_eval(transcript_names[n])
				for t in transcripts:
					if expression_matrix_4[r][0].startswith(t):
						tpm_4[n] += expression_matrix_4[r][5]

print(tpm_4)

##### ENCFF305QBE.tsv #####
for n in range(0, len(ensembl_names)):
	# print(n)
	for r in range(0, len(expression_matrix_5)):
		if pd.isna(ensembl_names[n]):
			continue
		if expression_matrix_5[r][1].startswith(ensembl_names[n]):
			# check if only one transcript is present
			if len(transcript_names[n]) == 15:
				if expression_matrix_5[r][0].startswith(transcript_names[n]):
					tpm_5[n] += expression_matrix_5[r][5]
			else:
				# get a list of transcripts that are present
				transcripts = ast.literal_eval(transcript_names[n])
				for t in transcripts:
					if expression_matrix_5[r][0].startswith(t):
						tpm_5[n] += expression_matrix_5[r][5]

print(tpm_5)

##### Write to .CSV file ###

# currently not including transcripts
data = np.array([TF_names, gene_names, ensembl_names, tpm_1, tpm_2, tpm_3, tpm_4, tpm_5])
data = data.T
np.savetxt("../expression_levels.csv", data, delimiter='\t', header="MA\tgene_name\tensembl_name\tENCFF297CNO_TPM\tENCFF879WBJ_TPM\tENCFF285HUZ_TPM\tENCFF853TRI_TPM\tENCFF305QBE_TPM", fmt="%s")

