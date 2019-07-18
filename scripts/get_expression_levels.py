import pandas as pd
import numpy as np
#import mygene
import ast
import argparse

### This file extracts the expression levels from different .tsv files, so that we have multiple replicates per cell line. ###

parser = argparse.ArgumentParser()
parser.add_argument(
    '-c', '--converted-id',
    help="The .csv file containing the converted IDs (e.g. converted_id.csv)"
)
parser.add_argument(
    '-l', '--list-matrix-path',
    help="A list of paths to the target expression matrices"
)
parser.add_argument(
    '-o', '--output',
    help="Output file name"
)
args = parser.parse_args()
fn_id = args.converted_id
fn_list_matrix_path = args.list_matrix_path
fn_output = args.output

# get data from different TSV files
df_converted_id = pd.read_csv(fn_id, sep='\t', header=0, index_col=0)
df_list_matrix_path = pd.read_csv(fn_list_matrix_path, header=0, sep=',')

list_matrix = []
list_accession = [] #: extracts only accession names
for i in range(df_list_matrix_path.shape[0]):
    fn_expr_mat = df_list_matrix_path.iloc[i,0]
    df_expr_mat = pd.read_csv(fn_expr_mat, sep='\t')
    list_matrix.append(df_expr_mat)
    fn_expr_mat = fn_expr_mat.split('/')[-1]
    fn_expr_mat = fn_expr_mat.split('.')[0]
    list_accession.append(df_list_matrix_path.iloc[i,1] + '-' + fn_expr_mat)

tf_names = df_converted_id['MA']
gene_names = df_converted_id['gene_name']
ensembl_names = df_converted_id['ensembl-gene']
transcript_names = df_converted_id['ensembl-transcript']

df_out = pd.DataFrame()
df_out['MA'] = tf_names
df_out['Gene Name'] = gene_names
df_out['Ensembl ID'] = ensembl_names
#df_out['Transcript Name'] = transcript_names

for i, mat in enumerate(list_matrix):
    tpm = np.zeros((tf_names.shape[0]))
    #mat = np.matrix(mat)
    mat = mat.to_numpy()
    for n in range(0, len(ensembl_names)):
        for r in range(0, len(mat)):
            if pd.isna(ensembl_names[n]):
                print ('Skip {} : NA'.format(ensembl_names[n]))
                continue
            if mat[r][1].startswith(ensembl_names[n]):
                # check if only one transcript is present
                if len(transcript_names[n]) == 15:
                    # print (transcript_names[n])
                    if mat[r][0].startswith(transcript_names[n]):
                        tpm[n] += mat[r][5]
                else:
                    # get a list of transcripts that are present
                    transcripts = ast.literal_eval(transcript_names[n])
                    for t in transcripts:
                        if mat[r][0].startswith(t):
                            tpm[n] += mat[r][5]
    df_out[list_accession[i]] = tpm
    print ('Finished {}'.format(list_accession[i]))

##### Write to .CSV file ###
df_out.to_csv(fn_output, sep='\t', index=None)

# currently not including transcripts
# data = np.array([TF_names, gene_names, ensembl_names, tpm_1, tpm_2, tpm_3, tpm_4, tpm_5])
# data = data.T
# np.savetxt("../expression_levels.csv", data, delimiter='\t', header="MA\tgene_name\tensembl_name\tENCFF297CNO_TPM\tENCFF879WBJ_TPM\tENCFF285HUZ_TPM\tENCFF853TRI_TPM\tENCFF305QBE_TPM", fmt="%s")

