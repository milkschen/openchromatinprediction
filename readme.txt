Prediction of Pioneer Transcription Factors Using Chromatin Accessibility and Gene Expression

Nae-Chyun Chen, Arun Das, Xinyu Feng

This is our project README - you can find a markdown version in "README.md".


##########################

In the next section, we specify how to download the relevant data files. The data for this project was obtained from the ENCODE project and from the JASPAR dataset. The code to obtain the latter is currently not provided (as it is made up of numerous path sensitive bash scripts) - instead, here is the link to the data:

TODO

###########################

Here is a step-by-step walkthrough of our project:

Download the data:
wget https://www.encodeproject.org/files/ENCFF285HUZ/@@download/ENCFF285HUZ.tsv
wget https://www.encodeproject.org/files/ENCFF297CNO/@@download/ENCFF297CNO.tsv
wget https://www.encodeproject.org/files/ENCFF853TRI/@@download/ENCFF853TRI.tsv
wget https://www.encodeproject.org/files/ENCFF305QBE/@@download/ENCFF305QBE.tsv
wget https://www.encodeproject.org/files/ENCFF879WBJ/@@download/ENCFF879WBJ.tsv

Convert gene JASPAR IDs into transcript and gene Ensembl IDs:
python3 convert_gene_id.py human_pwm_ids_sorted.txt converted_id.csv

Get expression levels from the replicates:
python3 get_expression_levels.py converted_id.csv ENCFF297CNO.tsv ENCFF879WBJ.tsv ENCFF285HUZ.tsv ENCFF853TRI.tsv ENCFF305QBE.tsv

Train a logistic classifier
./11-classifier2.py -r ENCFF342EGB_ENCFF305QBE_TPM_matrix.tsv 56659 ../human_pwm_ids_sorted.txt lr ENCFF342EGB_ENCFF305QBE_TPM_lr_scaled
./11-classifier2.py -r ENCFF342EGB_ENCFF853TRI_TPM_matrix.tsv 56659 ../human_pwm_ids_sorted.txt lr ENCFF342EGB_ENCFF853TRI_TPM_lr_scaled
./11-classifier2.py -r K562_rep1_ENCFF285HUZ_TPM_matrix.tsv 263576 ../human_pwm_ids_sorted.txt lr K562_rep1_ENCFF285HUZ_TPM_lr_scaled
./11-classifier2.py -r K562_rep1_ENCFF297CNO_TPM_matrix.tsv 263576 ../human_pwm_ids_sorted.txt lr K562_rep1_ENCFF297CNO_TPM_lr_scaled


Get the corrected p-values of each of the features:
python3 GWAS_selection.py K562_rep1_labels.csv K562_rep1_ENCFF285HUZ_TPM_matrix.tsv K562_rep1_ENCFF297CNO_TPM_matrix.tsv
python3 GWAS_selection.py ENCFF342EGB_labels.csv ENCFF342EGB_ENCFF305QBE_TPM_matrix.tsv ENCFF342EGB_ENCFF853TRI_TPM_matrix.tsv

###########################

Here is how we generated our plots:

Plotting the expression levels (Figure 2a):
python3 plot_expression.py expression_levels.csv

Plotting the expression heat map:
python3 plot_expression_heatmap.py expression_levels.csv

Plotting the correlation heatmap (Figure 2b):
python3 plot_correlation_heatmap.py ENCFF342EGB_ENCFF853TRI_TPM_matrix.tsv

Plotting p-values (Figure 3a):
python3 plot_pvalues.py GWAS_pvalues_ENCFF342EGB_ENCFF853TRI.csv

Plotting coefficients:
python3 plot_coeff.py ENCFF342EGB_ENCFF853TRI_TPM_lr_ranks.tsv
