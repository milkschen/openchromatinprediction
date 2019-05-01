# openchromatinprediction

# openchromatinprediction

## Step-by-step walkthrough of our project
### Download the data:
```shell
wget https://www.encodeproject.org/files/ENCFF285HUZ/@@download/ENCFF285HUZ.tsv

wget https://www.encodeproject.org/files/ENCFF297CNO/@@download/ENCFF297CNO.tsv

wget https://www.encodeproject.org/files/ENCFF853TRI/@@download/ENCFF853TRI.tsv

wget https://www.encodeproject.org/files/ENCFF305QBE/@@download/ENCFF305QBE.tsv

wget https://www.encodeproject.org/files/ENCFF879WBJ/@@download/ENCFF879WBJ.tsv
```
### TODO - add intermediate steps

### Convert Gene IDs
Convert JASPAR Gene IDs into Ensembl IDs:
```shell
python3 convert_gene_id.py human_pwm_ids_sorted.txt converted_id.csv
```

### Get expression levels
Get expression levels from the replicates:
```shell
python3 get_expression_levels.py converted_id.csv ENCFF297CNO.tsv ENCFF879WBJ.tsv ENCFF285HUZ.tsv ENCFF853TRI.tsv ENCFF305QBE.tsv
```

### TODO - add intermediate steps


### Significance Analysis


### Plotting

Plotting the expression levels:
```shell
python3 plot_expression.py expression_levels.csv
```
Plotting the expression heat map:
```shell
python3 plot_expression_heatmap.py expression_levels.csv
```
