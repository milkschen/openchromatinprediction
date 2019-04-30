import pandas as pd
import mygene

### This file does the conversion from MA IDs (in the JASPAR database) to Ensembl IDs###

# read in the PWM IDs
df_id = pd.read_csv('../human_pwm_ids_sorted.txt', sep=' ', header=None, index_col=0)
df_id.index = range(0, df_id.shape[0])
# get the relevant columns - MA and Gene Name
df_id.columns = ['MA', 'gene_name']


l = []
for d in df_id['gene_name']:
	if d.count(':') == 0 and d.count('(') == 0:
		l.append(d)
	else:
		d = d.split(':')[0]
		d = d.split('(')[0]
		l.append(d)
# filtered gene names
df_id['gene_name2'] = l


# we are using the mygene package to convert gene IDs from MA to the Ensembl IDs
mg = mygene.MyGeneInfo()

# get the ensembl gene and transcript IDs
ensbl_gene_id = []
ensbl_transcript_id = []

# for each ID, we will query the mygene database, and obtain the relevant translated ID
for i, d in enumerate(df_id['gene_name2']):
	qry = mg.query(d, scopes="symbol", fields=["ensembl"], species="human", verbose=False)

	# no matches found
	if len(qry) == 0:
		ensbl_gene_id.append('')
		ensbl_transcript_id.append('')
		print (i, 'empty query')
		continue
	hits = qry['hits']

	# no matches found
	if len(hits) == 0:
		ensbl_gene_id.append('')
		ensbl_transcript_id.append('')
		print (i, 'empty hits')
		continue
	ensbl = hits[0]['ensembl']

	try:
		ensbl_gene_id.append(ensbl['gene'])
		ensbl_transcript_id.append(ensbl['transcript'])
	except:
		print (i, 'multiple ensbl')
		ensbl_gene_id.append(ensbl[0]['gene'])
		ensbl_transcript_id.append(ensbl[0]['transcript'])

# store the IDs
df_id['ensembl-gene'] = ensbl_gene_id
df_id['ensembl-transcript'] = ensbl_transcript_id

# export CSV
df_id.to_csv('../converted_id.csv', sep='\t')