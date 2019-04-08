import pandas as pd
import mygene

df_id = pd.read_csv('../human_pwm_ids_sorted.txt', sep=' ', header=None, index_col=0)
df_id.index = range(0, df_id.shape[0])
df_id.columns = ['MA', 'gene_name']
l = []
for d in df_id['gene_name']:
    if d.count(':') == 0 and d.count('(') == 0:
        l.append(d)
    else:
        d = d.split(':')[0]
        d = d.split('(')[0]
        l.append(d)
df_id['gene_name2'] = l

mg = mygene.MyGeneInfo()

ensbl_gene_id = []
ensbl_transcript_id = []
# log = []
for i, d in enumerate(df_id['gene_name2']):
    qry = mg.query(d, scopes="symbol", fields=["ensembl"], species="human", verbose=False)
    if len(qry) == 0:
        ensbl_gene_id.append('')
        ensbl_transcript_id.append('')
        # log.append('empty query')
        print (i, 'empty query')
        continue
    hits = qry['hits']
    if len(hits) == 0:
        ensbl_gene_id.append('')
        ensbl_transcript_id.append('')
        # log.append('empty hits')
        print (i, 'empty hits')
        continue
    ensbl = hits[0]['ensembl']
    try:
        ensbl_gene_id.append(ensbl['gene'])
        ensbl_transcript_id.append(ensbl['transcript'])
        # log.append('1')
    except:
        print (i, 'multiple ensbl')
        ensbl_gene_id.append(ensbl[0]['gene'])
        ensbl_transcript_id.append(ensbl[0]['transcript'])
        # log.append('{0} ensbls'.format(len(ensbl)))

df_id['ensembl-gene'] = ensbl_gene_id
df_id['ensembl-transcript'] = ensbl_transcript_id
# df_id['ensembl-log'] = log


# df_rna = pd.read_csv('../ENCFF297CNO.tsv', sep='\t')