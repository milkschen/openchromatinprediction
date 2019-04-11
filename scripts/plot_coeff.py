import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

df = pd.read_csv('ENCFF342EGB_ENCFF853TRI_TPM_lr_ranks.tsv', sep='\t', header=None)
co_const = []
co_cntonly = []
co_cnt_and_tpm = []
for i, c in enumerate(df[2]):
    if df[0][i] == 'const':
        co_const.append((df[0][i], df[1][i], c))
    elif df[0][i].count('TPM') > 0:
        co_cnt_and_tpm.append((df[0][i], df[1][i], c))
    else:
        co_cntonly.append((df[0][i], df[1][i], c))
df_cntonly = pd.DataFrame(co_cntonly)
df_cntonly.columns = ['id', 'gene', 'coeff']

df_cnt_and_tmp = pd.DataFrame(co_cnt_and_tpm)
ma_id = []
for i in df_cnt_and_tmp[0]:
    ma_id.append(i[:i.find('*')])
df_cnt_and_tmp['id'] = ma_id
df_cnt_and_tmp.columns = ['raw_id', 'gene', 'coeff', 'id']

df_merged = pd.merge(df_cntonly, df_cnt_and_tmp, on=['id'])

print (pearsonr(df_merged['coeff_y'], df_merged['coeff_x']))

x = np.arange(df_merged.shape[0])
plt.plot(x, df_merged['coeff_y'], label='with expr')
plt.plot(x, df_merged['coeff_x'], label='cnt only')
plt.legend()
plt.ylabel('coefficient')
plt.xlabel('TF')
plt.show()
plt.clf()