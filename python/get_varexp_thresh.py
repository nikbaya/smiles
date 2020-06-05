#!/usr/bin/env python3

import pandas as pd
import scipy.stats as stats
import numpy as np

wd = '/Users/nbaya/Downloads'

df = pd.read_csv('{wd}/clumped_gwas_metadata - Sheet1.tsv',sep='\t')
df['varexp_thresh'] = stats.chi2.isf(q=df['pval thresh'], df=1)/df['n_bar']
df.loc[df['pval thresh']==1,'varexp_thresh'] = np.nan
df.to_csv(f'{wd}/clumped_gwas_metadata.tsv',sep='\t',index=False)