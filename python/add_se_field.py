#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 16:41:29 2020

Get standard errors for arabidopsis data

@author: nbaya
"""

from glob import glob
import pandas as pd
from scipy import stats

smiles_dir = '/Users/nbaya/Documents/lab/smiles/data'

files = []

for pheno in ['GxE','FSHL','FSLL','SLHL','SLLL','t50HL','t50LL']:
    files += glob(smiles_dir+f'/{pheno}*.tsv.gz')
    
for file in files:
    df = pd.read_csv(file,sep='\t',compression='gzip')
    df = df.rename(columns={'Chr':'chr',
                            'Pos':'pos',
                            'Pval':'pval'})    
    pre_ct = df.shape[0]
    df = df[df.MAF>0]
    print(f'{pre_ct-df.shape[0]} variants removed due to MAF=0')
    
    pre_ct = df.shape[0]
    df = df.dropna(axis=0, subset=['beta','pval','MAF'])
    print(f'{pre_ct-df.shape[0]} variants removed due to missing beta, pval, maf')
    
    df['se'] = abs(df.beta/stats.norm.ppf(df.pval/2))
    
    df['eaf'] = df.AC_1/(df.AC_0+df.AC_1) # AC_1 is effect allele
    
    df['raf'] = df.eaf*(df.beta>0) + (1-df.eaf)*(df.beta<0)
    
    df['rbeta'] = abs(df.beta)
    
    df['var_exp'] = 2*(df.eaf)*(1-df.eaf)*df.beta**2
    
    df = df[['chr','pos','eaf','raf','beta','se','pval']]
    
    prefix = file.split('.')[0]
    
    df.to_csv(f'{prefix}.cleaned.tsv.gz',index=False,sep='\t',compression='gzip')