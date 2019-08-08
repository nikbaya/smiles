#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 08:20:11 2019

For cleaning up summary statistics before plotting

@author: nbaya
"""

import pandas as pd
import numpy as np

wd_data = "/Users/nbaya/Documents/lab/smiles/data/"


## Cleaning up T2D data and removing unnecessary columns to save space
#t2d = pd.read_csv(wd_data+'Mahajan.NatGenet2018b.T2Dbmiadj.European.txt.gz',compression='gzip',delimiter='\t')
#t2d = t2d.rename(columns={'Chr':'chr','Beta':'beta','Pvalue':'pval','Neff':'n_complete_samples'})
#t2d['chr'] = t2d['chr'].astype(str)
#t2d = t2d[['chr','EAF','beta','pval','n_complete_samples']]
#t2d.to_csv(wd_data+'Mahajan.NatGenet2018b.T2Dbmiadj.European.txt.gz',compression='gzip',sep='\t',index=False)

## Cleaning up SCZ1 data and removing unnecessary columns to save space
scz = pd.read_csv(wd_data+'pgc.scz.full.2012-04.txt.gz',compression='gzip',delim_whitespace=True)
scz = scz.rename(columns={'hg18chr':'chr','or':'beta'})
scz['n_complete_samples'] = 150064 #based on abstract (https://www.ncbi.nlm.nih.gov/pubmed/25056061)
scz = scz[scz.CEUaf!='.'] #filter out SNPs with missing CEUaf
scz['beta'] = np.log10(scz['beta']) # convert OR to correct scale for effect size
scz.loc[scz.beta>0,'EAF'] = scz.loc[scz.beta>0,'CEUaf'].astype(float)
scz.loc[scz.beta<0,'EAF'] = 1-scz.loc[scz.beta<0,'CEUaf'].astype(float)
scz = scz[['chr','EAF','beta','pval','n_complete_samples']]
scz.to_csv(wd_data+'pgc.scz.full.2012-04.tsv.gz',compression='gzip',sep='\t',index=False)

## Clean up IBD data
ibd = pd.read_csv(wd_data+'EUR.IBD.gwas_info03_filtered.assoc.gz',compression='gzip',delim_whitespace=True)
n_cas = int([x for x in ibd.columns.values if 'FRQ_A' in x][0].split('_')[2])
n_con = int([x for x in ibd.columns.values if 'FRQ_U' in x][0].split('_')[2])
ibd = ibd.rename(columns={f'FRQ_U_{n_con}':'EAF','CHR':'chr','P':'pval','BP':'pos'})
ibd['beta'] = np.log10(ibd.OR)
ibd['n'] = n_cas + n_con
ibd = ibd[['chr','EAF','beta','pval','n','pos']]
ibd.to_csv(wd_data+'EUR.IBD.gwas_info03_filtered.assoc.tsv.gz',compression='gzip',sep='\t',index=False)

## Clean up CD (Crohn's disease) data
df = pd.read_csv(wd_data+'EUR.CD.gwas_info03_filtered.assoc.gz',compression='gzip',delim_whitespace=True)
n_cas = int([x for x in df.columns.values if 'FRQ_A' in x][0].split('_')[2])
n_con = int([x for x in df.columns.values if 'FRQ_U' in x][0].split('_')[2])
df = df.rename(columns={f'FRQ_U_{n_con}':'EAF','CHR':'chr','P':'pval','BP':'pos'})
df['beta'] = np.log10(df.OR)
df['n'] = n_cas + n_con
df = df[['chr','EAF','beta','pval','n','pos']]
df.to_csv(wd_data+'EUR.CD.gwas_info03_filtered.assoc.tsv.gz',compression='gzip',sep='\t',index=False)

## Clean up UC (ulcerative colitis) data
df = pd.read_csv(wd_data+'EUR.UC.gwas_info03_filtered.assoc.gz',compression='gzip',delim_whitespace=True)
n_cas = int([x for x in df.columns.values if 'FRQ_A' in x][0].split('_')[2])
n_con = int([x for x in df.columns.values if 'FRQ_U' in x][0].split('_')[2])
df = df.rename(columns={f'FRQ_U_{n_con}':'EAF','CHR':'chr','P':'pval','BP':'pos'})
df['beta'] = np.log10(df.OR)
df['n'] = n_cas + n_con
df = df[['chr','EAF','beta','pval','n','pos']]
df.to_csv(wd_data+'EUR.UC.gwas_info03_filtered.assoc.tsv.gz',compression='gzip',sep='\t',index=False)

## Clean up daner PGC SCZ data
df = pd.read_csv(wd_data+'daner_PGC_SCZ43_mds9.gz.hq2.gz',compression='gzip',delim_whitespace=True)
n_cas = int([x for x in df.columns.values if 'FRQ_A' in x][0].split('_')[2])
n_con = int([x for x in df.columns.values if 'FRQ_U' in x][0].split('_')[2])
df = df.rename(columns={f'FRQ_U_{n_con}':'EAF','CHR':'chr','P':'pval','BP':'pos'})
df['beta'] = np.log10(df.OR)
df['n'] = n_cas + n_con
df = df[['chr','EAF','beta','pval','n','pos']]
df.to_csv(wd_data+'daner_PGC_SCZ43_mds9.tsv.gz',compression='gzip',sep='\t',index=False)



df.shape

df.OR.max()


plt.hist(df.OR,50)
plt.xlabel('odds ratio')
plt.ylabel('density')

plt.hist(df.EAF,50)
plt.xlabel('EAF')
plt.ylabel('density')
plt.title('EAF distribution')

plt.hist(df[df.OR>1].EAF,50)
plt.xlabel('EAF')
plt.ylabel('density')
plt.title('EAF distribution for SNPs with OR>1')

plt.hist(df[df.OR<1].EAF,50)
plt.xlabel('EAF')
plt.ylabel('density')
plt.title('EAF distribution for SNPs with OR<1')


pval_threshold=5e-8
plt.hist(df[(df.OR>1)&(df.pval<pval_threshold)].EAF,50)
plt.xlabel('EAF')
plt.ylabel('density')
plt.title(f'EAF distribution for SNPs with OR>1 and pval<{pval_threshold}')

pval_threshold=5e-8
plt.hist(1-df[(df.OR<1)&(df.pval<pval_threshold)].EAF,50)
plt.xlabel('EAF')
plt.ylabel('density')
plt.title(f'EAF distribution for SNPs with OR<1 and pval<{pval_threshold}')



plt.hist(ss0[ss0.alt_a].alt_af)
plt.title('RAF distribution for UKB SCZ data\n(without low confidence variants)')
#plt.hist(ss0[(ss0.pval<1)&(~ss0.low_confidence_variant)].raf.dropna())
plt.title('Odds ratio distribution for UKB SCZ data\n(without low confidence variants)')
plt.hist(ss0[(~ss0.low_confidence_variant)].OR,50)
plt.xlabel('odds ratio')


plt.title('Odds ratio distribution for UKB SCZ data\n(with low confidence variants)')
plt.hist(ss0.beta.dropna(),50)
plt.xlabel('odds ratio')



plt.xlabel('RAF')
plt.ylabel('density')
plt.xlim([0.99,1.01])

ss0['OR'] = np.exp(ss0.beta)

plt.hist(df[df.OR<1].EAF,50)
plt.hist(df[df.pval<1e-8].EAF,50)


ss0.columns.values


