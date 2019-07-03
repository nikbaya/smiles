#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 08:20:11 2019

For cleaning up summary statistics before plotting

@author: nbaya
"""

import pandas as pd
import numpy as np

wd = "/Users/nbaya/Documents/lab/smiles/"


## Cleaning up T2D data and removing unnecessary columns to save space
#t2d = pd.read_csv(wd+'data/'+'Mahajan.NatGenet2018b.T2Dbmiadj.European.txt.gz',compression='gzip',delimiter='\t')
#t2d = t2d.rename(columns={'Chr':'chr','Beta':'beta','Pvalue':'pval','Neff':'n_complete_samples'})
#t2d['chr'] = t2d['chr'].astype(str)
#t2d = t2d[['chr','EAF','beta','pval','n_complete_samples']]
#t2d.to_csv(wd+'data/'+'Mahajan.NatGenet2018b.T2Dbmiadj.European.txt.gz',compression='gzip',sep='\t',index=False)

## Cleaning up SCZ1 data and removing unnecessary columns to save space
scz = pd.read_csv(wd+'data/'+'pgc.scz.full.2012-04.txt.gz',compression='gzip',delim_whitespace=True)
scz = scz.rename(columns={'hg18chr':'chr','or':'beta'})
scz['n_complete_samples'] = 150064 #based on abstract (https://www.ncbi.nlm.nih.gov/pubmed/25056061)
scz = scz[scz.CEUaf!='.'] #filter out SNPs with missing CEUaf
scz['beta'] = np.log10(scz['beta']) # convert OR to correct scale for effect size
scz.loc[scz.beta>0,'EAF'] = scz.loc[scz.beta>0,'CEUaf'].astype(float)
scz.loc[scz.beta<0,'EAF'] = 1-scz.loc[scz.beta<0,'CEUaf'].astype(float)
scz = scz[['chr','EAF','beta','pval','n_complete_samples']]
scz.to_csv(wd+'data/'+'pgc.scz.full.2012-04.tsv.gz',compression='gzip',sep='\t',index=False)

## Clean up IBD data
ibd = pd.read_csv(wd+'data/'+'EUR.IBD.gwas_info03_filtered.assoc.gz',compression='gzip',delim_whitespace=True)
ibd = ibd.rename(columns={'FRQ_U_21770':'EAF','CHR':'chr','P':'pval'})
ibd['beta'] = np.log10(ibd.OR)
ibd['n'] = 12882+21770
ibd = ibd[['chr','EAF','beta','pval','n']]
ibd.to_csv(wd+'data/'+'EUR.IBD.gwas_info03_filtered.assoc.tsv.gz',compression='gzip',sep='\t')

## Clean up CD (Crohn's disease) data
cd = pd.read_csv(wd+'data/'+'EUR.CD.gwas_info03_filtered.assoc.gz',compression='gzip',delim_whitespace=True)
n_cas = int([x for x in cd.columns.values if 'FRQ_A' in x][0].split('_')[2])
n_con = int([x for x in cd.columns.values if 'FRQ_U' in x][0].split('_')[2])
cd = cd.rename(columns={f'FRQ_U_{n_con}':'EAF','CHR':'chr','P':'pval'})
cd['beta'] = np.log10(cd.OR)
cd['n'] = n_cas + n_con
cd = cd[['chr','EAF','beta','pval','n']]
cd.to_csv(wd+'data/'+'EUR.CD.gwas_info03_filtered.assoc.tsv.gz',compression='gzip',sep='\t')

## Clean up UC (ulcerative colitis) data
uc = pd.read_csv(wd+'data/'+'EUR.UC.gwas_info03_filtered.assoc.gz',compression='gzip',delim_whitespace=True)
n_cas = int([x for x in uc.columns.values if 'FRQ_A' in x][0].split('_')[2])
n_con = int([x for x in uc.columns.values if 'FRQ_U' in x][0].split('_')[2])
uc = uc.rename(columns={f'FRQ_U_{n_con}':'EAF','CHR':'chr','P':'pval'})
uc['beta'] = np.log10(uc.OR)
uc['n'] = n_cas + n_con
uc = uc[['chr','EAF','beta','pval','n']]
uc.to_csv(wd+'data/'+'EUR.UC.gwas_info03_filtered.assoc.tsv.gz',compression='gzip',sep='\t')