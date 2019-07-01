#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:22:48 2019

Create "smiles" plots -- Risk allele effect size as a function of risk allele freq

@author: nbaya
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

wd = "/Users/nbaya/Documents/lab/smiles/"

phen_dict = {'50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz':'Standing height',
             '21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz':'BMI'}
filename='21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz'
phen='BMI'
for filename, phen in phen_dict.items():
    
    ss0 = pd.read_csv(wd+filename, sep='\t',compression='gzip')
            
    for pval_threshold in [1e-8,5e-8,1e-5]:
        
        ss = ss0[ss0.pval < pval_threshold].reset_index()
        ss
        n = int(ss.n_complete_samples.mean())
        
        ss.loc[ss.index,'rbeta'] = np.abs(ss['beta'])
        
        ss['raf'] = ss['minor_AF']
        ss.loc[ss.beta<0,'raf'] = 1-ss.loc[ss.beta<0,'raf']
        
        #effect size beta
        fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))
        ax.plot(ss.raf, ss.rbeta,'.',ms=2,alpha=1)
        plt.xlabel('Risk allele frequency')
        plt.ylabel('Estimated effect size')
        plt.title(f'AF vs. Effect Size\nphen: {phen}, n: {n}, pval_threshold: {pval_threshold}')
        fig=plt.gcf()
        fig.savefig(wd+"plots/"+f'{phen}.pval_{pval_threshold}.raf_effectsize.png',dpi=600)
    
            
        #variance explained
        fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))            
        ax.plot(ss.raf, 2*ss.raf*(1-ss.raf)*ss.rbeta**2,'.',ms=2)
        plt.xlabel('Risk allele frequency')
        plt.ylabel('Variance explained')
        plt.title(f'AF vs. Variance Explained\nphen: {phen}, n: {n}, pval_threshold: {pval_threshold}')
        fig.savefig(wd+"plots/"+f'{phen}.pval_{pval_threshold}.raf_varianceexplained.png',dpi=600)

        #variance explained, colored by chromosome
        fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))
        ss['chr'] = ss.variant.str.split(':',n=1,expand=True).iloc[:,0]
        for ch in list(range(22))+['X']:
            raf = ss[ss.chr==str(ch)]['raf']
            rbeta = ss[ss.chr==str(ch)]['rbeta']
            ax.plot(raf, 2*raf*(1-raf)*rbeta**2,'.',ms=2,alpha=1)
        plt.xlabel('Risk allele frequency')
        plt.ylabel('Variance explained')
        plt.title(f'AF vs. Variance Explained, colored by chromosome\nphen: {phen}, n: {n}, pval_threshold: {pval_threshold}')
        fig.savefig(wd+"plots/"+f'{phen}.pval_{pval_threshold}.raf_varianceexplained.coloredbychr.png',dpi=600)
