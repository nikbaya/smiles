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

phen_dict = {#'50_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'Standing height',
#             '21001_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'BMI',
#             'Mahajan.NatGenet2018b.T2Dbmiadj.European.txt.gz':'T2D'} #this is the edited version of the original data, some unnecessary columns were removed
             '2443.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'Diabetes diagnosed'} #Directly from UKB
#        'pgc.scz.full.2012-04.tsv.gz':'SCZ'} #cases: 36,989, controls: 113,075, total: 150064
#        'EUR.IBD.gwas_info03_filtered.assoc.coding.tsv.gz':'IBD',#EUR IBD from transethnic ancestry meta-analysis
#        'EUR.CD.gwas_info03_filtered.assoc.coding.tsv.gz':'CD', #EUR CD from transethnic ancestry meta-analysis
#        'EUR.UC.gwas_info03_filtered.assoc.coding.tsv.gz':'UC'} #EUR UC from transethnic ancestry meta-analysis

highlight_coding = True
ld_clumping = False
ld_window = int(100e3)

for filename, phen in phen_dict.items():
    print(f'Starting phen {phen}, using file {filename}')
    
    ss0 = pd.read_csv(wd+'data/'+filename, sep='\t',compression='gzip')
    
    if 'chr' not in ss0.columns.values:
        ss0['chr'] = ss0.variant.str.split(':',n=1,expand=True).iloc[:,0]
    ss0['chr'] = ss0['chr'].astype(str)
    
    if 'pos' not in ss0.columns.values:
        ss0['pos'] = ss0.variant.str.split(':',n=2,expand=True).iloc[:,1]
    ss0['pos'] = ss0['pos'].astype(int)
    
    if ld_clumping:
        if f'ld_index_{ld_window}' not in ss0.columns.values:
            ss0[f'ld_index_{ld_window}'] = [f'{entry.chr}-{int(entry.pos/ld_window)}' for id,entry in ss0.iterrows()]
            ss0.to_csv(wd+'data/'+filename, sep='\t',compression='gzip')
        ss0 = ss0.loc[ss0.groupby(f'ld_index_{ld_window}')['pval'].idxmin()]
            
            
    for pval_threshold in [1e-5,5e-8,1e-8]:
        
        ss = ss0[ss0.pval < pval_threshold].reset_index()
        if 'n_complete_samples' in ss.columns.values:
            print('renaming field n_complete_samples as n')
            ss = ss.rename(columns={'n_complete_samples':'n'})
        n = int(ss.n.mean())
        if ss.n.std() != 0:
            print('WARNING: Number of samples varies across SNPs')
        
        ss.loc[ss.index,'rbeta'] = np.abs(ss['beta']) #beta is transformed to risk (or trait-increasing) allele effect size
        
        if 'EAF' not in ss.columns.values:
            ss['ref'] = ss.variant.str.split(':',n=3,expand=True).iloc[:,2]
            ss.loc[ss.minor_allele!=ss.ref,'alt_af'] = ss.loc[ss.minor_allele!=ss.ref,'minor_AF']
            ss.loc[ss.minor_allele==ss.ref,'alt_af'] = 1-ss.loc[ss.minor_allele==ss.ref,'minor_AF']
            ss.loc[ss.beta>0,'raf'] = ss.loc[ss.beta>0,'alt_af']
            ss.loc[ss.beta<0,'raf'] = 1-ss.loc[ss.beta<0,'alt_af']
        else:
            ss = ss.rename(columns={'EAF':'raf'}) # where raf is "risk allele frequency"
        for maf in [0, 0.01]:
            if maf != 0:
                ss = ss[(ss.raf>maf)&(ss.raf<1-maf)]
            
            file_prefix = f'{phen}.pval_{pval_threshold}'+(f'.maf_{maf}' if maf!=0 else '')
            
            #effect size beta
            fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))
            if 'coding' in ss.columns.values and highlight_coding:
                ax.plot(ss[~ss.coding].raf, ss[~ss.coding].rbeta,'.',ms=2,c='#1f77b4') #plot noncoding variants
                ax.plot(ss[ss.coding].raf, ss[ss.coding].rbeta, 'o',markerfacecolor='none', markeredgewidth=0.75) #plot coding variants
                plt.legend(['non-coding','coding'])
            else:
                ax.plot(ss.raf, ss.rbeta,'.',ms=2,alpha=1)
            plt.xlabel('Risk allele frequency')
            plt.ylabel('Estimated effect size')
            plt.title(f'AF vs. Effect Size\nphen: {phen}, n: {n}, pval threshold: {pval_threshold}'
                      +(f', maf>{maf}' if maf!=0 else '')+(f', ld block: {int(ld_window/1000)}kb' if ld_clumping else ''))
            fig=plt.gcf()
            suffix = ''
            if 'coding' in ss.columns.values  and highlight_coding:
                suffix = '.highlightcoding'
            if ld_clumping:
                suffix += f'.ldclumping{ld_window}'
            fig.savefig(wd+f'smiles/plots/{file_prefix}.raf_effectsize{suffix}.png',dpi=600)
            plt.close()
                
            #variance explained
            fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))    
            if 'coding' in ss.columns.values  and highlight_coding:
                raf = ss[~ss.coding]['raf']
                rbeta = ss[~ss.coding]['rbeta']
                raf_c = ss[ss.coding]['raf'] #raf for coding variants
                rbeta_c = ss[ss.coding]['rbeta'] #raf for coding variants
                ax.plot(raf, 2*raf*(1-raf)*rbeta**2,'.',ms=2,c='#1f77b4')
                ax.plot(raf_c, 2*raf_c*(1-raf_c)*rbeta_c**2, 'o',markerfacecolor='none', markeredgewidth=0.75,c='#1f77b4')
                plt.legend(['non-coding','coding'])
            else:
                ax.plot(ss.raf, 2*ss.raf*(1-ss.raf)*ss.rbeta**2,'.',ms=2)
            plt.xlabel('Risk allele frequency')
            plt.ylabel('Variance explained')
            plt.title(f'AF vs. Variance Explained\nphen: {phen}, n: {n}, pval threshold: {pval_threshold}'
                      +(f', maf>{maf}' if maf!=0 else '')+(f', ld block: {int(ld_window/1000)}kb' if ld_clumping else ''))
            suffix = ''
            if 'coding' in ss.columns.values  and highlight_coding:
                suffix = '.highlightcoding'
            if ld_clumping:
                suffix += f'.ldclumping{ld_window}'
            fig.savefig(wd+f'smiles/plots/{file_prefix}.raf_varianceexplained{suffix}.png',dpi=600)
            plt.close()
    
            #variance explained, colored by chromosome
            fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            for i, ch in enumerate([*range(1,23)]+['X']):
                if 'coding' in ss.columns.values and highlight_coding:
                    raf = ss[(ss.chr==str(ch))&(~ss.coding)]['raf']
                    rbeta = ss[(ss.chr==str(ch))&(~ss.coding)]['rbeta']
                    raf_c = ss[(ss.chr==str(ch))&(ss.coding)]['raf'] #raf for coding variants
                    rbeta_c = ss[(ss.chr==str(ch))&(ss.coding)]['rbeta'] #raf for coding variants
                    ax.plot(raf, 2*raf*(1-raf)*rbeta**2,'.',ms=2,c=colors[i%10])
                    ax.plot(raf_c, 2*raf_c*(1-raf_c)*rbeta_c**2, 'o',markerfacecolor='none', markeredgewidth=0.75,c=colors[i%10])
                else:
                    raf = ss[ss.chr==str(ch)]['raf']
                    rbeta = ss[ss.chr==str(ch)]['rbeta']
                    ax.plot(raf, 2*raf*(1-raf)*rbeta**2,'.',ms=2)
            plt.xlabel('Risk allele frequency')
            plt.ylabel('Variance explained')
            plt.title(f'AF vs. Variance Explained, colored by chromosome\nphen: {phen}, n: {n}, pval threshold: {pval_threshold}'
                      +(f', maf>{maf}' if maf!=0 else '')+(f', ld block: {int(ld_window/1000)}kb' if ld_clumping else ''))
            suffix = ''
            if 'coding' in ss.columns.values  and highlight_coding:
                suffix = '.highlightcoding'
            if ld_clumping:
                suffix += f'.ldclumping{ld_window}'
            fig.savefig(wd+f'smiles/plots/{file_prefix}.raf_varianceexplained.coloredbychr{suffix}.png',dpi=600)
            plt.close()