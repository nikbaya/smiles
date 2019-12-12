#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 07:30:10 2019

Plot results of LD of simulated smiles

@author: nbaya
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

smiles_wd = '/Users/nbaya/Documents/lab/smiles/'

phen_dict = {'50_irnt':'Standing height',
             '21001_irnt':'BMI'}

phen = '50_irnt'
n_top_loci = 20
ld_window = 300e3 #1e6
block_mhc=True
random_betas = False

## Prepare data
df = pd.read_csv(smiles_wd+'data/'+phen+f'.gwas.corr.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb{".block_mhc" if block_mhc else ""}{".random_betas" if random_betas else ""}.tsv.gz',
                 compression='gzip',
                 sep='\t')

df.loc[df.beta>0,'sim_raf'] = df.loc[df.beta>0,'alt_af']
df.loc[df.beta<0,'sim_raf'] = 1-df.loc[df.beta<0,'alt_af']
df['sim_rbeta'] = abs(df.beta)

if 'chr' not in df.columns.values or 'pos' not in df.columns.values:
    df['chr'] = df.locus.str.split(':',n=2,expand=True).iloc[:,0].astype(int)
    df['pos'] = df.locus.str.split(':',n=2,expand=True).iloc[:,1].astype(int)
    
if 'sim_varexp' not in df.columns.values:
    df['sim_varexp'] = 2*df['minor_AF']*(1-df['minor_AF'])*df.beta**2

# get locus idx of sentinel variants
if 'locus_idx' not in df.columns.values:
    
    ss_tmp = df.copy() #dataframe from which top loci will be removed
    ss_tmp = ss_tmp.rename(columns={'locus.contig':'chr','locus.position':'pos',
                                    'real_pval':'pval'})
    top_variants = []
    ct = 0
    ss_keep = ss_tmp[[False]*len(ss_tmp)]
    i = 0
    print(f'filtering to top {n_top_loci} loci with window of {ld_window/1e3} kb')
    while len(ss_tmp)>0 and i<n_top_loci:
        ch, pos = ss_tmp[ss_tmp.pval==ss_tmp.pval.min()][['chr','pos']].values[0]
        top_variants.append(f'{ch}:{pos}')
        ss_w_top_locus = ss_tmp[(ss_tmp['chr'].astype(str)==str(ch))&(ss_tmp['pos']>=pos-ld_window/2)&(ss_tmp['pos']<=pos+ld_window/2)].copy() #extract rows around most significant hit
        ss_w_top_locus['loci_rank'] = i
        print(f'locus {i} ({ch}:{pos}) n variants: {len(ss_w_top_locus)}')
        ct += len(ss_w_top_locus)
        ss_keep = ss_keep.append(ss_w_top_locus)
        ss_tmp = ss_tmp[~((ss_tmp['chr']==ch)&(ss_tmp['pos']>=pos-ld_window/2)&(ss_tmp['pos']<=pos+ld_window/2))].copy() #keep rows not around most significant hit
        i += 1
        if len(ss_tmp)==0 and i<n_top_loci:
            print(f'\n############\nWARNING: Ran out of SNPs after reading through top {i} loci\n############')
    for locus_idx, locus in enumerate(top_variants):
        ch, pos = locus.split(':')
        df.loc[(df.chr==int(ch))&(df.pos==int(pos)),'locus_idx'] = locus_idx
    df.to_csv(smiles_wd+'data/'+phen+f'.gwas.corr.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb{".block_mhc" if block_mhc else ""}{".random_betas" if random_betas else ""}.tsv.gz',
                 compression='gzip',
                 sep='\t',
                 index=False)
    

        
## Plot smiles

r2_threshold_ls = [0.5]#[0,0.3,0.6]

#prop_raf_gt_sentinel_ls = [None]*len(r2_threshold_ls)

for idx, r2_threshold in enumerate(r2_threshold_ls):
    
    for locus_idx in range(1,5):
        df_locus = df[~df[f'r_locus{locus_idx}'].isna()]
        sentinel_chr, sentinel_pos = df_locus[df_locus.locus_idx==locus_idx][['chr','pos']].values.flatten()
        df_locus[f'r2_locus{locus_idx}']  = df_locus[f'r_locus{locus_idx}']**2
        df_locus = df_locus.loc[df_locus[f'r2_locus{locus_idx}']>r2_threshold] # only keep the variants with r2 to the sentinel greater than r2_threshold
        df_locus = df_locus.loc[~(df_locus.locus_idx==locus_idx)]
        
        for gwas in ['real','sim']:
            fig,ax=plt.subplots(figsize=(6,4))
            ax.scatter(x=df_locus[f'{gwas}_raf'], y=df_locus[f'{gwas}_rbeta'],s=50,c=df_locus[f'r2_locus{locus_idx}'],
                       marker='.',cmap='viridis',vmin=r2_threshold,vmax=1)
            ax.plot(df.loc[df.locus_idx==locus_idx][f'{gwas}_raf'],df.loc[df.locus_idx==locus_idx][f'{gwas}_rbeta'],'rx',ms=10)
            plt.xlim([0,1])
            plt.ylim([0.95*df.loc[(~df[f'r_locus{locus_idx}'].isna())&(df[f'r_locus{locus_idx}']**2>r2_threshold)][['real_rbeta','sim_rbeta']].min().min(),
                      1.05*df.loc[(~df[f'r_locus{locus_idx}'].isna())&(df[f'r_locus{locus_idx}']**2>r2_threshold)][['real_rbeta','sim_rbeta']].max().max()])
            plt.title(f'Risk AF vs. Effect Size ({gwas}) \nlocus {locus_idx+1}: {sentinel_chr}:{sentinel_pos} , r2>{r2_threshold}')
                    
        
        for gwas in ['real','sim']:
            fig,ax=plt.subplots(figsize=(6,4))
            ax.scatter(x=df_locus[f'{gwas}_raf'], y=df_locus[f'{gwas}_varexp'],s=50,c=df_locus[f'r2_locus{locus_idx}'],
                       marker='.',cmap='viridis',vmin=r2_threshold,vmax=1)
            ax.plot(df.loc[df.locus_idx==locus_idx][f'{gwas}_raf'],df.loc[df.locus_idx==locus_idx][f'{gwas}_varexp'],'rx',ms=10)
            plt.xlim([0,1])
            plt.ylim([0.95*df.loc[(~df[f'r_locus{locus_idx}'].isna())&(df[f'r_locus{locus_idx}']**2>r2_threshold)][['real_varexp','sim_varexp']].min().min(),
                      1.05*df.loc[(~df[f'r_locus{locus_idx}'].isna())&(df[f'r_locus{locus_idx}']**2>r2_threshold)][['real_varexp','sim_varexp']].max().max()])
            plt.title(f'Risk AF vs. Variance Explained ({gwas}) \nlocus {locus_idx+1}: {sentinel_chr}:{sentinel_pos} , r2>{r2_threshold}')
    
    
#    prop_raf_gt_sentinel = [] #proportion of variants in a given locus with raf > sentinel variant raf
#    for locus_idx in range(int(df.locus_idx.max())+1):
#        df_locus = df[~df[f'r_locus{locus_idx}'].isna()]
#        sentinel_chr, sentinel_pos = df_locus[df_locus.locus_idx==locus_idx][['chr','pos']].values.flatten()
#        df_locus[f'r2_locus{locus_idx}']  = df_locus[f'r_locus{locus_idx}']**2
#        df_locus = df_locus.loc[df_locus[f'r2_locus{locus_idx}']>r2_threshold] # only keep the variants with r2 to the sentinel greater than r2_threshold
#        df_locus = df_locus.loc[~(df_locus.locus_idx==locus_idx)]
#        if len(df_locus>0):
#            sentinel_raf = df[df.locus_idx==locus_idx].sim_raf.values[0]
#            prop = (df_locus.sim_raf>sentinel_raf).mean()
#            ct = df_locus.shape[0]
#            print(f'locus {locus_idx} proportion ({ct} snps): {prop}')
#        else:
#            print(f'WARNING: locus {locus_idx} has no variants with r2 > {r2_threshold} with the sentinel variant')
#            prop, ct = float('NaN'), float('NaN')
#        prop_raf_gt_sentinel.append([prop, ct])
#    
#    prop_raf_gt_sentinel_ls[idx] = prop_raf_gt_sentinel
    
    
for idx, r2_threshold in enumerate(r2_threshold_ls):
#    fig, ax = plt.subplots(figsize=(6*1.2,4*1.2))
#    ax.hist(prop_raf_gt_sentinel,20)
#    plt.xlim([0,1])
#    plt.xlabel('Prop. of variants with RAF > sentinel RAF')
#    plt.ylabel('Density')
#    plt.title(f'Distribution of proportion of variants in a locus \nwith RAF > sentinel RAF (phen: {phen_dict[phen]}, top {int(df.locus_idx.max())+1} loci, r2>{r2_threshold})')
#    plt.savefig(smiles_wd+'smiles/plots/'+phen+f'.prop_raf_gt_sentinel.r2_{r2_threshold}.top{n_top_loci}loci.simulated.png',dpi=300)

    prop_raf_gt_sentinel = prop_raf_gt_sentinel_ls[idx]
    prop_raf_gt_sentinel = np.asarray(prop_raf_gt_sentinel)
    prop_raf_gt_sentinel = prop_raf_gt_sentinel[~np.isnan(prop_raf_gt_sentinel[:,0]),:]
    n, bins, patches = plt.hist(prop_raf_gt_sentinel[:,0].astype(float),bins=np.linspace(0,1,21),alpha=0)
    fig, ax = plt.subplots(figsize=(6*1.2,4*1.2))
    prop_raf_gt_sentinel = np.asarray(prop_raf_gt_sentinel)
    ax.plot(prop_raf_gt_sentinel[:,0], prop_raf_gt_sentinel[:,1],'o')
    ax.bar((bins[:-1] + bins[1:]) / 2, n/max(n)*max(prop_raf_gt_sentinel[:,1]),width=0.05,color='grey',alpha=0.1)
    plt.xlim([-0.05,1.05])
    plt.xlabel('Prop. of variants in locus with RAF > sentinel RAF')
    plt.ylabel('Number of variants in locus')
    plt.title(f'Proportion of variants with RAF > sentinel RAF vs. number of variants\n(phen: {phen_dict[phen]}{", random_betas" if random_betas else ""}, top {len(prop_raf_gt_sentinel_ls[idx])} loci, r2>{r2_threshold})')
    plt.savefig(smiles_wd+'smiles/plots/'+phen+f'.prop_raf_gt_sentinel_vs_n_snps.r2_{r2_threshold}.top{n_top_loci}loci.simulated{".random_betas" if random_betas else ""}.png',dpi=300)
    