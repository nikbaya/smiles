#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 08:05:03 2020

Script for clumping human GWAS sumstats
Originally part of plot_smiles.py

@author: nbaya
"""

import argparse
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
#import multiprocessing as mp
#from functools import partial
from time import time
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

if os.path.isdir('/Users/nbaya/'): # if working locally
    smiles_dir = "/Users/nbaya/Documents/lab/smiles"
else:
    smiles_dir = "/stanley/genetics/users/nbaya/smiles"
    
phen_fname_dict = {
             'standing_height':['50_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz','Standing height'], # UKB
             'bmi':['21001_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz','BMI'], # UKB
             't2d_bmiadj':['Mahajan.NatGenet2018b.T2Dbmiadj.European.coding.tsv.gz','T2D_bmiadj'], #this is the edited version of the original BMI-adjusted data, some unnecessary columns were removed
             'ibd':['EUR.IBD.gwas_info03_filtered.assoc.coding.tsv.gz','IBD'], #,#EUR IBD from transethnic ancestry meta-analysis
             'cd':['EUR.CD.gwas_info03_filtered.assoc.coding.tsv.gz','CD'], #EUR CD from transethnic ancestry meta-analysis
             'uc':['EUR.UC.gwas_info03_filtered.assoc.coding.tsv.gz','UC'], #EUR UC from transethnic ancestry meta-analysis
             'scz':['daner_PGC_SCZ43_mds9.coding.tsv.gz', 'SCZ'], #PGC data for SCZ (NOTE: The SNP effects were flipped when converting from odds ratio because daner files use odds ratios based on A1, presumably the ref allele, unlike UKB results which have betas based on the alt allele)
             'ad':['AD_sumstats_Jansenetal_2019sept.coding.tsv.gz','AD'], #Alzheimer's disease meta-analysis
             'breast_cancer':['breastcancer.michailidou2017.b37.cleaned.coding.tsv.gz','Breast cancer'],
             'cad':['UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.cleaned.coding.tsv.gz','CAD'],
             'ldl':['30780_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz','LDL'],
             'hdl':['30760_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz','HDL'],
             'wbc_count':['30000_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz','WBC count'],
             'rbc_count':['30010_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz', 'RBC count'],
             'urate':['30880_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz','Urate'],
             'systolic_bp':['4080_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz','Systolic BP'],
             'diastolic_bp':['4079_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz','Diastolic BP'],
             'triglycerides':['30870_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz','Triglycerides'],
             }


def get_fname_suffix(mixcompdist, betahat, pointmass):
    suffix= f'{mixcompdist}.{f"{betahat}." if betahat!="beta" else ""}{"" if pointmass else "no_pointmass."}' 
    return suffix

def get_genmap():
    r'''
    Gets genetic map and splits into list of dataframes for faster searching later.
    '''
    genmap = pd.read_csv(f'{smiles_dir}/data/genetic_map_combined_b37.txt.gz',
                         sep=' ', names=['chr','pos','rate','cm'], compression='gzip')
    
    genmap_chr_list = [genmap[genmap.chr==chr].sort_values(by='pos') for chr in range(1,23)]

    return genmap_chr_list
    
def impute_cm(genmap_chr, pos):
    r'''
    Finds genetic distance in cM based on base pair position `pos` or extrapolates 
    genetic distance in cM on a given chromosome `chr`, using a recombination 
    map DataFrame `genmap_chr`
    '''
    if pos in genmap_chr.pos.values:
        return genmap_chr[genmap_chr.pos==pos].cm.values[0]
    else:
        a_pos, a_cm = genmap_chr[genmap_chr.pos<pos].tail(1)[['pos','cm']].values[0] # throws an IndexError if pos < min(genmap_chr.pos)            
        b_pos, b_cm = genmap_chr[genmap_chr.pos>pos].head(1)[['pos','cm']].values[0] # throws an IndexError if pos > max(genmap_chr.pos)
        cm = a_cm + (pos-a_pos)*(b_cm-a_cm)/(b_pos-a_pos)
        return cm

    
def impute_pos(genmap_chr, cm):
    r'''
    Extrapolates position in base pairs based on the genetic distance `cm` in cM, 
    using a single-chromosome recombination map DataFrame `genmap_chr`.
    '''
    if cm in genmap_chr.cm.values:
        return genmap_chr[genmap_chr.cm==cm].pos.values[0]
    else:
        try:
            a_pos, a_cm = genmap_chr[genmap_chr.cm<cm].tail(1)[['pos','cm']].values[0] # throws an IndexError if cm < min(genmap_chr.cm)
        except IndexError:
            pos = 1
            return pos
        try:
            b_pos, b_cm = genmap_chr[genmap_chr.cm>cm].head(1)[['pos','cm']].values[0] # throws an IndexError if pos is > max(genmap_chr.pos)
        except IndexError: 
            pos = np.Inf
            return pos        
        pos = a_pos + (b_pos-a_pos)*(cm-a_cm)/(b_cm-a_cm)
        return pos
    


def pre_clump_qc(ss0, phen, filter_by_varexp, use_ash, block_mhc, mixcompdist, 
                 betahat, pointmass, maf=0.01):
    r'''
    Pre-clumping QC
    '''
    assert 'se' in ss0.columns.values
    if 'chr' not in ss0.columns.values:
        ss0['chr'] = ss0.variant.str.split(':',n=1,expand=True).iloc[:,0]
    else:
        ss0['chr'] = ss0['chr'].astype(str)
    if not ('A1' in ss0.columns.values)&('A2' in ss0.columns.values):
        if 'variant' in ss0.columns.values:
            ss0['A1'] = ss0.variant.str.split(':',expand=True)[2]
            ss0['A2'] = ss0.variant.str.split(':',expand=True)[3]
    assert ('A1' in ss0.columns.values)&('A2' in ss0.columns.values)
    print(f'Filtering to autosomes...(variant ct: {ss0.shape[0]})')
    ss0 = ss0[ss0.chr.isin([str(x) for x in range(1,23)])] # only need to go up to chr 22 because genetic map only covers 22 chr
    print(f'Filtered to autosomes...(variant ct: {ss0.shape[0]})')
    ss0['chr'] = ss0['chr'].astype(int) # previous filtering step is necessary for converting to int for cases of chr='X' or 'MT'
        
    if 'pos' not in ss0.columns.values:
        ss0['pos'] = ss0.variant.str.split(':',n=2,expand=True).iloc[:,1]
    ss0['pos'] = ss0['pos'].astype(int)
        
    if 'n_complete_samples' in ss0.columns.values:
        print('renaming field "n_complete_samples" to "n"')
        ss0 = ss0.rename(columns={'n_complete_samples':'n'})
    sample_ct_col_ls = list({'n_complete_samples','n','Neff'}.intersection(ss0.columns.values)) 
    if len(sample_ct_col_ls)>0: 
        sample_ct_col = sample_ct_col_ls[0] #take first element of list of intersection in the very unlikely case that there are multiple sample ct fields in the intersection
        if ss0[sample_ct_col].std() != 0:
            print('WARNING: Number of samples (using column "{sample_ct_col}") varies across SNPs, these are likely meta-analyzed sumstats.')
    else: 
        print('WARNING: phen {phen} is missing a column indicating sample count.')
        
    if 'EAF' in ss0.columns.values:
        ss0 = ss0.rename(columns={'EAF':'eaf'})
    if 'eaf' not in ss0.columns.values:
        if all(x in ss0.columns.values for x in ['variant','minor_AF']):
            ss0['ref'] = ss0.variant.str.split(':',n=3,expand=True).iloc[:,2] #take the first allele listed in the variant ID string
            ss0.loc[ss0.minor_allele!=ss0.ref,'alt_af'] = ss0.loc[ss0.minor_allele!=ss0.ref,'minor_AF']
            ss0.loc[ss0.minor_allele==ss0.ref,'alt_af'] = 1-ss0.loc[ss0.minor_allele==ss0.ref,'minor_AF']
            ss0 = ss0.rename(columns={'alt_af':'eaf'}) # this is true for UKB, which is the only dataset without EAF
        else:
            assert False, 'insufficient information to calculate risk allele frequency'
    
    # remove SNPs with beta=NA
    snp_ct_before = ss0.shape[0]
    ss0 = ss0.dropna(axis=0, subset=['beta'])
    snp_ct_after = ss0.shape[0]
    print(f'SNPs removed with missing beta: {snp_ct_before-snp_ct_after}')
    
    # remove SNPs with EAF=NA
    snp_ct_before = ss0.shape[0]
    ss0 = ss0.dropna(axis=0, subset=['eaf'])
    snp_ct_after = ss0.shape[0]
    print(f'SNPs removed with missing EAF: {snp_ct_before-snp_ct_after}')
          
    if maf!=None:
        snp_ct_before = ss0.shape[0]
        ss0 = ss0[(ss0.eaf>maf)&(ss0.eaf<1-maf)]
        snp_ct_after = ss0.shape[0]
        print(f'# of SNPs removed with MAF<{maf}: {snp_ct_before-snp_ct_after}')
    else:
        # remove invariant SNPs (MAF=0)
        snp_ct_before = ss0.shape[0]
        ss0 = ss0[(ss0.eaf>0) & (ss0.eaf<1)]
        snp_ct_after = ss0.shape[0]
        print(f'SNPs removed with MAF=0: {snp_ct_before-snp_ct_after}')
    
    if 'low_confidence_variant' in ss0.columns.values:
        snp_ct_before = ss0.shape[0]
        ss0 = ss0[~ss0.low_confidence_variant]
        snp_ct_after = ss0.shape[0]
        print(f'Low-confidence SNPs removed: {snp_ct_before-snp_ct_after}')
              
    if phen=='Breast cancer':
        pval_threshold = 1e-5
        snp_ct_before = ss0.shape[0]
        ss0 = ss0[ss0.pval<=pval_threshold]
        print(f'NOTE: Removed {snp_ct_before-ss0.shape[0]} SNPs in "{phen}" with pval > {pval_threshold}')
        
    # note the number of SNPs with beta=0
    print(f'SNPs with beta=0: {ss0[ss0.beta==0].shape[0]}')
    
    ss0.loc[ss0.beta>=0,'raf'] = ss0.loc[ss0.beta>=0,'eaf']
    ss0.loc[ss0.beta<0,'raf'] = 1-ss0.loc[ss0.beta<0,'eaf']

    ss0.loc[ss0.index,'rbeta'] = np.abs(ss0['beta']) #beta is transformed to risk (or trait-increasing) allele effect size
    
    if filter_by_varexp or use_ash: # do if use_ash=True because we may wish to compare var_exp vs var_exp_ash
        ss0['var_exp'] = 2*ss0.raf*(1-ss0.raf)*ss0.rbeta**2
        
    if use_ash:
        ash_fname = f'{smiles_dir}/data/ash'
        ash_fname += f'.{phen}.{mixcompdist}.{"block_mhc." if block_mhc else ""}'
        ash_fname += f'{f"{betahat}." if betahat!= "beta" else ""}{"" if pointmass else "no_pointmass."}'
        ash_fname += f'{f"maf_{maf}." if maf!=None else ""}tsv.gz'
        ash = pd.read_csv(ash_fname, sep='\t',compression='gzip')
        if ('A1' in ss0.columns.values)&('A2' in ss0.columns.values):
            ss0 = ss0.merge(ash, on=['chr','pos', 'A1', 'A2'])
        else:
            print(f'WARNING: Joining ash results without allele columns')
            ss0 = ss0.merge(ash, on=['chr','pos'])
        print(f'SNPs after merge with ash results: {ss0.shape[0]}')
        if betahat=='var_exp':
            ss0['var_exp_ash'] = ss0.PosteriorMean
        else:
            ss0['var_exp_ash'] = 2*ss0.raf*(1-ss0.raf)*ss0.PosteriorMean**2
    
    return ss0

def get_n_hat_mean(phen, pval = 1e-5, ):
    r'''
    `phen`: phenotype
    `pval`: p-val thresh of clumped file to use for estimating n_hat_mean    
    '''
    clumped_fname = f'clumped_gwas.{phen}.ld_wind_cm_{ld_wind_cm}.block_mhc.{f"pval_{pval}." if pval < 1 else ""}tsv.gz'
    clumped = pd.read_csv(f'{smiles_dir}/data/{clumped_fname}',compression='gzip', sep='\t')
#    if use_ash:
#        clumped['pval_ash'] = 2*stats.norm.cdf(-abs(clumped.PosteriorMean/clumped.PosteriorSD)) # 2*clumped.raf*(1-clumped.raf)*clumped.rbeta**2
#    clumped.loc[clumped[f'pval{"_ash" if use_ash else ""}']==0, f'pval{"_ash" if use_ash else ""}'] = 5e-324
    clumped.loc[clumped['pval']==0, 'pval'] = 5e-324
    clumped = clumped[(clumped.raf!=0)&(clumped.raf!=1)]
    clumped['var_exp'] = 2*clumped.raf*(1-clumped.raf)*clumped.rbeta**2
    clumped = clumped[clumped[f'var_exp']!=0] # remove variants with varexp=0 to avoid infinite n_hat
#    clumped['n_hat'] = stats.chi2.isf(q=clumped[f'pval{"_ash" if use_ash else ""}'],df=1)/clumped[f'var_exp{"_ash" if use_ash else ""}']
    clumped['n_hat'] = stats.chi2.isf(q=clumped['pval'],df=1)/clumped[f'var_exp']
    n_hat_mean = clumped.n_hat.mean() # also called n_bar
    
    return n_hat_mean

def get_varexp_threshold(phen, pval_threshold, return_n_hat_mean=False):
    n_hat_mean = get_n_hat_mean(phen=phen)
    varexp_thresh = stats.chi2.isf(q=pval_threshold,df=1)/n_hat_mean
    if return_n_hat_mean:
        return varexp_thresh, n_hat_mean
    else:
        return varexp_thresh

def threshold(phen, ss, pval_threshold, filter_by_varexp, use_ash=False):
    r'''
    Keeps variants with pval<`pval_threshold` or variance-explained greater than 
    a varexp threshold, if `filter_by_varexp`=True.
    '''
    print(f'pval threshold: {pval_threshold}')
    if filter_by_varexp and pval_threshold<1: # only filter when pval_threshold<1
        print(f'filtering by {"ash " if use_ash else ""}variance-explained, converting pval to varexp')
        varexp_thresh, n_hat_mean = get_varexp_threshold(phen=phen, pval_threshold=pval_threshold, return_n_hat_mean=True)
        print(f'estimated n_eff for {phen}: {n_hat_mean}\nvarexp thresh ({phen}): {varexp_thresh}')
        ss = ss[ss[f'var_exp{"_ash" if use_ash else ""}']>varexp_thresh].reset_index(drop=True)
    elif not filter_by_varexp and pval_threshold<1:        
        ss = ss[ss.pval < pval_threshold].reset_index(drop=True) # filter by p-value
    else:
        ss = ss
    print(f'Post-thresholding num of variants: {ss.shape[0]}\n')    
    return ss

def get_blocked_mhc(ss):
    r'''
    Removes all but the most significant variant in the MHC region
    '''
    genes = pd.read_csv(f'{smiles_dir}/data/cytoBand.txt',delim_whitespace=True,header=None,names=['chr','start','stop','region','gene'])
    mhc_region = genes[(genes.chr=='chr6')&(genes.region.str.contains('p21.'))]
    start = min(mhc_region.start)
    stop = max(mhc_region.stop)
    mhc = ss[(ss.chr==6)&(ss.pos>=start)&(ss.pos<=stop)]
    if len(mhc)>0: #if there are variants in MHC
        print(f'\nNumber of variants in MHC: {len(mhc)}')
        non_mhc = ss[~((ss.chr==6)&(ss.pos>=start)&(ss.pos<=stop))]
        ss = non_mhc.append(mhc[mhc.pval==mhc.pval.min()])
    print(f'Post-filter # of variants keeping only one variant in MHC: {ss.shape[0]}')
    
    return ss

def get_clumped_ss(ss, ld_wind_cm, pval_threshold, use_ash):
    r'''
    Performs clumping on sumstats dataframe `ss`. Iterates through most 
    significant variants, removing all but the sentinel (most significant) variant 
    in a given window `ld_wind_cm` defined in cM.
    '''
    start = time()
    ss_noncoding = ss[~ss.coding]
    ss_coding = ss[ss.coding]
    ss_final = None
    print(f'\nGetting loci for LD window={ld_wind_cm}cM ...')
    print(f'Pre-clump number of variants (pval={pval_threshold}, use_ash={use_ash}): {ss.shape[0]}')
    for ss_tmp in [ss_noncoding, ss_coding]:
#        pool = mp.Pool(2) #mp.cpu_count()-1)
#        func = partial(clump_chrom, ss_tmp, ld_wind_cm, use_ash)
#        ss_keep_list = pool.map(func, range(1,23))
#        pool.close()
        ss_keep_list = [clump_chrom(ss_tmp=ss_tmp,ld_wind_cm=ld_wind_cm,use_ash=use_ash,chrom=chrom) for chrom in range(1,23)]
        ss_keep = pd.concat(ss_keep_list)
        ss_final = ss_keep if ss_final is None else ss_final.append(ss_keep)
        
    print(f'Time to clump: {round((time()-start)/60,3)} min')
    print(f'Post-clump number of variants (LD window={ld_wind_cm}cM): {ss_final.shape[0]}\n')

    return ss_final

def clump_chrom(ss_tmp, ld_wind_cm, use_ash, chrom):
    r'''
    Clumps variants on a single chromosome `chrom` for a given sumstats 
    dataframe `ss_tmp`
    '''
#    print(f'chr {chrom}')
    ss_keep = ss_tmp[[False]*len(ss_tmp)]
    genmap_chr = genmap_chr_list[chrom-1]
    if ss_tmp.chr.dtype!=str:
        ss_tmp['chr'] = ss_tmp.chr.astype(str)
        chrom = str(chrom)
    ss_tmp = ss_tmp[ss_tmp.chr==chrom]
    while len(ss_tmp)>0:
#        idx = ss_tmp.var_exp_ash==ss_tmp.var_exp_ash.max() if use_ash else ss_tmp.pval==ss_tmp.pval.min()
#        print(ss_tmp[idx][['chr','pos']])
#        chr, pos = ss_tmp[idx][['chr','pos']].values[0] # sentinel variant is most significant non-clumped variant
        idx = ss_tmp.var_exp_ash.idxmax() if use_ash else ss_tmp.pval.idxmin()
        chr, pos = ss_tmp.loc[idx,['chr','pos']]
#        print(pos)
        if pos < genmap_chr.pos.min():
            a_pos = 1 # start of window around sentinel
            b_pos = impute_pos(genmap_chr=genmap_chr, cm=ld_wind_cm/2) # end of window around sentinel
        elif pos > genmap_chr.pos.max():
            a_pos = impute_pos(genmap_chr=genmap_chr, cm=genmap_chr.cm.max()-ld_wind_cm/2) # start of window around sentinel
            b_pos = np.Inf # end of window around sentinel
        else:
            cm = impute_cm(genmap_chr=genmap_chr, pos=pos) # genetic distance of sentinel in cM
            a_pos = impute_pos(genmap_chr=genmap_chr, cm=cm-ld_wind_cm/2) # start of window around sentinel
            b_pos = impute_pos(genmap_chr=genmap_chr, cm=cm+ld_wind_cm/2) # end of window around sentinel
#        print(f'[{a_pos}, {b_pos}]')
        ss_w_top_locus = ss_tmp.loc[idx,:].copy() #extract sentinel variant, NOTE: This step must come before extracting rows in window around sentinel
        ss_tmp = ss_tmp[~((ss_tmp['chr']==chr)&(ss_tmp['pos']>=a_pos)&(ss_tmp['pos']<=b_pos))].copy() #extract rows not in window around current sentinel variant
        ss_keep = ss_keep.append(ss_w_top_locus)
    
    return ss_keep

def plot_varexp_original_vs_ash(phen, ss, block_mhc, mixcompdist, betahat, pointmass, 
                                maf, filter_by_varexp=True, clumped=False, pval_threshold=None, 
                                logscale=True, color_by='pval'):
    assert all([f in ss.columns.values for f in ['var_exp','var_exp_ash']]), "fields missing from sumstats file"
    print(f'Number of variants: {ss.shape[0]}')
#    if logscale:
#        ss = ss[ss.var_exp != 0]
#        ss = ss[ss.var_exp_ash != 0]
#        print(f'Number of variants with var_exp or var_exp_ash != 0: {ss.shape[0]}')
#    plt.plot(ss['var_exp'],ss['var_exp_ash'],'.',alpha=0.05)
    if color_by=='pval':
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        pval_color_mapping = lambda x: 'lightgrey' if x.pval > 1e-5 else (
                colors[0] if x.pval > 1e-6 else (
                        colors[1] if x.pval > 1e-7 else (
                                colors[2] if x.pval > 5e-8 else (
                                        colors[3] if x.pval > 1e-8 else colors[4]))))
    elif color_by=='maf' and not 'maf' in ss.columns:
        calc_maf = lambda x: min(x.eaf, 1-x.eaf)
        ss['maf'] = ss.apply(calc_maf, axis=1)
    plt.figure(figsize=(6*1.5,4*1.5))
    plt.scatter(ss['var_exp'],ss['var_exp_ash'], 
                c=ss[color_by] if color_by!='pval' else ss.apply(pval_color_mapping, axis=1), 
                cmap='cool' if color_by!='pval' else None,  # other options: bwr, Spectral, coolwarm
                alpha=0.1)
    if color_by!='pval':
        plt.colorbar()
    else:
        legend_elements = [Line2D([0], [0], marker='.', color='w', markerfacecolor='lightgrey', label='>1e-5', markersize=15),
                           Line2D([0], [0], marker='.', color='w', markerfacecolor=colors[0], label='<1e-5', markersize=15),
                           Line2D([0], [0], marker='.', color='w', markerfacecolor=colors[1], label='<1e-6', markersize=15),
                           Line2D([0], [0], marker='.', color='w', markerfacecolor=colors[2], label='<1e-7', markersize=15),
                           Line2D([0], [0], marker='.', color='w', markerfacecolor=colors[3], label='<5e-8', markersize=15),
                           Line2D([0], [0], marker='.', color='w', markerfacecolor=colors[4], label='<1e-8', markersize=15),
                           ]
        plt.legend(handles=legend_elements, loc=4)
    plt.plot(*[[0,ss[['var_exp','var_exp_ash']].max().min()]]*2,'k--')
    if logscale:
        plt.xscale('symlog', linthreshx=1e-10)
        plt.yscale('symlog', linthreshy=1e-10)
        plt.xlim([-1e-11,ss.var_exp.max()*2])
        plt.ylim([-1e-11,ss.var_exp_ash.max()*2])
    else:
        plt.xlim([-ss.var_exp.max()*0.01,ss.var_exp.max()*1.01])
        plt.ylim([-ss.var_exp_ash.max()*0.01,ss.var_exp_ash.max()*1.01])
    plt.xlabel('variance explained from original GWAS')
    plt.ylabel('variance explained from ash')
    plt.title(f'{phen_fname_dict[phen][1]}{f", pval={pval_threshold}" if pval_threshold not in [None, 1] else ""}'+
              f'{" (var-exp thresh.)" if filter_by_varexp&(pval_threshold not in [None, 1]) else ""}'+
              f'{", clumped" if clumped else ""}'+
              f'{", maf>{maf}" if maf!=None else ""}'+
              f', color={color_by}'+
              f'\n{mixcompdist}, betahat={betahat}, pointmass={pointmass}, block_mhc={block_mhc} ({ss.shape[0]} SNPs)')
    plt.tight_layout()
    plot_fname = f'{smiles_dir}/plots/varexp_comparison.'
    plot_fname += "logscale." if logscale else ""
    plot_fname += phen
    plot_fname += f'{".block_mhc" if block_mhc else ""}.'
    plot_fname += f'{get_fname_suffix(mixcompdist, betahat, pointmass)}'
    plot_fname += f"maf_{maf}." if maf!=None else ""
    plot_fname += f'{"clumped." if clumped else ""}{f"pval_{pval_threshold}." if pval_threshold!=None else ""}color_by_{color_by.replace(" ","_")}.png'
    plt.savefig(plot_fname, dpi=300)
    
def plot_smiles_original_vs_ash(phen, ss, block_mhc, mixcompdist, betahat, pointmass, maf,
                                pval_threshold=None, clumped=False, highlight_coding=False):
    assert all([f in ss.columns.values for f in ['var_exp','var_exp_ash']]), "fields missing from sumstats file"
    plt.figure(figsize=(6*1.5,4*1.5))
    print(f'Number of variants: {ss.shape[0]}')
#    ss = ss[ss.var_exp != 0]
#    ss = ss[ss.var_exp_ash != 0]
    print(f'Number of variants with var_exp or var_exp_ash != 0: {ss.shape[0]}')
    if highlight_coding:
        assert False, 'Plotting with `highlight_coding`=True not implemented'
    else:
        plt.plot(ss['raf'],ss['var_exp'],'.',alpha=0.1)
        plt.plot(ss['raf'],ss['var_exp_ash'],'.',alpha=0.1)
        plt.legend(['original','ash'])
    ticks = list(map(lambda x: float('%s' % float('%.2g' % x)), np.logspace(np.log10(ss.var_exp_ash.min()), np.log10(ss.var_exp_ash.max()), 3)))
    plt.yticks(ticks)
#    plt.ylim([ss.var_exp_ash.min(), ss.var_exp_ash.max()*1.01])
    plt.xlabel('raf')
    plt.ylabel('variance explained')
    plt.yscale('symlog', linthreshy=1e-8)
    plt.title(f'{phen_fname_dict[phen][1]} (variants: {ss.shape[0]})\n{mixcompdist}, betahat={betahat}, pointmass={pointmass}, block_mhc={block_mhc}, clumped={clumped}')
    plt.tight_layout()
    plot_fname = f'{smiles_dir}/plots/{"clumped_" if clumped else ""}'
    plot_fname += f'smiles_original_vs_ash.{phen}'
    plot_fname += f'{mixcompdist}.{"block_mhc." if block_mhc else ""}{f"pval_{pval_threshold}." if pval_threshold!=None else ""}'
    plot_fname += f'{get_fname_suffix(mixcompdist, betahat, pointmass)}png'
    plt.savefig(plot_fname, dpi=300)
    plt.close()
    
def manhattan_plot(ss, phen, use_varexp, maf, use_ash=False, subtitle=None, show=True,
                   chroms = range(1,23)):
    if use_ash: use_varexp = True
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    if ss.chr.dtype != str:
        ss['chr'] = ss.chr.astype(int).astype(str)
    plt.figure(figsize=(12,6))
    if not use_varexp:
        ss['nlog10p'] = -np.log10(ss.pval)
    prev_pos = 0

    for chrom in chroms:
        ss_chr = ss[ss.chr==str(chrom)]
#        print(ss_chr.shape[0])
        min_pos = genmap_chr_list[chrom-1].pos.min() #min(ss_chr.pos.min(), genmap_chr_list[chrom-1].pos.min())
        max_pos = genmap_chr_list[chrom-1].pos.max() #max(ss_chr.pos.max(), genmap_chr_list[chrom-1].pos.max())
        if (ss_chr.pos<min_pos).any():
            print(f'(chr{chrom}) Number of variants with pos<min_pos: {(ss_chr.pos<min_pos).sum()}')
            ss_chr.loc[ss_chr.pos<min_pos,'pos'] = min_pos
        if (ss_chr.pos>max_pos).any():
            print(f'(chr{chrom}) Number of variants with pos>max_pos: {(ss_chr.pos>max_pos).sum()}')
            ss_chr.loc[ss_chr.pos>max_pos,'pos'] = max_pos
        plt.plot(ss_chr.pos-min_pos+prev_pos, 
                 ss_chr.nlog10p if not use_varexp else (ss_chr.var_exp_ash if use_ash else ss_chr.var_exp),
                 '.', c = colors[chrom%2])#, c=[ colors[chrom%2] for chrom in ss.chr])
        prev_pos = prev_pos+max_pos
    min_threshold = None
    max_threshold = None
    for pval in [1e-5, 1e-6, 1e-7, 1e-8]:
        threshold = get_varexp_threshold(phen=phen, pval_threshold=pval) if use_varexp else -np.log10(pval)
        plt.plot([0, prev_pos], [threshold]*2, 'k--')
        min_threshold = threshold if min_threshold==None else min(threshold,min_threshold)
        max_threshold = threshold if max_threshold==None else max(threshold,max_threshold)
    if use_varexp:
        plt.yscale('symlog')
#        ticks = list(map(lambda x: float('%s' % float('%.2g' % x)), np.linspace(ss.var_exp.min(), ss.var_exp.max(), 5)))
        ticks = list(map(lambda x: float('%s' % float('%.2g' % x)), np.logspace(np.log10(ss[ss.var_exp>0].var_exp.min()), np.log10(ss.var_exp.max()), 2)))
#        print(ticks)
        plt.yticks(ticks)
    plt.ylim([min_threshold, max_threshold])
    plt.ylabel('variance explained' if use_varexp else '-log10(p)' )
    plt.title(f'Manhattan plot for {phen_fname_dict[phen][1]} (variants: {ss.shape[0]})'+
              ("\n"+subtitle if subtitle!=None else ""))
    plt.savefig(f'{smiles_dir}/plots/manhattan.{phen}{"."+subtitle.replace(", ",".").replace("=","_").lower() if subtitle is not None else ""}.png',
                dpi=300)
    if not show:
        plt.close()

def main(phen):
    fname, phen_str = phen_fname_dict[phen]
    print(f'Starting phen {phen_str}, using file {fname}')
    print(f'LD clumping: {ld_clumping}')
    print(f'Block MHC: {block_mhc}')
    print(f'Filter by var-exp: {filter_by_varexp}')
    print(f'Use ash results: {use_ash}'+
          f'\nmixcompdist: {mixcompdist}\nash betahat: {betahat}\npointmass: {pointmass}' if use_ash else '')
    if ld_clumping:
#            print(f'LD window: {ld_wind_kb/1e3}mb, {ld_wind_kb}kb')
        print(f'LD window: {ld_wind_cm}cM')
    print(f'MAF: {maf}' if maf!=None else '')
    
    print(f'Loading dataframe for {phen_str}')
    ss0 = pd.read_csv(f'{smiles_dir}/data/{fname}', sep='\t',compression='gzip')
    ss1 = pre_clump_qc(ss0=ss0,
                       phen=phen,
                       filter_by_varexp=filter_by_varexp,
                       use_ash=use_ash,
                       block_mhc=block_mhc,
                       mixcompdist=mixcompdist,
                       betahat=betahat,
                       pointmass=pointmass,
                       maf=maf)
        
    if use_ash and make_plots:
        plot_varexp_original_vs_ash(phen=phen, 
                                    ss=ss1, 
                                    block_mhc=block_mhc, 
                                    mixcompdist=mixcompdist, 
                                    betahat=betahat, 
                                    pointmass=pointmass, 
                                    maf=maf,
                                    filter_by_varexp=filter_by_varexp, 
                                    clumped=False, 
                                    pval_threshold=None, 
                                    logscale=True, 
                                    color_by='pval')
                    
        plot_smiles_original_vs_ash(phen=phen,
                                    ss=ss1,
                                    block_mhc=block_mhc,
                                    mixcompdist=mixcompdist,
                                    betahat=betahat,
                                    pointmass=pointmass,
                                    maf=maf)
        
    ss1 = get_blocked_mhc(ss1) if block_mhc else ss1
    
    if ld_clumping:
        ss_clumped = None
        pre_clumped = False
        
        for pval_threshold in pval_thresholds:
            # TODO: Account for case where ld_clumping=False and pre_clumped=False
            if not pre_clumped:
                ss = threshold(phen=phen,
                               ss=ss1, 
                               pval_threshold=pval_threshold, 
                               filter_by_varexp=filter_by_varexp,
                               use_ash=use_ash)
            
            else:
                ss_clumped = threshold(phen=phen,
                                       ss=ss_clumped, 
                                       pval_threshold=pval_threshold, 
                                       filter_by_varexp=filter_by_varexp,
                                       use_ash=use_ash)
            if ld_clumping and not pre_clumped:
                ss_clumped = get_clumped_ss(ss=ss,
                                            ld_wind_cm=ld_wind_cm,
                                            pval_threshold=pval_threshold,
                                            use_ash=use_ash)
                pre_clumped = True # computational speed up to take advantage of clumped version for less stringent pval thresholds

            if use_ash and make_plots:
                for color_by in ['pval','maf']:
                    plot_varexp_original_vs_ash(phen=phen, 
                                                ss=ss_clumped, 
                                                block_mhc=block_mhc, 
                                                mixcompdist=mixcompdist, 
                                                betahat=betahat, 
                                                pointmass=pointmass, 
                                                maf=maf,
                                                filter_by_varexp=filter_by_varexp, 
                                                clumped=True, 
                                                pval_threshold=pval_threshold, 
                                                logscale=True, 
                                                color_by='pval')
                
                plot_smiles_original_vs_ash(phen=phen, 
                                            ss=ss_clumped, 
                                            block_mhc=block_mhc,
                                            pval_threshold=pval_threshold,
                                            mixcompdist=mixcompdist, 
                                            clumped=ld_clumping,
                                            betahat=betahat,
                                            pointmass=pointmass,
                                            maf=maf)
            
#                if make_plots:
#                    manhattan_plot(ss=ss_clumped, phen=phen, 
#                                   use_varexp=filter_by_varexp, 
#                                   use_ash=use_ash, 
#                                   chroms=range(1,23),
#                                   subtitle=f'clumped, pval_thresh={pval_threshold}, varexp_thresholded={filter_by_varexp}, use_ash={use_ash}')
            
            if save_files:
                out_fname = f'{"clumped_" if ld_clumping else ""}'
                out_fname += "ash" if use_ash else "gwas"
                out_fname += f'.{phen}.'
                out_fname += f'ld_wind_cm_{ld_wind_cm}.' if ld_clumping else ''
                out_fname += f'{"block_mhc." if block_mhc else ""}'
                out_fname += f'pval_{pval_threshold}.{"varexp_thresh." if filter_by_varexp else ""}' if pval_threshold <1 else ''
                out_fname += get_fname_suffix(mixcompdist, betahat, pointmass) if use_ash else ''
                out_fname += f'maf_{maf}.' if maf!=None else ''
                out_fname += 'tsv.gz'
                ss_clumped.to_csv(f'{smiles_dir}/data/'+out_fname,sep='\t',index=False, compression='gzip')


#           for fname, phen in fname_dict.items():
#            for pval_threshold in [1]: #pval_thresholds:
#                phen_str = get_phen_str(phen=phen)
#                out_fname = f'{"clumped_" if ld_clumping else ""}'
#                out_fname += "ash" if use_ash else "gwas"
#                out_fname += f'.{phen_str}.'
#                out_fname += f'ld_wind_cm_{ld_wind_cm}.' if ld_clumping else ''
#                out_fname += f'{"block_mhc." if block_mhc else ""}'
#                out_fname += f'pval_{pval_threshold}.{"varexp_thresh." if filter_by_varexp else ""}' if pval_threshold <1 else ''
#                out_fname += get_fname_suffix(mixcompdist, betahat, pointmass)
#                out_fname += 'tsv.gz'
#                ss_clumped = pd.read_csv(f'{smiles_dir}/data/{out_fname}',sep='\t', compression='gzip')
    
    
if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--phen', type=str, default=None, help='If included, only run this phenotype')
    parser.add_argument('--filter_by_varexp', type=str, default='True', help='whether to filter by variance-explained threshold, converted from a p-value threshold')
    parser.add_argument('--use_ash', type=str, default='True', help='whether to clump ash results instead of GWAS results')
    parser.add_argument('--betahat', type=str, help='what ash is fit on (options: beta, var_exp, abs_beta)')
    parser.add_argument('--mixcompdist', type=str, help='ash-specific parameter for mixture distribution (options: halfnormal, +uniform, normal, uniform)')
    parser.add_argument('--pointmass', type=str, help='whether a point mass at zero was included in the ash prior (options: True, False)')
    args = parser.parse_args()    

    ld_clumping = True #only show the top hit (variant with lowest p-value) in a window of size `ld_wind_kb` kb or `ld_wind_cm` cm
#    ld_wind_kb = 300           #int(300e3) if get_top_loci else int(1000e3) #measured in kb; default for when ld_clumping is true: 100e3; default for when get_top_loci is true: 300e3
    ld_wind_cm = 0.6            # 1 Mb = 1 cM -> 600 Kb = 0.6 cM
    block_mhc = True            # remove all but the sentinel variant in the MHC region
    maf = 0.01
    
    filter_by_varexp = bool(args.filter_by_varexp) # whether to filter by variance-explained threshold, converted from a p-value threshold
    use_ash = bool(args.use_ash) # if True, use ash posterior mean effect sizes to get var_exp and clump on var_exp instead of p-val
    betahat = args.betahat # what betahat ash was fit on. options: beta (default), abs_beta, var_exp
    mixcompdist = args.mixcompdist # ash-specific parameter for mixture distribution
    pointmass= bool(args.pointmass) # whether a point mass at zero was included in the ash prior (default: true)
    
    make_plots=True # whether to make plots
    save_files = True# whether to save summary stats dataframe post filtering
    
    if not filter_by_varexp:
        pval_thresholds = [1e-5]
    elif use_ash:
        pval_thresholds = [1, 1e-5]
    
    pval_thresholds = sorted(pval_thresholds,reverse=True)
    genmap_chr_list = get_genmap()
    
    if args.phen==None:
        for phen, _ in phen_fname_dict.items():
            main(phen=phen)
    else:
        main(phen=args.phen)
    