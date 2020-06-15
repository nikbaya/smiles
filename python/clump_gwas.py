#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 08:05:03 2020

Script for clumping human GWAS sumstats
Originally part of plot_smiles.py

@author: nbaya
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
import multiprocessing as mp
from functools import partial
from time import time

smiles_wd = "/Users/nbaya/Documents/lab/smiles"

fname_dict = {
             '50_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'Standing height', # UKB
             '21001_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'BMI', # UKB
             'Mahajan.NatGenet2018b.T2Dbmiadj.European.coding.tsv.gz':'T2D_bmiadj', #this is the edited version of the original BMI-adjusted data, some unnecessary columns were removed
             'EUR.IBD.gwas_info03_filtered.assoc.coding.tsv.gz':'IBD',#,#EUR IBD from transethnic ancestry meta-analysis
             'EUR.CD.gwas_info03_filtered.assoc.coding.tsv.gz':'CD', #EUR CD from transethnic ancestry meta-analysis
             'EUR.UC.gwas_info03_filtered.assoc.coding.tsv.gz':'UC', #EUR UC from transethnic ancestry meta-analysis
             'daner_PGC_SCZ43_mds9.coding.tsv.gz':'SCZ', #PGC data for SCZ (NOTE: The SNP effects were flipped when converting from odds ratio because daner files use odds ratios based on A1, presumably the ref allele, unlike UKB results which have betas based on the alt allele)
             'AD_sumstats_Jansenetal_2019sept.coding.tsv.gz':'AD', #Alzheimer's disease meta-analysis
             }

def get_genmap():
    r'''
    Gets genetic map and splits into list of dataframes for faster searching later.
    '''
    genmap = pd.read_csv(f'{smiles_wd}/data/genetic_map_combined_b37.txt.gz',
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

def pre_clump_qc(ss0, phen, filter_by_varexp, use_ash):
    r'''
    Pre-clumping QC
    '''
    if 'chr' not in ss0.columns.values:
        ss0['chr'] = ss0.variant.str.split(':',n=1,expand=True).iloc[:,0]
    ss0 = ss0[ss0.chr.isin([str(x) for x in range(1,23)])] # only need to go up to chr 22 because genetic map only covers 22 chr
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
        print('WARNING: phen {phen} is missing number of samples.')
        
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
    print(f'# of SNPs removed with missing beta: {snp_ct_before-snp_ct_after}')
    
    # remove SNPs with EAF=NA
    snp_ct_before = ss0.shape[0]
    ss0 = ss0.dropna(axis=0, subset=['eaf'])
    snp_ct_after = ss0.shape[0]
    print(f'# of SNPs removed with missing EAF: {snp_ct_before-snp_ct_after}')
          
    # remove invariant SNPs (MAF=0)
    snp_ct_before = ss0.shape[0]
    ss0 = ss0[(ss0.eaf>0) & (ss0.eaf<1)]
    snp_ct_after = ss0.shape[0]
    print(f'# of SNPs removed with MAF=0: {snp_ct_before-snp_ct_after}\n')
    
    ss0.loc[ss0.beta>0,'raf'] = ss0.loc[ss0.beta>0,'eaf']
    ss0.loc[ss0.beta<0,'raf'] = 1-ss0.loc[ss0.beta<0,'eaf']

    ss0.loc[ss0.index,'rbeta'] = np.abs(ss0['beta']) #beta is transformed to risk (or trait-increasing) allele effect size
    
    if 'low_confidence_variant' in ss0.columns.values:
        ss0 = ss0[~ss0.low_confidence_variant]
    
    if filter_by_varexp:
        ss0['var_exp'] = 2*ss0.raf*(1-ss0.raf)*ss0.rbeta**2
        
    if use_ash:
        phen_str = phen.replace(' ','_').lower()
        ash = pd.read_csv(f'{smiles_wd}/data/ash.{phen_str}.tsv.gz', sep='\t',compression='gzip')
        ss0 = ss0.merge(ash, on=['chr','pos'])
        ss0['var_exp_ash'] = 2*ss0.PosteriorMean*(1-ss0.raf)*ss0.rbeta**2
        
    return ss0

def threshold(ss, pval_threshold, filter_by_varexp):
    r'''
    Keeps variants with pval<`pval_threshold` or variance-explained greater than 
    a varexp threshold, if `filter_by_varexp`=True.
    '''
    print(f'pval threshold: {pval_threshold}')
    if filter_by_varexp and pval_threshold<1: # only filter when pval_threshold<1
        print('filtering by variance-explained, converting pval to varexp')
        phen_str = phen.replace(' ','_').lower()
        clumped = pd.read_csv(f'{smiles_wd}/data/clumped_gwas.{phen_str}.ld_wind_cm_{ld_wind_cm}.{"block_mhc." if block_mhc else ""}pval_1e-05.tsv.gz',
                              compression='gzip', sep='\t')
        clumped.loc[clumped.pval==0, 'pval'] = 5e-324
        clumped = clumped[(clumped.raf!=0)&(clumped.raf!=1)]
        clumped['var_exp'] = 2*clumped.raf*(1-clumped.raf)*clumped.rbeta**2
        clumped['n_hat'] = stats.chi2.isf(q=clumped.pval,df=1)/clumped.var_exp
        n_hat_mean = clumped.n_hat.mean() # also called n_bar
        varexp_thresh = stats.chi2.isf(q=pval_threshold,df=1)/n_hat_mean
        print(f'estimated n_eff: {n_hat_mean}\nvarexp thresh: {varexp_thresh}')
        ss = ss[ss.var_exp>varexp_thresh].reset_index(drop=True)
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
    genes = pd.read_csv(f'{smiles_wd}/data/cytoBand.txt',delim_whitespace=True,header=None,names=['chr','start','stop','region','gene'])
    mhc_region = genes[(genes.chr=='chr6')&(genes.region.str.contains('p21.'))]
    start = min(mhc_region.start)
    stop = max(mhc_region.stop)
    mhc = ss[(ss.chr==6)&(ss.pos>=start)&(ss.pos<=stop)]
    if len(mhc)>0: #if there are variants in MHC
        print(f'Number of variants in MHC: {len(mhc)}')
        non_mhc = ss[~((ss.chr==6)&(ss.pos>=start)&(ss.pos<=stop))]
        ss = non_mhc.append(mhc[mhc.pval==mhc.pval.min()])
    print(f'Post-filter # of variants keeping only one variant in MHC: {ss.shape[0]}')
    
    return ss

def get_clumped_ss(ss, ld_wind_cm, use_ash):
    r'''
    Performs clumping on sumstats dataframe `ss`. Iterates through most 
    significant variants, removing all but the sentinel (most significant) variant 
    in a given window `ld_wind_cm` defined in cM.
    '''
    start = time()
    ss_noncoding = ss[~ss.coding]
    ss_coding = ss[ss.coding]
    ss_final = None
    print(f'Getting loci for LD window={ld_wind_cm}cM ...\nPre-filter # of variants (pval={pval_threshold}) = {ss.shape[0]}')
    for ss_tmp in [ss_noncoding, ss_coding]:
#        pool = mp.Pool(mp.cpu_count()-1)
#        func = partial(clump_chrom, ss_tmp, ld_wind_cm, use_ash)
#        ss_keep_list = pool.map(func, range(1,23))
#        pool.close()
        ss_keep_list = [clump_chrom(ss_tmp=ss_tmp,ld_wind_cm=ld_wind_cm,use_ash=True,chrom=chrom) for chrom in range(1,23)]
        ss_keep = pd.concat(ss_keep_list)
        ss_final = ss_keep if ss_final is None else ss_final.append(ss_keep, sort=False)
        
    print(f'Time to clump: {round((time()-start)/60,3)} min')
    print(f'Post-filter # of variants (LD window={ld_wind_cm}cM): {ss_final.shape[0]}\n')

    return ss_final
  
def clump_chrom(ss_tmp, ld_wind_cm, use_ash, chrom):
    r'''
    Clumps variants on a single chromosome `chrom` for a given sumstats 
    dataframe `ss_tmp`
    '''
    ss_keep = ss_tmp[[False]*len(ss_tmp)]
    genmap_chr = genmap_chr_list[chrom-1]
    ss_tmp = ss_tmp[ss_tmp.chr==chrom]
    while len(ss_tmp)>0:
        idx = ss_tmp.var_exp_ash==ss_tmp.var_exp_ash.max() if use_ash else ss_tmp.pval==ss_tmp.pval.min()
        chr, pos = ss_tmp[idx][['chr','pos']].values[0] # sentinel variant is most significant non-clumped variant
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
        ss_w_top_locus = ss_tmp[(ss_tmp['chr']==chr)&(ss_tmp['pos']==pos)].copy() #extract sentinel variant, NOTE: This step must come before extracting rows in window around sentinel
        ss_tmp = ss_tmp[~((ss_tmp['chr']==chr)&(ss_tmp['pos']>=a_pos)&(ss_tmp['pos']<=b_pos))].copy() #extract rows not in window around current sentinel variant
        ss_keep = ss_keep.append(ss_w_top_locus, sort=False)
    
    return ss_keep

if __name__=="__main__":
    ld_clumping = True          #only show the top hit (variant with lowest p-value) in a window of size `ld_wind_kb` kb or `ld_wind_cm` cm
#    ld_wind_kb = 300           #int(300e3) if get_top_loci else int(1000e3) #measured in kb; default for when ld_clumping is true: 100e3; default for when get_top_loci is true: 300e3
    ld_wind_cm = 0.6            # 1 Mb = 1 cM -> 600 Kb = 0.6 cM
    block_mhc = True            # remove all but the sentinel variant in the MHC region
    save = True                 # whether to save summary stats dataframe post filtering
    filter_by_varexp = False    # whether to filter by variance-explained threshold, converted from a p-value threshold
    use_ash = True              # if True, use ash posterior mean effect sizes to get var_exp and clump on var_exp instead of p-val

#    pval_thresholds = [1]
    pval_thresholds = [1e-5, 1e-6, 1e-7, 1e-8]
#    pval_thresholds = sorted([1, 1e-5, 1e-6, 1e-7, 5e-8, 1e-8],reverse=True)

    genmap_chr_list = get_genmap()
    
    for fname, phen in fname_dict.items():
        print(f'Starting phen {phen}, using file {fname}')
        print(f'LD clumping: {ld_clumping}')
        print(f'Block MHC: {block_mhc}')
        if ld_clumping:
#            print(f'LD window: {ld_wind_kb/1e3}mb, {ld_wind_kb}kb')
            print(f'LD window: {ld_wind_cm}cM')
        
        ss0 = pd.read_csv(f'{smiles_wd}/data/{fname}', sep='\t',compression='gzip')
        
        ss1 = pre_clump_qc(ss0=ss0,
                           phen=phen,
                           filter_by_varexp=filter_by_varexp,
                           use_ash=use_ash)

        ss1 = get_blocked_mhc(ss1) if block_mhc else ss1
        
        ## Use to check bug 
##        ss1 = ss1.sort_values(by=['pval','chr','pos'])
#        ss2 = ss1.drop_duplicates(subset='pval',keep='first')
#        ss2 = ss2[~ss2.coding]
        
        ss_clumped = None
        pre_clumped = False
        for pval_threshold in pval_thresholds:
            # TODO: Account for case where ld_clumping=False and pre_clumped=False
            if not pre_clumped:
                ss = threshold(ss=ss1, 
                               pval_threshold=pval_threshold, 
                               filter_by_varexp=filter_by_varexp)
            
            else:
                ss_clumped = threshold(ss=ss_clumped, 
                                       pval_threshold=pval_threshold, 
                                       filter_by_varexp=filter_by_varexp)
            if ld_clumping and not pre_clumped:
                ss_clumped = get_clumped_ss(ss=ss,
                                            ld_wind_cm=ld_wind_cm,
                                            use_ash=use_ash)
                pre_clumped = True # computational speed up to take advantage of clumped version for higher pval thresholds
                
            if save:
                phen_str = phen.replace(' ','_').lower()
                out_fname = f'{"clumped_" if ld_clumping else ""}'
                out_fname += "ash." if use_ash else "gwas."
                out_fname += f'{phen_str}.'
                out_fname += f'ld_wind_cm_{ld_wind_cm}.' if ld_clumping else ''
                out_fname += f'{"block_mhc." if block_mhc else ""}'
                out_fname += f'pval_{pval_threshold}.{"varexp_thresh." if filter_by_varexp else ""}' if pval_threshold <1 else ''
                out_fname += 'tsv.gz'
                ss_clumped.to_csv(f'{smiles_wd}/data/'+out_fname,sep='\t',index=False, compression='gzip')
                
            