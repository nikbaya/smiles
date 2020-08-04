#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 07:50:09 2020

Prune GWAS sumstats for ashR

@author: nbaya
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from time import time
import numpy as np

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

def pre_prune_qc(ss0, phen, maf=0.01):
    r'''
    Pre-pruning QC (does not include reading in ash results)
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
          
    if maf!=None:
        snp_ct_before = ss0.shape[0]
        ss0 = ss0[(ss0.eaf>maf)&(ss0.eaf<1-maf)]
        snp_ct_after = ss0.shape[0]
        print(f'# of SNPs removed with MAF<{maf}: {snp_ct_before-snp_ct_after}')
              
    if 'low_confidence_variant' in ss0.columns.values:
        ss0 = ss0[~ss0.low_confidence_variant]

    # note the number of SNPs with beta=0
    print(f'SNPs with beta=0: {ss0[ss0.beta==0].shape[0]}')
    
    ss0.loc[ss0.beta>=0,'raf'] = ss0.loc[ss0.beta>=0,'eaf']
    ss0.loc[ss0.beta<0,'raf'] = 1-ss0.loc[ss0.beta<0,'eaf']

    ss0.loc[ss0.index,'rbeta'] = ss0['beta'].abs() #beta is transformed to risk (or trait-increasing) allele effect size
        
    ss0['var_exp'] = 2*ss0.raf*(1-ss0.raf)*ss0.rbeta**2
    ss0['var_exp_se'] = 4*ss0.raf*(1-ss0.raf)*ss0.rbeta*ss0.se
        
    return ss0

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
        print(f'Number of variants in MHC: {len(mhc)}')
        non_mhc = ss[~((ss.chr==6)&(ss.pos>=start)&(ss.pos<=stop))]
        ss = non_mhc.append(mhc[mhc.pval==mhc.pval.min()])
    print(f'Post-filter # of variants keeping only one variant in MHC: {ss.shape[0]}')
    
    return ss    

def get_pruned_ss(ss, ld_wind_kb):
    r'''
    Prunes SNPs, agnostic to betas, p-values, etc. Only keeps first SNP in each
    window of `ld_wind_kb` kilobase pairs
    '''
    start = time()
    print(f'...Getting loci for LD window={ld_wind_kb}kb...\nPre-filter no. of variants = {ss.shape[0]}')
    ss[f'ld_index_{ld_wind_kb}'] = ss.chr.astype(str)+'-'+(ss.pos/(ld_wind_kb*1e3)).astype(int).astype(str)
    ss = ss.drop_duplicates(subset=f'ld_index_{ld_wind_kb}',keep='first') # if dataframe is sorted by position, this keeps the first SNP in terms of base pair position
    ss = ss.drop(columns=f'ld_index_{ld_wind_kb}')
    print(f'Time to prune: {round((time()-start)/60,3)} min')
    print(f'Post-filter # of variants (LD window={ld_wind_kb}kb): {ss.shape[0]}\n')
    return ss    

def plot_smile(ss, title, yaxis='var_exp'):
    plt.figure()
    plt.plot(ss.raf, ss[yaxis], '.')
    plt.xlabel('raf')
    plt.ylabel(yaxis)
    plt.yscale('log')
    plt.title(title)
    
def plot_varexp_hist(ss, phen, maf):
    r'''
    Plots histograms of variance explained and log(variance explained)
    '''
    if maf!=None:
        ss = ss[(ss.eaf>maf)&(ss.eaf<1-maf)]
        
    phen_str = phen_fname_dict[phen][1]
    
    for symlogy in [False, True]:
        plt.figure()
        plt.hist(ss.var_exp, 50)
        plt.title(f'{phen_str}\n({ss.shape[0]} variants{f" with MAF>{maf}" if maf!=None else ", no MAF filter"})')
        plt.xlabel('variance explained')
        plt.ylabel('density')    
        if symlogy:
            plt.yscale('symlog',linthreshy=1e0)
        plt.savefig(f'/Users/nbaya/Downloads/varexp_hist.{phen}.{f"maf_{maf}." if maf!=None else ""}{"symlogy." if symlogy else ""}png', dpi=300)

#    plt.figure()
#    plt.hist(ss.var_exp, 5000)
#    plt.xlim(left=0, right=1e-4)
#    plt.title(f'{phen_str}\n({ss.shape[0]} variants{f" with MAF>{maf}" if maf!=None else ", no MAF filter"})')
#    plt.xlabel('variance explained')
#    plt.ylabel('density')    
#    plt.savefig(f'/Users/nbaya/Downloads/tmp-varexp_hist.{phen}.{f"maf_{maf}." if maf!=None else ""}{"symlogy." if symlogy else ""}png', dpi=300)
    
    snp_ct_before = ss.shape[0]
    ss = ss[ss.var_exp>0]
    snp_ct_after = ss.shape[0]
    snp_ct_diff = snp_ct_before-snp_ct_after
    plt.figure()
    plt.hist(np.log10(ss.var_exp),50)
    plt.title(f'{phen_str}, {f"{snp_ct_diff} SNPs w/ varexp=0 removed" if snp_ct_diff>0 else ""})\n({ss.shape[0]} variants{f" with MAF>{maf}" if maf!=None else ", no MAF filter"})')
    plt.xlabel('log10(variance explained)')
    plt.ylabel('density')
    plt.savefig(f'/Users/nbaya/Downloads/log_varexp_hist.{phen}.{f"maf_{maf}." if maf!=None else ""}png', dpi=300)
        
def main(phen, maf):
    fname, phen_str = phen_fname_dict[phen]
    print(f'Starting phen {phen_str}, using file {fname}')
    print(f'Block MHC: {block_mhc}')
    print(f'LD window: {ld_wind_kb}kb')
    print(f'MAF: {maf}' if maf!=None else 'WARNING: No MAF filter')
    
    ss0 = pd.read_csv(f'{smiles_dir}/data/'+fname, sep='\t',compression='gzip')
    
    ss1 = pre_prune_qc(ss0=ss0, 
                       phen=phen, 
                       maf=maf)
    
#    for tmp_maf in ([maf] if maf!=None else [None, 0.01]):
#        plot_varexp_hist(ss=ss1, phen=phen, maf=tmp_maf)
    
    ss2 = get_blocked_mhc(ss1) if block_mhc else ss1
    
    if save_unpruned:
        unpruned_fname = f'gwas.{phen}.{"block_mhc." if block_mhc else ""}{f"maf_{maf}." if maf!=None else ""}tsv.gz'
        ss2[['chr','pos','A1','A2','beta','se', 'var_exp','var_exp_se']].to_csv(f'{smiles_dir}/data/{unpruned_fname}',compression='gzip',sep='\t', index=False)


    if phen=='breast_cancer':
        pval_threshold = 1e-5
        snp_ct_before = ss2.shape[0]
        ss2 = ss2[ss2.pval<=pval_threshold]
        print(f'NOTE: Removed {snp_ct_before-ss2.shape[0]} SNPs in "{phen}" with pval > {pval_threshold}')

    ss = get_pruned_ss(ss=ss2,
                       ld_wind_kb=ld_wind_kb)
    
#        plot_smile(ss, yaxis='var_exp', title=f'{phen}\npruned, {ld_wind_kb} kb window ({ss.shape[0]} SNPs)')
    
    if save_pruned:
        pruned_fname = f'pruned_gwas.{phen}.ld_wind_kb_{ld_wind_kb}.{"block_mhc." if block_mhc else ""}{f"maf_{maf}." if maf!=None else ""}tsv.gz'
        ss[['chr','pos','A1','A2','beta','se', 'var_exp','var_exp_se']].to_csv(f'{smiles_dir}/data/{pruned_fname}',compression='gzip',sep='\t', index=False)
    

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--phen', type=str, default=None, help='If included, only run this phenotype')
    args = parser.parse_args()
    
    ld_wind_kb = 500 # default: 500 kb
    block_mhc = True
    maf = 0.01 # to not include MAF filter, set `maf`=None
    
    save_unpruned = True
    save_pruned = True
    
    if args.phen==None:
        for phen, _ in phen_fname_dict.items():
            main(phen=phen, maf=maf)
    else:
        main(phen=args.phen, maf=maf)
