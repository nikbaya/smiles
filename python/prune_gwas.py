#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 07:50:09 2020

Prune GWAS sumstats for ashR

@author: nbaya
"""

import os
import pandas as pd
from time import time

if os.path.isdir('/Users/nbaya/'): # if working locally
    smiles_dir = "/Users/nbaya/Documents/lab/smiles"
else:
    smiles_dir = "/stanley/genetics/users/nbaya/smiles"

fname_dict = {
#             '50_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'Standing height', # UKB
#             '21001_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'BMI', # UKB
#             'Mahajan.NatGenet2018b.T2Dbmiadj.European.coding.tsv.gz':'T2D_bmiadj', #this is the edited version of the original BMI-adjusted data, some unnecessary columns were removed
#             'EUR.IBD.gwas_info03_filtered.assoc.coding.tsv.gz':'IBD',#,#EUR IBD from transethnic ancestry meta-analysis
#             'EUR.CD.gwas_info03_filtered.assoc.coding.tsv.gz':'CD', #EUR CD from transethnic ancestry meta-analysis
#             'EUR.UC.gwas_info03_filtered.assoc.coding.tsv.gz':'UC', #EUR UC from transethnic ancestry meta-analysis
#             'daner_PGC_SCZ43_mds9.coding.tsv.gz':'SCZ', #PGC data for SCZ (NOTE: The SNP effects were flipped when converting from odds ratio because daner files use odds ratios based on A1, presumably the ref allele, unlike UKB results which have betas based on the alt allele)
#             'AD_sumstats_Jansenetal_2019sept.coding.tsv.gz':'AD', #Alzheimer's disease meta-analysis
             'breastcancer.michailidou2017.b37.cleaned.coding.tsv.gz':'Breast cancer',
#             'MICAD.EUR.ExA.Consortium.PublicRelease.310517.cleaned.coding.tsv.gz': 'MICAD',
#             '30780_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz': 'LDL',
#             '30000_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz': 'WBC count',
#             '30010_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz': 'RBC count',
#             '30880_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz': 'Urate'
             '4080_irnt.gwas.imputed_v3.both_sexes.tsv.bgz':'Systolic BP',
             '4079_irnt.gwas.imputed_v3.both_sexes.tsv.bgz':'Diastolic BP',
             '30870_irnt.gwas.imputed_v3.both_sexes.tsv.bgz': 'Triglycerides',
             '30760_irnt.gwas.imputed_v3.both_sexes.tsv.bgz': 'HDL'
             }

def get_phen_str(phen):
    return phen.replace(' ','_').lower()

def pre_prune_qc(ss0, phen):
    r'''
    Pre-pruning QC (does not include MAF filter that pre_clump_qc has)
    '''
    
    assert 'se' in ss0.columns.values
    if 'chr' not in ss0.columns.values:
        ss0['chr'] = ss0.variant.str.split(':',n=1,expand=True).iloc[:,0]
    else:
        ss0['chr'] = ss0['chr'].astype(str)
    print(f'Filtering to autosomes...(variant qt: {ss0.shape[0]})')
    ss0 = ss0[ss0.chr.isin([str(x) for x in range(1,23)])] # only need to go up to chr 22 because genetic map only covers 22 chr
    print(f'Filtered to autosomes...(variant qt: {ss0.shape[0]})')
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
              
    if 'low_confidence_variant' in ss0.columns.values:
        ss0 = ss0[~ss0.low_confidence_variant]
    
    if phen=='Breast cancer':
        se_threshold = 1000
        snp_ct_before = ss0.shape[0]
        ss0 = ss0[ss0.se<se_threshold]
        print(f'NOTE: Removed {snp_ct_before-ss0.shape[0]} SNPs in "{phen}" with standard error > {se_threshold}')
        
    ss0.loc[ss0.beta>0,'raf'] = ss0.loc[ss0.beta>0,'eaf']
    ss0.loc[ss0.beta<0,'raf'] = 1-ss0.loc[ss0.beta<0,'eaf']

    ss0.loc[ss0.index,'rbeta'] = ss0['beta'].abs() #beta is transformed to risk (or trait-increasing) allele effect size
        
    ss0['var_exp'] = 2*ss0.raf*(1-ss0.raf)*ss0.rbeta**2
        
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


if __name__=="__main__":
    ld_wind_kb = 500 # default: 500 kb
    block_mhc = True
    save = True
    save_unpruned = True
    
    for fname, phen in fname_dict.items():
        print(f'Starting phen {phen}, using file {fname}')
        print(f'Block MHC: {block_mhc}')
        print(f'LD window: {ld_wind_kb}kb')
        
        ss0 = pd.read_csv(f'{smiles_dir}/data/'+fname, sep='\t',compression='gzip')
        
        ss1 = pre_prune_qc(ss0, phen)
        
        ss2 = get_blocked_mhc(ss1) if block_mhc else ss1
        
        if save_unpruned:
            phen_str = get_phen_str(phen)
            unpruned_fname = f'gwas.{get_phen_str(phen)}.{"block_mhc." if block_mhc else ""}tsv.gz'
            ss2[['chr','pos','beta','se']].to_csv(f'{smiles_dir}/data/{unpruned_fname}',compression='gzip',sep='\t', index=False)

        ss = get_pruned_ss(ss=ss2,
                           ld_wind_kb=ld_wind_kb)
        
        if save:
            pruned_fname = f'pruned_gwas.{get_phen_str(phen)}.ld_wind_kb_{ld_wind_kb}.{"block_mhc." if block_mhc else ""}tsv.gz'
            ss[['chr','pos','beta','se']].to_csv(f'{smiles_dir}/data/{pruned_fname}',compression='gzip',sep='\t', index=False)
