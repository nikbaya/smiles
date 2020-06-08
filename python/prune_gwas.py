#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 07:50:09 2020

Prune GWAS sumstats for ashR

@author: nbaya
"""

import pandas as pd
from time import time
import multiprocessing as mp
from functools import partial

smiles_wd = "/Users/nbaya/Documents/lab/smiles/"

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

def get_blocked_mhc(ss):
    r'''
    Removes all but the most significant variant in the MHC region
    '''
    genes = pd.read_csv(smiles_wd+'data/cytoBand.txt',delim_whitespace=True,header=None,names=['chr','start','stop','region','gene'])
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
    ss_final = None
    print(f'...Getting loci for LD window={ld_wind_kb}kb...\nPre-filter no. of variants = {ss.shape[0]}')
    pool = mp.Pool(mp.cpu_count()-1)
    func = partial(prune_chrom, ss, ld_wind_kb)
    ss_keep_list = pool.map(func, range(1,23))
    pool.close()
    ss_final = pd.concat(ss_keep_list)
    print(f'Time to clump: {round((time()-start)/60,3)} min')
    print(f'Post-filter # of variants (LD window={ld_wind_kb}kb): {ss_final.shape[0]}\n')

    return ss_final
  
def prune_chrom(ss, ld_wind_kb, chrom):
    r'''
    Prunes variants on a single chromosome `chrom` for a given sumstats 
    dataframe `ss`, such that only the first SNP in a `ld_wind_kb` kb window is kept.
    '''
    ss = ss[ss.chr==chrom]
    ss[f'ld_index_{ld_wind_kb}'] = [f'{entry.chr}-{int(entry.pos/(ld_wind_kb*1e3))}' for id,entry in ss.iterrows()]
    ss = ss.loc[ss.groupby(f'ld_index_{ld_wind_kb}')['pval'].idxmin()]    
    return ss



if __name__=="__main__":
    ld_wind_kb = 500 # default: 500 kb
    block_mhc = False
    
    for fname, phen in fname_dict.items():
        print(f'Starting phen {phen}, using file {fname}')

        print(f'Block MHC: {block_mhc}')
        print(f'LD window: {ld_wind_kb}kb')
        
        ss0 = pd.read_csv(smiles_wd+'data/'+fname, sep='\t',compression='gzip')
        
        ss = get_blocked_mhc(ss0) if block_mhc else ss0

        ss = get_pruned_ss(ss=ss,
                           ld_wind_kb=ld_wind_kb)
        
        break

