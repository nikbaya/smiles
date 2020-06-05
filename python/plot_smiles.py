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
from sklearn.linear_model import LinearRegression
import scipy.stats as stats
from time import time

smiles_wd = "/Users/nbaya/Documents/lab/smiles/"

phen_dict = {
#            '50_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'Standing height',
#             '21001_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'BMI',
####             '2443.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'Diabetes diagnosed'}#, #diabetes diagnosed by doctor
####             'Mahajan.NatGenet2018b.T2D.European.coding.tsv.gz':'T2D', #this is the edited version of the original data, some unnecessary columns were removed
             'Mahajan.NatGenet2018b.T2Dbmiadj.European.coding.tsv.gz':'T2D_bmiadj', #this is the edited version of the original BMI-adjusted data, some unnecessary columns were removed
####            '2443.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'Diabetes diagnosed'} #Directly from UKB
####            'pgc.scz.full.2012-04.tsv.gz':'SCZ'} #cases: 36,989, controls: 113,075, total: 150064
            'EUR.IBD.gwas_info03_filtered.assoc.coding.tsv.gz':'IBD',#,#EUR IBD from transethnic ancestry meta-analysis
#            'EUR.CD.gwas_info03_filtered.assoc.coding.tsv.gz':'CD', #EUR CD from transethnic ancestry meta-analysis
#            'EUR.UC.gwas_info03_filtered.assoc.coding.tsv.gz':'UC', #EUR UC from transethnic ancestry meta-analysis
#            'daner_PGC_SCZ43_mds9.coding.tsv.gz':'SCZ', #PGC data for SCZ (NOTE: The SNP effects were flipped when converting from odds ratio because daner files use odds ratios based on A1, presumably the ref allele, unlike UKB results which have betas based on the alt allele)
###             '20544_2.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'SCZ-UKB'} #UKB data for SCZ
#            'AD_sumstats_Jansenetal_2019sept.coding.tsv.gz':'AD',#,
#               '50_irnt.gwas.imputed_v3.male.coding.tsv.bgz':'Standing height male',
#               '50_irnt.gwas.imputed_v3.female.coding.tsv.bgz':'Standing height female'}
            }
## GWAS catalog
# myocardial infarction
# lung cancer
# breast cancer 
# multiple sclerosis

# other cancers?
# Stroke?
# Immunochip chip?
# RA (?) 

#for i in range(1,23):
#    genmap_chr = genmap[genmap.chr==i].sort_values(by='pos') # for comptuational speed improvement, save all genmap_chr to a list
#    plt.plot(genmap_chr.pos/genmap_chr.pos.max(), genmap_chr.cm/genmap_chr.cm.max())


highlight_coding = True
ld_clumping = True #only show the top hit (variant with lowest p-value) in a window of size {ld_wind_kb}
get_top_loci = False # only show top {n_top_loci} in the plot, and color by loci
n_top_loci = int(1e10)
ld_wind_kb = 300 #int(300e3) if get_top_loci else int(1000e3) #measured in kb; default for when ld_clumping is true: 100e3; default for when get_top_loci is true: 300e3
ld_wind_cm = 0.6 # 1 Mb = 1 cM -> 600 Kb = 0.6 cM
block_mhc = True
savefig = False #True
showfig = True
save_clumped = False #True #whether to save the post-clumping version of the summary stats dataframe
filter_by_varexp = True

#for i in range(1,23):
#    genmap = pd.read_csv(f'{smiles_wd}data/genetic_map_chr{i}_combined_b37.txt', 
#                         delim_whitespace=True, low_memory=False)
#    genmap['chr'] = i
#    genmap[['chr', 'position', 'Genetic_Map(cM)']].rename({'Genetic_Map(cM)':'cM', 'position':'pos'}, axis=1).to_csv(
#            f'{smiles_wd}data/genetic_map_chr{i}_combined_b37.txt.gz', 
#            sep=' ', compression='gzip', index=False)

genmap = pd.read_csv(f'{smiles_wd}data/genetic_map_combined_b37.txt.gz',
                     sep=' ', names=['chr','pos','rate','cm'], compression='gzip')

genmap_chr_list = [genmap[genmap.chr==chr].sort_values(by='pos') for chr in range(1,23)]

def impute_cm(genmap_chr, chr, pos):
    r'''
    Extrapolates genetic distance in cM based on position `pos` in base pairs 
    on a given chromosome `chr`, using a recombination map DataFrame `genmap_chr`
    '''
    if pos in genmap_chr.pos.values:
        return genmap_chr[genmap_chr.pos==pos].cm.values[0]
    else:
        a_pos, a_cm = genmap_chr[genmap_chr.pos<pos].tail(1)[['pos','cm']].values[0] # throws an IndexError if pos < min(genmap_chr.pos)
        b_pos, b_cm = genmap_chr[genmap_chr.pos>pos].head(1)[['pos','cm']].values[0] # throws an IndexError if pos > max(genmap_chr.pos)
        cm = a_cm + (b_cm-a_cm)*(pos-a_pos)/(b_pos-a_pos)
        return cm

    
def impute_pos(genmap_chr, chr, cm):
    r'''
    Extrapolates position in base pairs based on the genetic distance `cm` in cM
    along a given chromosome `chr`, using a recombination map DataFrame `genmap_chr`
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
    
    
for fname, phen in phen_dict.items():
    print(f'Starting phen {phen}, using file {fname}')
    print(f'Highlight coding: {highlight_coding}')
    print(f'LD clumping: {ld_clumping}')
    print(f'Block MHC: {block_mhc}')
    print(f'Get top {n_top_loci} loci: {get_top_loci}')
    if ld_clumping or get_top_loci:
        print(f'LD window: {ld_wind_cm}cM')
#        print(f'LD window: {ld_wind_kb/1e3}mb, {ld_wind_kb}kb')
    
    
    ss0 = pd.read_csv(smiles_wd+'data/'+fname, sep='\t',compression='gzip')
    
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
        n = round(ss0[sample_ct_col].mean())
        if ss0[sample_ct_col].std() != 0:
            print('WARNING: Number of samples (using column "{sample_ct_col}") varies across SNPs, these are likely meta-analyzed sumstats.')
    else: 
        print('WARNING: phen {phen} is missing number of samples.')
        n = np.nan
    
    if 'EAF' not in ss0.columns.values:
        if all(x in ss0.columns.values for x in ['variant','minor_AF']):
            ss0['ref'] = ss0.variant.str.split(':',n=3,expand=True).iloc[:,2] #take the first allele listed in the variant ID string
            ss0.loc[ss0.minor_allele!=ss0.ref,'alt_af'] = ss0.loc[ss0.minor_allele!=ss0.ref,'minor_AF']
            ss0.loc[ss0.minor_allele==ss0.ref,'alt_af'] = 1-ss0.loc[ss0.minor_allele==ss0.ref,'minor_AF']
    #        ss0.loc[ss0.beta>0,'raf'] = ss0.loc[ss0.beta>0,'alt_af'] #only necessary for UKB, which runs a GWAS on the alternate allele (betas and OR are wrt alt allele...?)
    #        ss0.loc[ss0.beta<0,'raf'] = 1-ss0.loc[ss0.beta<0,'alt_af']
            ss0 = ss0.rename(columns={'alt_af':'EAF'}) # this is true for UKB, which is the only dataset without EAF
        else:
            assert False, 'insufficient information to calculate risk allele frequency'
    
    # remove SNPs with beta=NA
    snp_ct_before = ss0.shape[0]
    ss0 = ss0.dropna(axis=0, subset=['beta'])
    snp_ct_after = ss0.shape[0]
    print(f'# of SNPs removed with missing beta: {snp_ct_before-snp_ct_after}')
    
    # remove SNPs with EAF=NA
    snp_ct_before = ss0.shape[0]
    ss0 = ss0.dropna(axis=0, subset=['EAF'])
    snp_ct_after = ss0.shape[0]
    print(f'# of SNPs removed with missing EAF: {snp_ct_before-snp_ct_after}')
          
    # remove invariant SNPs (MAF=0)
    snp_ct_before = ss0.shape[0]
    ss0 = ss0[(ss0.EAF>0) & (ss0.EAF<1)]
    snp_ct_after = ss0.shape[0]
    print(f'# of SNPs removed with MAF=0: {snp_ct_before-snp_ct_after}')
    
    ss0.loc[ss0.beta>0,'raf'] = ss0.loc[ss0.beta>0,'EAF']
    ss0.loc[ss0.beta<0,'raf'] = 1-ss0.loc[ss0.beta<0,'EAF']

    ss0.loc[ss0.index,'rbeta'] = np.abs(ss0['beta']) #beta is transformed to risk (or trait-increasing) allele effect size
    
    if 'low_confidence_variant' in ss0.columns.values:
        ss0 = ss0[~ss0.low_confidence_variant]
    
    if filter_by_varexp:
        ss0['var_exp'] = 2*ss0.raf*(1-ss0.raf)*ss0.rbeta**2
#    
#    if ld_clumping and not get_top_loci:
#        if f'ld_index_{ld_wind_kb}' not in ss0.columns.values:
#            ss0[f'ld_index_{ld_wind_kb}'] = [f'{entry.chr}-{int(entry.pos/ld_wind_kb)}' for id,entry in ss0.iterrows()]
#        ss0 = ss0.loc[ss0.groupby(f'ld_index_{ld_wind_kb}')['pval'].idxmin()]
    
#    pval_thresholds = sorted([1, 1e-5, 1e-6, 1e-7, 5e-8, 1e-8],reverse=True) # must be sorted largest to smallest to be able to take advantage of time-saving step where ss0 = ss.copy()
    pval_thresholds = sorted([1e-5, 1e-6, 1e-7, 5e-8, 1e-8],reverse=True)
    
    for pval_threshold in pval_thresholds:
        if filter_by_varexp and pval_threshold<1: # this filtering is unnecessary when pval_threshold=1
            print('filtering by variance-explained')
            phen_str = phen.replace(' ','_').lower()
            clumped = pd.read_csv(f'{smiles_wd}data/clumped_gwas.{phen_str}.ld_wind_cm_{ld_wind_cm}.{"block_mhc." if block_mhc else ""}pval_1e-05.tsv.gz',
                                  compression='gzip', sep='\t')
            clumped.loc[clumped.pval==0, 'pval'] = 5e-324
            clumped = clumped[(clumped.raf!=0)&(clumped.raf!=1)]
            clumped['var_exp'] = 2*clumped.raf*(1-clumped.raf)*clumped.rbeta**2
            clumped['n_hat'] = stats.chi2.isf(q=clumped.pval,df=1)/clumped.var_exp
            n_hat_mean = clumped.n_hat.mean()
            varexp_thresh = stats.chi2.isf(q=pval_threshold,df=1)/n_hat_mean
            print(f'estimated n_eff: {n_hat_mean}\nvarexp thresh: {varexp_thresh}')
            ss = ss0[ss0.var_exp>varexp_thresh].reset_index(drop=True)
            
        elif not filter_by_varexp and pval_threshold<1:        
            ss = ss0[ss0.pval < pval_threshold].reset_index(drop=True) # filter by p-value
        else:
            ss = ss0
            
        if block_mhc: #combine all variants in MHC
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
            
#            print(f'\nGetting loci for LD window={ld_wind_kb}...\nPre-filter # of variants (pval={pval_threshold}) = {ss.shape[0]}')
##            ss_ls = [None]*n_top_loci #list of dataframes with top loci
#            ss_tmp = ss.copy() #dataframe from which top loci will be removed
#            ss_keep = ss[[False]*len(ss)]
#            i = 0 
#            while len(ss_tmp)>0:
#                ch, pos = ss_tmp[ss_tmp.pval==ss_tmp.pval.min()][['chr','pos']].values[0]
#                ss_w_top_locus = ss_tmp[(ss_tmp['chr']==ch)&(ss_tmp['pos']==pos)&(ss_tmp['pos']==pos)].copy() #extract sentinel variant
#                ss_w_top_locus['loci_rank'] = i
#                ss_keep = ss_keep.append(ss_w_top_locus)
#                ss_tmp = ss_tmp[~((ss_tmp['chr']==ch)&(ss_tmp['pos']>=pos-ld_wind_kb*1e3/2)&(ss_tmp['pos']<=pos+ld_wind_kb*1e3/2))].copy() #remove rows not around most significant hit
#                i += 1
#            ss = ss_keep    
#            print(f'Post-filter # of variants (LD window={ld_wind_kb}): {ss.shape[0]}')

            
        ss_noncoding = ss[~ss.coding]
        ss_coding = ss[ss.coding]
        ss_final = None
        
#        if get_top_loci:
#            for coding, ss_tmp in [[False, ss_noncoding], [True, ss_coding]]:
#                ss_keep = ss[[False]*len(ss)]
#                i = 0 
#                while len(ss_tmp)>0 and i<n_top_loci:
#                    ch, pos = ss_tmp[ss_tmp.pval==ss_tmp.pval.min()][['chr','pos']].values[0]
#                    ss_w_top_locus = ss_tmp[(ss_tmp['chr']==ch)&(ss_tmp['pos']>=pos-ld_wind_kb*1e3/2)&(ss_tmp['pos']<=pos+ld_wind_kb*1e3/2)].copy() #extract rows around most significant hit
#                    ss_w_top_locus['loci_rank'] = i
#                    ss_keep = ss_keep.append(ss_w_top_locus)
#                    ss_tmp = ss_tmp[~((ss_tmp['chr']==ch)&(ss_tmp['pos']>=pos-ld_wind_kb*1e3/2)&(ss_tmp['pos']<=pos+ld_wind_kb*1e3/2))].copy() #remove rows not around most significant hit
#                    i += 1
#                ss_tmp = ss_keep
#                ss_tmp = ss_tmp.sort_values(by='index') #unnecessary step, but restores the order of variants
#                ss_final = ss_tmp if ss_final is None else ss_final.append(ss_tmp)                
#            ss = ss_final
        if ld_clumping:
            start = time()
            print(f'\nGetting loci for LD window={ld_wind_cm}cM ...\nPre-filter # of variants (pval={pval_threshold}) = {ss.shape[0]}')
            for coding, ss_tmp in [[False, ss_noncoding], [True, ss_coding]]:
#                success = 0
#                fail = 0
                ss_keep = ss[[False]*len(ss)]
                i = 0 
                while len(ss_tmp)>0:
                    chr, pos = ss_tmp[ss_tmp.pval==ss_tmp.pval.min()][['chr','pos']].values[0] # sentinel variant is most significant non-clumped variant
                    genmap_chr = genmap_chr_list[chr-1]
                    if pos < genmap_chr.pos.min():
                        a_pos = 1 # start of window around sentinel
                        b_pos = impute_pos(genmap_chr=genmap_chr, chr=chr, cm=ld_wind_cm/2) # end of window around sentinel
                    elif pos > genmap_chr.pos.max():
                        a_pos = impute_pos(genmap_chr=genmap_chr, chr=chr, cm=genmap_chr.cm.max()-ld_wind_cm/2) # start of window around sentinel
                        b_pos = np.Inf # end of window around sentinel
                    else:
                        cm = impute_cm(genmap_chr=genmap_chr, chr=chr, pos=pos) # genetic distance of sentinel in cM
                        a_pos = impute_pos(genmap_chr=genmap_chr, chr=chr, cm=cm-ld_wind_cm/2) # start of window around sentinel
                        b_pos = impute_pos(genmap_chr=genmap_chr, chr=chr, cm=cm+ld_wind_cm/2) # end of window around sentinel
                    ss_w_top_locus = ss_tmp[(ss_tmp['chr']==chr)&(ss_tmp['pos']==pos)].copy() #extract sentinel variant, NOTE: This step must come before extracting rows in window around sentinel
                    ss_w_top_locus['loci_rank'] = i
                    ss_tmp = ss_tmp[~((ss_tmp['chr']==chr)&(ss_tmp['pos']>=a_pos)&(ss_tmp['pos']<=b_pos))].copy() #extract rows not in window around current sentinel variant
                    ss_keep = ss_keep.append(ss_w_top_locus, sort=False)
                    i += 1
#                print(f'success rate (pruning w/ cM): {success/(success+fail)}')
                ss_tmp = ss_keep.sort_index() #unnecessary step, but restores the order of variants
                ss_final = ss_tmp if ss_final is None else ss_final.append(ss_tmp, sort=False)                
            ss = ss_final
            print(f'Time to clump: {round((time()-start)/60,3)} min')
        print(f'Post-filter # of variants (LD window={ld_wind_cm}cM): {ss.shape[0]}\n')
            
        if save_clumped:
            phen_str = phen.replace(' ','_').lower()
            clumped_fname = f'clumped_gwas.{phen_str}.ld_wind_cm_{ld_wind_cm}.{"block_mhc." if block_mhc else ""}pval_{pval_threshold}.{"varexp_thresh." if filter_by_varexp else ""}tsv.gz'
            ss.to_csv(smiles_wd+'data/'+clumped_fname,sep='\t',index=False, compression='gzip')
            
        if pval_threshold == 1e-8:
            for maf in [0]:
                
                ss = ss[(ss.raf>maf)&(ss.raf<1-maf)]
                
                file_prefix = f'{phen}.pval_{pval_threshold}'+(f'.maf_{maf}' if maf!=0 else '')
                colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
                
                #effect size beta
                fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))
                if 'coding' in ss.columns.values and highlight_coding:
                    if get_top_loci:
                        for loci_i in range(n_top_loci):
                            ax.plot(ss[(~ss.coding)&(ss.loci_rank==loci_i)].raf, ss[~ss.coding&(ss.loci_rank==loci_i)].rbeta,'.',ms=4,c=colors[loci_i%10]) #plot noncoding variants
                            ax.plot(ss[ss.coding&(ss.loci_rank==loci_i)].raf, ss[ss.coding&(ss.loci_rank==loci_i)].rbeta, 'o',ms=10,markerfacecolor='none', markeredgewidth=0.75,c=colors[loci_i%10]) #plot coding variants
                    else:
                        ax.plot(ss[~ss.coding].raf, ss[~ss.coding].rbeta,'.',ms=4,c='#1f77b4') #plot noncoding variants
                        ax.plot(ss[ss.coding].raf, ss[ss.coding].rbeta, 'o',ms=10,markerfacecolor='none', markeredgewidth=0.75) #plot coding variants
                        plt.legend(['non-coding','coding'], loc=9)
                else:
                    ax.plot(ss.raf, ss.rbeta,'.',ms=2,alpha=1)
                plt.xlabel('Risk allele frequency')
                plt.ylabel('Effect size')
                plt.xlim([-0.01, 1.01])
                plt.title(f'AF vs. Effect Size\nphen: {phen}, n: {n}, pval threshold: {pval_threshold}'
                          +(f', maf>{maf}' if maf!=0 else '')+(f', ld block: {int(ld_wind_kb)}kb' if ld_clumping else '')
                          +(f', loci window: {int(ld_wind_kb)}kb' if get_top_loci else ''))
                suffix = ''
                if 'coding' in ss.columns.values  and highlight_coding:
                    suffix = '.highlightcoding'
                if ld_clumping and not get_top_loci:
                    suffix += f'.ldclumping{ld_wind_kb}'
                if get_top_loci:
                    suffix += f'.top{n_top_loci}loci'
                if block_mhc:
                    suffix += f'.block_mhc'
                if savefig:
                    fig.savefig(smiles_wd+f'smiles/plots/{file_prefix}.raf_effectsize{suffix}.png',dpi=600)
                if not showfig:
                    plt.close()
                    
                #variance explained
                fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))    
                if 'coding' in ss.columns.values  and highlight_coding:
                    ss['var_exp'] = 2*ss.raf*(1-ss.raf)*ss.rbeta**2
                    if get_top_loci:
                        for loci_i in range(n_top_loci):
                            ax.plot(ss[~ss.coding&(ss.loci_rank==loci_i)].raf, ss[~ss.coding&(ss.loci_rank==loci_i)].var_exp,'.',ms=4,c=colors[loci_i%10]) #plot noncoding variants
                            ax.plot(ss[ss.coding&(ss.loci_rank==loci_i)].raf,  ss[ss.coding&(ss.loci_rank==loci_i)].var_exp, 'o',ms=10,markerfacecolor='none', markeredgewidth=0.75,c=colors[loci_i%10]) #plot coding variants
                    else:
                        ax.plot(ss[~ss.coding].raf, ss[~ss.coding].var_exp,'.',ms=4,c='#1f77b4')
                        ax.plot(ss[ss.coding].raf, ss[ss.coding].var_exp, 'o',ms=10,markerfacecolor='none', markeredgewidth=0.75,c='#1f77b4')
                        plt.legend(['non-coding','coding'],loc=9)
                else:
                    ax.plot(ss.raf, 2*ss.raf*(1-ss.raf)*ss.rbeta**2,'.',ms=2)
                plt.xlabel('Risk allele frequency')
                plt.ylabel('Variance explained')
                plt.xlim([-0.01, 1.01])
                plt.title(f'AF vs. Variance Explained\nphen: {phen}, n: {n}, pval threshold: {pval_threshold}'
                          +(f', maf>{maf}' if maf!=0 else '')+(f', ld block: {int(ld_wind_kb)}kb' if ld_clumping else '')
                          +(f', loci window: {int(ld_wind_kb)}kb' if get_top_loci else ''))
                suffix = ''
                if 'coding' in ss.columns.values  and highlight_coding:
                    suffix = '.highlightcoding'
                if ld_clumping:
                    suffix += f'.ldclumping{ld_wind_kb}'
                if get_top_loci:
                    suffix += f'.top{n_top_loci}loci'
                if block_mhc:
                    suffix += f'.block_mhc'
                if savefig:
                    fig.savefig(smiles_wd+f'smiles/plots/{file_prefix}.raf_varianceexplained{suffix}.png',dpi=600)
                if not showfig:
                    plt.close()
    
                #variance explained, colored by chromosome
                if not get_top_loci:
                    fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))
                    ss['var_exp'] = 2*ss.raf*(1-ss.raf)*ss.rbeta**2
                    for i, ch in enumerate([*range(1,23)]):
                        if 'coding' in ss.columns.values and highlight_coding:
                            ax.plot(ss[(ss.chr==ch)&(~ss.coding)].raf, ss[(ss.chr==ch)&(~ss.coding)].var_exp,'.',ms=4,c=colors[i%10])
                            ax.plot(ss[(ss.chr==ch)&(ss.coding)].raf, ss[(ss.chr==ch)&(ss.coding)].var_exp, 'o',ms=10,markerfacecolor='none', markeredgewidth=0.75,c=colors[i%10])
                        else:
                            raf = ss[ss.chr==str(ch)]['raf']
                            rbeta = ss[ss.chr==str(ch)]['rbeta']
                            ax.plot(raf, 2*raf*(1-raf)*rbeta**2,'.',ms=2)
                    plt.xlabel('Risk allele frequency')
                    plt.ylabel('Variance explained')
                    plt.xlim([-0.01, 1.01])
                    plt.title(f'AF vs. Variance Explained, colored by chromosome\nphen: {phen}, n: {n}, pval threshold: {pval_threshold}'
                              +(f', maf>{maf}' if maf!=0 else '')+(f', ld block: {int(ld_wind_kb)}kb' if ld_clumping else ''))
                    suffix = ''
                    if 'coding' in ss.columns.values  and highlight_coding:
                        suffix = '.highlightcoding'
                    if ld_clumping:
                        suffix += f'.ldclumping{ld_wind_kb}'
                    if block_mhc:
                        suffix += f'.block_mhc'
                    if savefig:
                        fig.savefig(smiles_wd+f'smiles/plots/{file_prefix}.raf_varianceexplained.coloredbychr{suffix}.png',dpi=600)
                    if not showfig:
                        plt.close()