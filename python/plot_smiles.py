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

phen_dict = {#'50_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'Standing height'}#,
#             '21001_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'BMI'}#,
             '2443.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'Diabetes diagnosed'}#, #diabetes diagnosed by doctor
#             'Mahajan.NatGenet2018b.T2D.European.coding.tsv.gz':'T2D', #this is the edited version of the original data, some unnecessary columns were removed
#             'Mahajan.NatGenet2018b.T2Dbmiadj.European.coding.tsv.gz':'T2D_bmiadj'} #this is the edited version of the original BMI-adjusted data, some unnecessary columns were removed
#             '2443.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'Diabetes diagnosed'} #Directly from UKB
#        'pgc.scz.full.2012-04.tsv.gz':'SCZ'} #cases: 36,989, controls: 113,075, total: 150064
#        'EUR.IBD.gwas_info03_filtered.assoc.coding.tsv.gz':'IBD'}#EUR IBD from transethnic ancestry meta-analysis
#        'EUR.CD.gwas_info03_filtered.assoc.coding.tsv.gz':'CD', #EUR CD from transethnic ancestry meta-analysis
#        'EUR.UC.gwas_info03_filtered.assoc.coding.tsv.gz':'UC', #EUR UC from transethnic ancestry meta-analysis
#        'daner_PGC_SCZ43_mds9.coding.tsv.gz':'SCZ'} #PGC data for SCZ (NOTE: The SNP effects were flipped when converting from odds ratio because daner files use odds ratios based on A1, presumably the ref allele, unlike UKB results which have betas based on the alt allele)
#        '20544_2.gwas.imputed_v3.both_sexes.coding.tsv.bgz':'SCZ-UKB'} #UKB data for SCZ

highlight_coding = True
ld_clumping = True #only show the top hit (variant with lowest p-value) in a window of size {ld_window}
get_top_loci = False # only show top {n_top_loci} in the plot, and color by loci
n_top_loci = 10
ld_window = int(300e3) if get_top_loci else int(1000e3) #default for when ld_clumping is true: 100e3; default for when get_top_loci is true: 300e3
block_mhc = True
savefig = True


for filename, phen in phen_dict.items():
    print(f'Starting phen {phen}, using file {filename}')
    print(f'Highlight coding: {highlight_coding}')
    print(f'LD clumping: {ld_clumping}')
    print(f'Block MHC: {block_mhc}')
    print(f'Get top {n_top_loci} loci: {get_top_loci}')
    if ld_clumping or get_top_loci:
        print(f'LD window: {ld_window} ({ld_window/1e6}mb, {ld_window/1e3}kb)')
    
    
    ss0 = pd.read_csv(wd+'data/'+filename, sep='\t',compression='gzip')
    
    if 'chr' not in ss0.columns.values:
        ss0['chr'] = ss0.variant.str.split(':',n=1,expand=True).iloc[:,0]
    ss0['chr'] = ss0['chr'].astype(str)
    
    if 'pos' not in ss0.columns.values:
        ss0['pos'] = ss0.variant.str.split(':',n=2,expand=True).iloc[:,1]
    ss0['pos'] = ss0['pos'].astype(int)
        
    if 'n_complete_samples' in ss0.columns.values:
        print('renaming field "n_complete_samples" to "n"')
        ss0 = ss0.rename(columns={'n_complete_samples':'n'})
    n = int(ss0.n.mean())
    if ss0.n.std() != 0:
        print('WARNING: Number of samples varies across SNPs')
        
    if 'EAF' not in ss0.columns.values and all(x in ss0.columns.values for x in ['variant','minor_AF']):
        ss0['ref'] = ss0.variant.str.split(':',n=3,expand=True).iloc[:,2] #take the first allele listed in the variant ID string
        ss0.loc[ss0.minor_allele!=ss0.ref,'alt_af'] = ss0.loc[ss0.minor_allele!=ss0.ref,'minor_AF']
        ss0.loc[ss0.minor_allele==ss0.ref,'alt_af'] = 1-ss0.loc[ss0.minor_allele==ss0.ref,'minor_AF']
        ss0.loc[ss0.beta>0,'raf'] = ss0.loc[ss0.beta>0,'alt_af']
        ss0.loc[ss0.beta<0,'raf'] = 1-ss0.loc[ss0.beta<0,'alt_af']
    elif 'EAF' in ss0.columns.values:
#        ss0 = ss0.rename(columns={'EAF':'raf'}) # where raf is "risk allele frequency"
        ss0.loc[ss0.beta>0,'raf'] = ss0.loc[ss0.beta>0,'EAF']
        ss0.loc[ss0.beta<0,'raf'] = 1-ss0.loc[ss0.beta<0,'EAF']
    else:
        assert False, 'insufficient information to calculate risk allele frequency'

    ss0.loc[ss0.index,'rbeta'] = np.abs(ss0['beta']) #beta is transformed to risk (or trait-increasing) allele effect size
    
    if 'low_confidence_variant' in ss0.columns.values:
        ss0 = ss0[~ss0.low_confidence_variant]
#    
#    if ld_clumping and not get_top_loci:
#        if f'ld_index_{ld_window}' not in ss0.columns.values:
#            ss0[f'ld_index_{ld_window}'] = [f'{entry.chr}-{int(entry.pos/ld_window)}' for id,entry in ss0.iterrows()]
#        ss0 = ss0.loc[ss0.groupby(f'ld_index_{ld_window}')['pval'].idxmin()]
    
    
    for pval_threshold in [1e-5,5e-8,1e-8]:
        
        ss = ss0[ss0.pval < pval_threshold].reset_index()
        
        if ld_clumping:
            print(f'\nGetting loci for LD window={ld_window}...\nPre-filter # of variants (pval={pval_threshold}) = {ss.shape[0]}')
#            ss_ls = [None]*n_top_loci #list of dataframes with top loci
            ss_tmp = ss.copy() #dataframe from which top loci will be removed
            ss_keep = ss[[False]*len(ss)]
            i = 0 
            while len(ss_tmp)>0:
                ch, pos = ss_tmp[ss_tmp.pval==ss_tmp.pval.min()][['chr','pos']].values[0]
                ss_w_top_locus = ss_tmp[(ss_tmp['chr']==ch)&(ss_tmp['pos']==pos)&(ss_tmp['pos']==pos)].copy() #extract sentinel variant
                ss_w_top_locus['loci_rank'] = i
                ss_keep = ss_keep.append(ss_w_top_locus)
                ss_tmp = ss_tmp[~((ss_tmp['chr']==ch)&(ss_tmp['pos']>=pos-ld_window/2)&(ss_tmp['pos']<=pos+ld_window/2))].copy() #remove rows not around most significant hit
                i += 1
            ss = ss_keep    
            print(f'Post-filter # of variants (LD window={ld_window}): {ss.shape[0]}')

        if block_mhc: #combine all variants in MHC
            genes = pd.read_csv(wd+'data/cytoBand.txt',delim_whitespace=True,header=None,names=['chr','start','stop','region','gene'])
            mhc_region = genes[(genes.chr=='chr6')&(genes.region.str.contains('p21.'))]
            start = min(mhc_region.start)
            stop = max(mhc_region.stop)
            mhc = ss[(ss.chr=='6')&(ss.pos>=start)&(ss.pos<=stop)]
            if len(mhc)>0: #if there are variants in MHC
                print(f'Number of variants in MHC: {len(mhc)}')
                non_mhc = ss[~((ss.chr=='6')&(ss.pos>=start)&(ss.pos<=stop))]
                ss = non_mhc.append(mhc[mhc.pval==mhc.pval.min()])
            print(f'Post-filter # of variants keeping only one variant in MHC: {ss.shape[0]}')
            

        if get_top_loci:
#            ss_ls = [None]*n_top_loci #list of dataframes with top loci
            ss_tmp = ss.copy() #dataframe from which top loci will be removed
            ss_keep = ss[[False]*len(ss)]
            i = 0 
            while len(ss_tmp)>0 and i<n_top_loci:
                ch, pos = ss_tmp[ss_tmp.pval==ss_tmp.pval.min()][['chr','pos']].values[0]
                ss_w_top_locus = ss_tmp[(ss_tmp['chr']==ch)&(ss_tmp['pos']>=pos-ld_window/2)&(ss_tmp['pos']<=pos+ld_window/2)].copy() #extract rows around most significant hit
                ss_w_top_locus['loci_rank'] = i
                ss_keep = ss_keep.append(ss_w_top_locus)
                ss_tmp = ss_tmp[~((ss_tmp['chr']==ch)&(ss_tmp['pos']>=pos-ld_window/2)&(ss_tmp['pos']<=pos+ld_window/2))].copy() #remove rows not around most significant hit
                i += 1
            ss = ss_keep
            ss = ss.sort_values(by='index') #unnecessary step, but restores the order of variants

            
        for maf in [0, 0.01]:
            if maf != 0:
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
            plt.ylabel('Estimated effect size')
            plt.xlim([-0.01, 1.01])
            plt.title(f'AF vs. Effect Size\nphen: {phen}, n: {n}, pval threshold: {pval_threshold}'
                      +(f', maf>{maf}' if maf!=0 else '')+(f', ld block: {int(ld_window/1000)}kb' if ld_clumping else '')
                      +(f', loci window: {int(ld_window/1000)}kb' if get_top_loci else ''))
            suffix = ''
            if 'coding' in ss.columns.values  and highlight_coding:
                suffix = '.highlightcoding'
            if ld_clumping and not get_top_loci:
                suffix += f'.ldclumping{ld_window}'
            if get_top_loci:
                suffix += f'.top{n_top_loci}loci'
            if block_mhc:
                suffix += f'.block_mhc'
            if savefig:
                fig.savefig(wd+f'smiles/plots/{file_prefix}.raf_effectsize{suffix}.png',dpi=600)
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
                      +(f', maf>{maf}' if maf!=0 else '')+(f', ld block: {int(ld_window/1000)}kb' if ld_clumping else '')
                      +(f', loci window: {int(ld_window/1000)}kb' if get_top_loci else ''))
            suffix = ''
            if 'coding' in ss.columns.values  and highlight_coding:
                suffix = '.highlightcoding'
            if ld_clumping:
                suffix += f'.ldclumping{ld_window}'
            if get_top_loci:
                suffix += f'.top{n_top_loci}loci'
            if block_mhc:
                suffix += f'.block_mhc'
            if savefig:
                fig.savefig(wd+f'smiles/plots/{file_prefix}.raf_varianceexplained{suffix}.png',dpi=600)
            plt.close()
    
            #variance explained, colored by chromosome
            if not get_top_loci:
                fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))
                ss['var_exp'] = 2*ss.raf*(1-ss.raf)*ss.rbeta**2
                for i, ch in enumerate([*range(1,23)]+['X']):
                    if 'coding' in ss.columns.values and highlight_coding:
                        ax.plot(ss[(ss.chr==str(ch))&(~ss.coding)].raf, ss[(ss.chr==str(ch))&(~ss.coding)].var_exp,'.',ms=4,c=colors[i%10])
                        ax.plot(ss[(ss.chr==str(ch))&(ss.coding)].raf, ss[(ss.chr==str(ch))&(ss.coding)].var_exp, 'o',ms=10,markerfacecolor='none', markeredgewidth=0.75,c=colors[i%10])
                    else:
                        raf = ss[ss.chr==str(ch)]['raf']
                        rbeta = ss[ss.chr==str(ch)]['rbeta']
                        ax.plot(raf, 2*raf*(1-raf)*rbeta**2,'.',ms=2)
                plt.xlabel('Risk allele frequency')
                plt.ylabel('Variance explained')
                plt.xlim([-0.01, 1.01])
                plt.title(f'AF vs. Variance Explained, colored by chromosome\nphen: {phen}, n: {n}, pval threshold: {pval_threshold}'
                          +(f', maf>{maf}' if maf!=0 else '')+(f', ld block: {int(ld_window/1000)}kb' if ld_clumping else ''))
                suffix = ''
                if 'coding' in ss.columns.values  and highlight_coding:
                    suffix = '.highlightcoding'
                if ld_clumping:
                    suffix += f'.ldclumping{ld_window}'
                if block_mhc:
                    suffix += f'.block_mhc'
                if savefig:
                    fig.savefig(wd+f'smiles/plots/{file_prefix}.raf_varianceexplained.coloredbychr{suffix}.png',dpi=600)
                plt.close()