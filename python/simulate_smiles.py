#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:16:51 2019

Simulates loci by loci to understand LD structure

@author: nbaya
"""


import hail as hl
#import pandas as pd
#import matplotlib.pyplot as plt
import numpy as np
import subprocess
import argparse
import requests
#from hail.experimental.ldscsim import calculate_phenotypes
url = 'https://raw.githubusercontent.com/nikbaya/ldscsim/master/ldscsim.py'
r = requests.get(url).text
exec(r)
calculate_phenotypes=calculate_phenotypes
hl.init(log='/tmp/foo.log')


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--phen', type=str, required=True, help="phenotype. Options: 50_irnt, 21001_irnt")

args = parser.parse_args()

wd = 'gs://nbaya/smiles/'


def get_sim_all_loci_mt(phen, n_top_loci, ld_window, pval_threshold, maf, block_mhc, 
                        random_betas=False, n_top_loci_sim=None):
    r'''
    If `random_betas`=True, new top loci are created by randomly choosing new
    sentinel variants. If `random_betas`=False, loci are defined by sentinel 
    variants from GWAS.
    n_top_loci_sim : Number of top loci to include in simulation
    '''
    n_top_loci_sim = n_top_loci if n_top_loci_sim== None else n_top_loci_sim
    print('\n####################')
    print(f'Phenotype: {phen}')
    print(f'pval_threshold: {pval_threshold}')
    print(f'ld_window: {ld_window}')
    print(f'maf: {maf}')
    print(f'block_mhc: {block_mhc}')
    print(f'n_top_loci: {n_top_loci}')
    print(f'n_top_loci_sim: {n_top_loci_sim}')
    print(f'random_betas: {random_betas}')
    print('####################\n')
    suffix = f'ldwindow{int(ld_window/1e3)}kb{".block_mhc" if block_mhc else ""}{".random_betas" if random_betas else ""}'
    sim_all_loci_path = wd +phen + f'.simulation.top{n_top_loci_sim}loci.{suffix}.mt'
    try:
        subprocess.check_output(
            [
                f'gsutil',
                'ls',
                sim_all_loci_path+'/_SUCCESS']) is not None
        sim_all_loci = hl.read_matrix_table(sim_all_loci_path)
        sim_all_loci.describe()
    except BaseException:
        print('\n... Getting sim_all_loci matrix table ...')
        print(f'Path for sim_all_loci matrix table: {sim_all_loci_path}\n')
        mt_all_loci_path = wd +phen + f'.genotypes.top{n_top_loci}loci.{suffix}.mt'
        try:
            subprocess.check_output(
                [
                    f'gsutil',
                    'ls',
                    mt_all_loci_path+'/_SUCCESS']) is not None
            mt_all_loci = hl.read_matrix_table(mt_all_loci_path)
        except BaseException:
            print('\n... Getting mt_all_loci matrix table ...')
            print(f'Path for mt_all_loci matrix table: {mt_all_loci_path}\n')

            mt0 = hl.read_matrix_table(f'gs://nbaya/split/ukb31063.qc_pos_variants.gwas_samples_repart.mt')
            
            # restrict individuals to only those with the phenotype defined
            phen_tb_path = 'gs://ukb31063/ukb31063.PHESANT_January_2019.both_sexes.tsv.bgz'
            phen_tb_all = hl.import_table(phen_tb_path,missing='',impute=True,types={'s': hl.tstr}, key='s')
            phen_tb = phen_tb_all.select(phen)
            mt1 = mt0.annotate_cols(phen_str = hl.str(phen_tb[mt0.s][phen]))
            mt1 = mt1.filter_cols(mt1.phen_str == '', keep=False)
            if phen_tb[phen].dtype == hl.dtype('bool'):
                mt1 = mt1.annotate_cols(phen = hl.bool(mt1.phen_str)).drop('phen_str')
            else:
                mt1 = mt1.annotate_cols(phen = hl.float64(mt1.phen_str)).drop('phen_str')     
                
            # remove withdrawn samples
            withdrawn = hl.import_table('gs://nbaya/w31063_20181016.csv',missing='',no_header=True).key_by('f0')
            mt1 = mt1.filter_cols(hl.is_defined(withdrawn[mt1.s]),keep=False) #Remove withdrawn samples
    
            # get gwas of real phenotype
#            real_gwas_path = wd+f'data/{phen}.gwas.qc_pos_variants.ht'
#            try:
#                subprocess.check_output(
#                    [
#                        f'gsutil',
#                        'ls',
#                        real_gwas_path]) is not None
#            except BaseException:
#                print('\n###################')
#                print(f'Starting GWAS for real phenotype {phen}')
#                print('###################\n')
#            
#                cov_list = ['isFemale', 'age', 'age_squared', 'age_isFemale',
#                            'age_squared_isFemale'] + ['PC{:}'.format(i) for i in range(1, 21)]
#                cov_list = list(
#                    map(lambda x: mt1[x] if isinstance(x, str) else x, cov_list))
#                cov_list += [1]
#            
#                gwas_ht = hl.linear_regression_rows(y=mt1.phen,
#                                                    x=mt1.dosage,
#                                                    covariates=cov_list)
#                mt2 = mt1.annotate_rows(alt_AF = hl.agg.mean(mt1.dosage)/2)
#                mt2 = mt2.annotate_rows(minor_AF = hl.min(mt2.alt_AF, 1-mt2.alt_AF))
#                mt2 = mt2.annotate_rows(minor_allele = mt2.alleles[hl.int(mt2.minor_AF==mt2.alt_AF)])
#                mt2_rows = mt2.rows()
#                gwas_ht = gwas_ht.annotate(alt_AF = mt2_rows[gwas_ht.locus, gwas_ht.alleles].alt_AF,
#                                                minor_AF = mt2_rows[gwas_ht.locus, gwas_ht.alleles].minor_AF,
#                                                minor_allele = mt2_rows[gwas_ht.locus, gwas_ht.alleles].minor_allele)
#                ss0 = hl.import_table(wd+f'data/{phen}.gwas.imputed_v3.both_sexes.coding.tsv.bgz',impute=True)
#                ss0 = ss0.annotate(variant = hl.parse_variant(ss0.variant,reference_genome='GRCh37'))
##                ss0 = ss0.rename({'beta':'real_betahat','se':'real_se','tstat':'real_tstat','pval':'real_pval'})
#                ss0 = ss0.annotate(locus = ss0.variant.locus,
#                                   alleles = ss0.variant.alleles)
#                ss0 = ss0.key_by(ss0.locus,ss0.alleles)
#                gwas_ht = gwas_ht.annotate(low_confidence_variant = ss0[gwas_ht.locus, gwas_ht.alleles].low_confidence_variant,
#                                           coding = ss0[gwas_ht.locus, gwas_ht.alleles].coding,
#                                           )
#                gwas_ht.write(real_gwas_path, overwrite=True)
                
#            ss0 = hl.read_table(real_gwas_path).to_pandas()
#            ss0 = ss0.rename(columns={'p_value':'pval'})
            
            ss0 = hl.import_table(wd+f'data/{phen}.gwas.imputed_v3.both_sexes.coding.tsv.bgz',impute=True).to_pandas()
            
            
            # add chromosome and position fields
            if 'chr' not in ss0.columns.values:
                if 'variant' in ss0.columns.values:
                    ss0['chr'] = ss0.variant.str.split(':',n=1,expand=True).iloc[:,0]
                elif 'locus.contig' in ss0.columns.values:
                    ss0 = ss0.rename(columns={'locus.contig':'chr'})
            ss0['chr'] = ss0['chr'].astype(str)
            if 'pos' not in ss0.columns.values:
                if 'variant' in ss0.columns.values:
                    ss0['pos'] = ss0.variant.str.split(':',n=2,expand=True).iloc[:,1]
                elif 'locus.position' in ss0.columns.values:
                    ss0 = ss0.rename(columns={'locus.position':'pos'})
            ss0['pos'] = ss0['pos'].astype(int)
    
            if 'n_complete_samples' in ss0.columns.values:
                print('renaming field n_complete_samples as n')
                ss0 = ss0.rename(columns={'n_complete_samples':'n'})
            if ss0.n.std() != 0:
                print('WARNING: Number of samples varies across SNPs')
    
            if 'EAF' not in ss0.columns.values:
                if 'variant' in ss0.columns.values:
                    ss0['ref'] = ss0.variant.str.split(':',n=3,expand=True).iloc[:,2]
                elif 'alleles' in ss0.columns.values:
                    ss0['ref'] = ss0.alleles.str[0]
                ss0.loc[ss0.minor_allele!=ss0.ref,'alt_af'] = ss0.loc[ss0.minor_allele!=ss0.ref,'minor_AF']
                ss0.loc[ss0.minor_allele==ss0.ref,'alt_af'] = 1-ss0.loc[ss0.minor_allele==ss0.ref,'minor_AF']
                ss0.loc[ss0.beta>0,'raf'] = ss0.loc[ss0.beta>0,'alt_af']
                ss0.loc[ss0.beta<0,'raf'] = 1-ss0.loc[ss0.beta<0,'alt_af']
            else:
                ss0 = ss0.rename(columns={'EAF':'raf'}) # where raf is "risk allele frequency"
    
            ss0.loc[ss0.index,'rbeta'] = np.abs(ss0['beta']) #beta is transformed to risk (or trait-increasing) allele effect size
    
            ss0['var_exp'] = 2*ss0.raf*(1-ss0.raf)*ss0.rbeta**2
    
            ss = ss0[((ss0.pval < pval_threshold)&
                      (~ss0.low_confidence_variant))].reset_index()
    
            if block_mhc: #combine all variants in MHC (restrict to only one hit)
                start = 30400000 # from plot_smiles.py and cytoBand.txt
                stop = 46200000 # from plot_smiles.py and cytoBand.txt
                mhc = ss[(ss.chr=='6')&(ss.pos>=start)&(ss.pos<=stop)]
                if len(mhc)>0: #if there are variants in MHC
                    print(f'Number of variants in MHC (start=6:{start},stop=6:{stop}): {len(mhc)}')
                    non_mhc = ss[~((ss.chr=='6')&(ss.pos>=start)&(ss.pos<=stop))]
                    ss = non_mhc.append(mhc[mhc.pval==mhc.pval.min()])
                print(f'\n... Post-filter # of variants keeping only one variant in MHC: {ss.shape[0]} ...')
    
            # get top variants and top loci (all SNPs within ld_window around sentinel variant)
            top_variants = []
            ct = 0
            if get_top_loci:
                ss_tmp = ss.copy() #dataframe from which top loci will be removed
                ss_keep = ss[[False]*len(ss)]
                i = 0
                if not random_betas:
                    print(f'\n... Filtering to top {n_top_loci} loci with window of {ld_window/1e3} kb ...')
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
                elif random_betas:
                    print(f'\n... Randomly creating top {n_top_loci} loci with window of {ld_window/1e3} kb ...\n')
                    while len(ss_tmp)>0 and i<n_top_loci:
                        ch, pos = np.random.choice(ss_tmp['chr'].astype(str)+':'+ss_tmp['pos'].astype(str)).split(':')
                        pos=int(pos)
                        top_variants.append(f'{ch}:{pos}')
                        ss_w_top_locus = ss_tmp[(ss_tmp['chr'].astype(str)==str(ch))&(ss_tmp['pos']>=pos-ld_window/2)&(ss_tmp['pos']<=pos+ld_window/2)].copy() #extract rows around most significant hit
                        ss_w_top_locus['loci_rank'] = i
                        print(f'locus {i} ({ch}:{pos}) n variants: {len(ss_w_top_locus)}')
                        ct += len(ss_w_top_locus)
                        ss_keep = ss_keep.append(ss_w_top_locus,sort=True)
                        ss_tmp = ss_tmp[~((ss_tmp['chr']==ch)&(ss_tmp['pos']>=pos-ld_window/2)&(ss_tmp['pos']<=pos+ld_window/2))].copy() #keep rows not around most significant hit
                        i += 1
                        if len(ss_tmp)==0 and i<n_top_loci:
                            print(f'\n############\nWARNING: Ran out of SNPs after reading through top {i} loci\n############')
                ss = ss_keep
                ss = ss.sort_values(by='index') #unnecessary step, but restores the order of variants
                
#            ht = hl.Table.from_pandas(ss0)
#            hl.
    
            print(f'\ntotal number of SNPs in top {n_top_loci} loci for {phen} (block_mhc={block_mhc}): {ct}')
            mt2 = mt1.annotate_rows(in_top_loci = 0) #-1 for loci rank indicates that SNP is not in a locus
            for idx, variant in enumerate(top_variants): #WARNING: Loci rank is not a perfect indicator of which top locus a SNP was in if the locus windows are overlapping.
                ch, pos = variant.split(':')
                mt2 = mt2.annotate_rows(in_top_loci = (hl.max(mt2.in_top_loci, ((hl.str(mt2.locus.contig)==ch)*
                                        (hl.int(mt2.locus.position)>=int(pos)-ld_window/2)*
                                        (hl.int(mt2.locus.position)<=int(pos)+ld_window/2)))))
#            mt2 = mt2.annotate_rows(loci_rank = mt2.loci_rank-1)    
#            mt_all_loci = mt2.filter_rows(mt2.loci_rank>=0)  #loci rank is zero indexed 
            mt_all_loci = mt2.filter_rows(mt2.in_top_loci==1)  
            
            count = mt_all_loci.count()
            
            print(f'\ntotal number of SNPs in top {n_top_loci} loci for {phen} in matrix table (block_mhc={block_mhc}): {count}')
            
            
            assert count[0] >= ct
            
            mt_all_loci = mt_all_loci.annotate_globals(**{f'{phen}_top{n_top_loci}loci': top_variants}) #add list of top loci (sentinel variants) to globals
            
            
            ss0 = hl.import_table(wd+f'data/{phen}.gwas.imputed_v3.both_sexes.coding.tsv.bgz',impute=True)
#            ss0 = hl.read_table(real_gwas_path)

            pre = ss0.count()
            print(f'\npre-low confidence filter count: {pre}')
            ss0 = ss0.filter(~ss0.low_confidence_variant)
            post = ss0.count()
            print(f'\npost-low confidence filter count: {post}')
            
#            ss0 = ss0.rename({'beta':'real_betahat','standard_error':'real_se','t_stat':'real_tstat','pval':'real_pval'})
            ss0 = ss0.rename({'beta':'real_betahat','se':'real_se','tstat':'real_tstat','pval':'real_pval'})
#            if 'variant' in ss0.columns.values
            ss0 = ss0.annotate(variant = hl.parse_variant(ss0.variant,reference_genome='GRCh37'))
            ss0 = ss0.annotate(locus = ss0.variant.locus,
                            alleles = ss0.variant.alleles)
            ss0 = ss0.key_by(ss0.locus,ss0.alleles)
    
            ss1 = ss0.annotate(ref = ss0.alleles[0])
            ss1 = ss1.annotate(alt_AF = (ss1.minor_AF*(ss1.minor_allele!=ss1.ref)+
                                         (1-ss1.minor_AF)*(ss1.minor_allele==ss1.ref)))
            ss1 = ss1.annotate(real_raf = ss1.alt_AF*(ss1.real_betahat>0)+(1-ss1.alt_AF)*(ss1.real_betahat<0), #risk allele frequency
                               real_rbeta = hl.abs(ss1.real_betahat)) #beta of risk allele
            ss1 = ss1.annotate(real_varexp = 2*ss1.real_raf*(1-ss1.real_raf)*ss1.real_rbeta**2)
    
            # set only sentinel variants in top hits to have beta = real_betahat, otherwise = 0
            mt_all_loci = mt_all_loci.annotate_rows(sim_truebeta = (hl.literal(top_variants).contains(hl.str(mt_all_loci.locus))*
                                                                    ss1[mt_all_loci.locus,mt_all_loci.alleles].real_betahat),
                                                    real_rbeta = ss1[mt_all_loci.locus,mt_all_loci.alleles].real_rbeta,
                                                    real_raf = ss1[mt_all_loci.locus,mt_all_loci.alleles].real_raf,
                                                    real_varexp = ss1[mt_all_loci.locus,mt_all_loci.alleles].real_varexp,
                                                    real_pval = ss1[mt_all_loci.locus,mt_all_loci.alleles].real_pval,
                                                    minor_AF = ss1[mt_all_loci.locus,mt_all_loci.alleles].minor_AF,
                                                    minor_allele = ss1[mt_all_loci.locus,mt_all_loci.alleles].minor_allele,
                                                    alt_af = ss1[mt_all_loci.locus,mt_all_loci.alleles].alt_AF,
                                                    coding = ss1[mt_all_loci.locus,mt_all_loci.alleles].coding)
            
            mt_all_loci = mt_all_loci.filter_rows(hl.is_defined(ss1[mt_all_loci.locus,mt_all_loci.alleles]))
            
            count = mt_all_loci.count()
            
            print('\n####################')
            print(f'Phenotype: {phen}')
            print(f'Dimensions of matrix table of top {n_top_loci} loci: {count}')
            print(f'pval_threshold: {pval_threshold}')
            print(f'ld_window: {ld_window}')
            print(f'maf: {maf}')
            print(f'block_mhc: {block_mhc}')
            print('####################\n')
    
            mt_all_loci.write(mt_all_loci_path,overwrite=True)

        mt_all_loci = hl.read_matrix_table(mt_all_loci_path)
        
#        if get_top_loci and n_top_loci_sim!=n_top_loci:
#            if 'loci_rank' not in mt_all_loci.row:
#                top_variants = hl.eval(mt_all_loci[f'{phen}_top{n_top_loci}loci'])
#                top_variants = top_variants[:n_top_loci_sim]
#                mt_all_loci = mt_all_loci.annotate_rows(loci_rank = -1) #-1 for loci rank indicates that SNP is not in a locus
#                for idx, variant in enumerate(top_variants):
#                    ch, pos = variant.split(':')
#                    mt_all_loci = mt_all_loci.annotate_rows(loci_rank = idx*(hl.max(mt_all_loci.loci_rank, (hl.str(mt_all_loci.locus.contig)==ch)*
#                                            (hl.int(mt_all_loci.locus.position)>=int(pos)-ld_window/2)*
#                                            (hl.int(mt_all_loci.locus.position)<=int(pos)+ld_window/2))))
#            mt_all_loci = mt_all_loci.filter_rows((mt_all_loci.loci_rank>=0)&(mt_all_loci.loci_rank<n_top_loci_sim)) #loci rank is zero indexed
#            assert hl.agg.max(mt_all_loci.loci_rank)==n_top_loci_sim-1
        
#        h2 = mt_all_loci.aggregate_rows(hl.agg.sum(mt_all_loci.real_varexp*(mt_all_loci.sim_truebeta!=0))) #total variance explained by n_top_loci_sim sentinel variants
        h2 = 1
        


        count = mt_all_loci.count()
        print('\n####################')
        print(f'Phenotype: {phen}')
        print(f'Dimensions of matrix table of top {n_top_loci} loci: {count}')
        print(f'pval_threshold: {pval_threshold}')
        print(f'ld_window: {ld_window}')
        print(f'maf: {maf}')
        print(f'block_mhc: {block_mhc}')
        print(f'n_top_loci_sim: {n_top_loci_sim}')
        print(f'h2: {h2}')
        print('####################\n')
        sim_all_loci = calculate_phenotypes(
            mt=mt_all_loci,
            genotype=mt_all_loci.dosage,
            beta=mt_all_loci.sim_truebeta,
            h2=h2,
            exact_h2=False)

        sim_all_loci = sim_all_loci.checkpoint(sim_all_loci_path,overwrite=True)

    return sim_all_loci


def gwas_all_loci(sim_all_loci_mt, phen, n_top_loci, ld_window, block_mhc, random_betas=False):
    
#    withdrawn = hl.import_table('gs://nbaya/w31063_20181016.csv',missing='',no_header=True).key_by('f0')
#    sim_all_loci_mt = sim_all_loci_mt.filter_cols(hl.is_defined(withdrawn[sim_all_loci_mt.s]),keep=False) #Remove withdrawn samples
    
    suffix = f'.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb{".block_mhc" if block_mhc else ""}{".random_betas" if random_betas else ""}'
    try:
        subprocess.check_output(
            [
                f'gsutil',
                'ls',
                wd +
                phen +
                f'.gwas{suffix}.ht/_SUCCESS']) is not None
    except BaseException:
        print('\n###################')
        print(f'Starting GWAS for simulated phenotype {phen}')
        print('###################\n')
    
        cov_list = ['isFemale', 'age', 'age_squared', 'age_isFemale',
                    'age_squared_isFemale'] + ['PC{:}'.format(i) for i in range(1, 21)]
        cov_list = list(
            map(lambda x: sim_all_loci_mt[x] if isinstance(x, str) else x, cov_list))
        cov_list += [1]
    
        gwas_all_loci_ht = hl.linear_regression_rows(y=sim_all_loci_mt.y,
                                                     x=sim_all_loci_mt.dosage,
                                                     covariates=cov_list,
                                                     pass_through=['sim_truebeta', 'real_rbeta',
                                                                   'real_raf', 'real_pval',
                                                                   'coding', 'real_varexp',
                                                                   'minor_AF', 'minor_allele',
                                                                   'alt_af'])
        gwas_all_loci_ht.write(wd + phen +f'.gwas{suffix}.ht',overwrite=True)


def get_loci_ld(phen, n_top_loci, ld_window, block_mhc, random_betas=False):
    suffix = f'.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb{".block_mhc" if block_mhc else ""}{".random_betas" if random_betas else ""}'
    sim = hl.read_matrix_table(
                wd +
                phen +
                f'.simulation{suffix}.mt')
    ht = hl.read_table(
                wd +
                phen +
                f'.gwas{suffix}.ht')
    df = ht.to_pandas()
    ss = df.copy()
    ss = ss.rename(columns={'locus.contig':'chr','locus.position':'pos',
                           'real_pval':'pval'})
#    if block_mhc: #combine all variants in MHC (restrict to only one hit)
#        start = 30400000 # from plot_smiles.py and cytoBand.txt
#        stop = 46200000 # from plot_smiles.py and cytoBand.txt
#        mhc = ss[(ss.chr=='6')&(ss.pos>=start)&(ss.pos<=stop)]
#        if len(mhc)>0: #if there are variants in MHC
#            print(f'Number of variants in MHC (start=6:{start},stop=6:{stop}): {len(mhc)}')
#            non_mhc = ss[~((ss.chr=='6')&(ss.pos>=start)&(ss.pos<=stop))]
#            ss = non_mhc.append(mhc[mhc.pval==mhc.pval.min()])
    print(f'\n... Post-filter # of variants keeping a max of one variant in MHC: {ss.shape[0]} ...')
          
    top_variants = []
    ct = 0
    
    ss_tmp = ss.copy() #dataframe from which top loci will be removed
    ss_keep = ss[[False]*len(ss)]
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
    for s_id, sentinel in enumerate(top_variants):
        print(f'\n############\nStarting sentinel {sentinel} ({s_id+1} of {len(top_variants)})\n############\n')
        chrom, pos = sentinel.split(':')
        mt_sentinel = sim.filter_rows(sim.locus==hl.parse_locus(sentinel))
#        gt_at_sentinel = mt_sentinel.dosage.collect()
#        mt1 = sim.add_col_index('tmp_col_idx')
#        mt1 = mt1.annotate_cols(gt_at_sentinel = hl.literal(gt_at_sentinel)[hl.int32(mt1.tmp_col_idx)])
        mt_sentinel = mt_sentinel.annotate_cols(gt_at_sentinel = hl.agg.collect(mt_sentinel.dosage)[0])
        mt1 = sim.annotate_cols(gt_at_sentinel = mt_sentinel.cols()[sim.s].gt_at_sentinel)
        mt_locus = hl.filter_intervals(mt1, [hl.parse_locus_interval(f'{chrom}:{int(pos)-int(ld_window/2)}-{int(pos)+int(ld_window/2)}')])
        mt_locus_row_count = mt_locus.count_rows()
        print(f'\n############\nRows in mt_locus: {mt_locus_row_count}\n############\n')
        mt_locus = mt_locus.annotate_rows(r = hl.agg.corr(mt_locus.gt_at_sentinel, mt_locus.dosage))
        ht = ht.annotate(**{f'r_locus{s_id}': mt_locus.rows()[ht.locus, ht.alleles].r})
        
    print('\nStarting to write out table...\n')
    
    ht.export(wd +
              phen +
              f'.gwas.corr{suffix}.tsv.gz')
        
    
if __name__ == '__main__':
    get_top_loci = True
    n_top_loci = 20
    n_top_loci_sim = 20
    ld_window = int(300e3)
    pval_threshold = 1
    maf = 0
    block_mhc = True
    random_betas=False
    
    phen = args.phen
    if phen == '50':
        phen='50_irnt'
    if phen=='21001':
        phen='21001_irnt'
        
    sim_all_loci_mt = get_sim_all_loci_mt(phen=phen,
                                          n_top_loci=n_top_loci,
                                          ld_window=ld_window,
                                          pval_threshold=pval_threshold,
                                          maf=maf,
                                          block_mhc=block_mhc,
                                          random_betas=random_betas,
                                          n_top_loci_sim=n_top_loci_sim)
    
    gwas_all_loci(sim_all_loci_mt=sim_all_loci_mt,
                  phen=phen,
                  n_top_loci=n_top_loci,
                  ld_window=ld_window,
                  block_mhc=block_mhc,
                  random_betas=random_betas)
    
    get_loci_ld(phen=phen, 
                n_top_loci=n_top_loci, 
                ld_window=ld_window, 
                block_mhc=block_mhc,
                random_betas=random_betas)