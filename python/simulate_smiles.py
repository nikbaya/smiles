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
import requests
import subprocess
import argparse
#from hail.experimental.ldscsim import calculate_phenotypes
url = 'https://raw.githubusercontent.com/nikbaya/ldscsim/master/ldscsim.py'
r = requests.get(url).text
exec(r)
hl.init(log='/tmp/foo.log')


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--phen', type=str, required=True, help="phenotype")

args = parser.parse_args()

wd = 'gs://nbaya/smiles/'


def get_sim_all_loci_mt(phen, n_top_loci, ld_window, pval_threshold, maf):
    try:
        subprocess.check_output(
            [
                f'gsutil',
                'ls',
                wd +
                phen +
                f'.simulation.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb.mt/_SUCCESS']) is not None
        sim_all_loci = hl.read_matrix_table(
            wd + phen + f'.simulation.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb.mt')

    except BaseException:
        print('Getting sim_all_loci matrix table...')
#        mt0 = hl.read_matrix_table(f'gs://nbaya/split/ukb31063.qc_pos_variants.gwas_samples_repart.mt')
#        withdrawn = hl.import_table('gs://nbaya/w31063_20181016.csv',missing='',no_header=True).key_by('f0')
#        mt0 = mt0.filter_cols(hl.is_defined(withdrawn[mt0.s]),keep=False) #Remove withdrawn samples
#
#        phen_tb_path = 'gs://ukb31063/ukb31063.PHESANT_January_2019.both_sexes.tsv.bgz'
#        phen_tb_all = hl.import_table(phen_tb_path,missing='',impute=True,types={'s': hl.tstr}, key='s')
#        phen_tb = phen_tb_all.select(phen)
#        mt1 = mt0.annotate_cols(phen_str = hl.str(phen_tb[mt0.s][phen]))
#        mt1 = mt1.filter_cols(mt1.phen_str == '', keep=False)
#
#        ss0 = hl.import_table(wd+f'data/{phen}.gwas.imputed_v3.both_sexes.coding.tsv.bgz',impute=True).to_pandas()
#        if 'chr' not in ss0.columns.values:
#            ss0['chr'] = ss0.variant.str.split(':',n=1,expand=True).iloc[:,0]
#        ss0['chr'] = ss0['chr'].astype(str)
#
#        if 'pos' not in ss0.columns.values:
#            ss0['pos'] = ss0.variant.str.split(':',n=2,expand=True).iloc[:,1]
#        ss0['pos'] = ss0['pos'].astype(int)
#
#        if 'n_complete_samples' in ss0.columns.values:
#            print('renaming field n_complete_samples as n')
#            ss0 = ss0.rename(columns={'n_complete_samples':'n'})
##        n = int(ss0.n.mean())
#        if ss0.n.std() != 0:
#            print('WARNING: Number of samples varies across SNPs')
#
#        if 'EAF' not in ss0.columns.values:
#            ss0['ref'] = ss0.variant.str.split(':',n=3,expand=True).iloc[:,2]
#            ss0.loc[ss0.minor_allele!=ss0.ref,'alt_af'] = ss0.loc[ss0.minor_allele!=ss0.ref,'minor_AF']
#            ss0.loc[ss0.minor_allele==ss0.ref,'alt_af'] = 1-ss0.loc[ss0.minor_allele==ss0.ref,'minor_AF']
#            ss0.loc[ss0.beta>0,'raf'] = ss0.loc[ss0.beta>0,'alt_af']
#            ss0.loc[ss0.beta<0,'raf'] = 1-ss0.loc[ss0.beta<0,'alt_af']
#        else:
#            ss0 = ss0.rename(columns={'EAF':'raf'}) # where raf is "risk allele frequency"
#
#        ss0.loc[ss0.index,'rbeta'] = np.abs(ss0['beta']) #beta is transformed to risk (or trait-increasing) allele effect size
#
#        ss0['var_exp'] = 2*ss0.raf*(1-ss0.raf)*ss0.rbeta**2
#
#        ss = ss0[(ss0.pval < pval_threshold)&(~ss0.low_confidence_variant)].reset_index()
#
#
#        top_variants = []
#
#        ct = 0
#
#        if get_top_loci:
# ss_ls = [None]*n_top_loci #list of dataframes with top loci
#            ss_tmp = ss.copy() #dataframe from which top loci will be removed
#            ss_keep = ss[[False]*len(ss)]
#            i = 0
#            print(f'filtering to top {n_top_loci} loci with window of {ld_window/1e3} kb')
#            while len(ss_tmp)>0 and i<n_top_loci:
#                ch, pos = ss_tmp[ss_tmp.pval==ss_tmp.pval.min()][['chr','pos']].values[0]
#                top_variants.append(f'{ch}:{pos}')
#                ss_w_top_locus = ss_tmp[(ss_tmp['chr'].astype(str)==str(ch))&(ss_tmp['pos']>=pos-ld_window/2)&(ss_tmp['pos']<=pos+ld_window/2)].copy() #extract rows around most significant hit
#                ss_w_top_locus['loci_rank'] = i
#                print(f'locus {i} ({ch}:{pos}) n variants: {len(ss_w_top_locus)}')
#                ct += len(ss_w_top_locus)
#                ss_keep = ss_keep.append(ss_w_top_locus)
#                ss_tmp = ss_tmp[~((ss_tmp['chr']==ch)&(ss_tmp['pos']>=pos-ld_window/2)&(ss_tmp['pos']<=pos+ld_window/2))].copy() #keep rows not around most significant hit
#                i += 1
#                if len(ss_tmp)==0 and i<n_top_loci:
#                    print(f'\n############\nWARNING: Ran out of SNPs after reading through top {i} loci\n############')
#            ss = ss_keep
#            ss = ss.sort_values(by='index') #unnecessary step, but restores the order of variants
#
#        print(f'total number of SNPs in top {n_top_loci} loci for {phen}: {ct}')
#
#        mt2 = mt1.annotate_rows(in_top_loci = 0)
#        for variant in top_variants:
#            ch, pos = variant.split(':')
#            mt2 = mt2.annotate_rows(in_top_loci = hl.max(mt2.in_top_loci, (hl.str(mt2.locus.contig)==ch)*
#                                    (hl.int(mt2.locus.position)>=int(pos)-ld_window/2)*
#                                    (hl.int(mt2.locus.position)<=int(pos)+ld_window/2)))
#        mt_top_loci = mt2.filter_rows(mt2.in_top_loci==1)
#
#
#        ss0 = hl.import_table(wd+f'data/{phen}.gwas.imputed_v3.both_sexes.coding.tsv.bgz',impute=True)
#        pre = ss0.count()
#        print(f'\npre-filter count: {pre}')
#        ss0 = ss0.filter(~ss0.low_confidence_variant)
#        post = ss0.count()
#        print(f'\npost-filter count: {post}')
#        ss0 = ss0.annotate(variant = hl.parse_variant(ss0.variant,reference_genome='GRCh37'))
#        ss0 = ss0.rename({'beta':'real_betahat','se':'real_se','tstat':'real_tstat','pval':'real_pval'})
#        ss0 = ss0.annotate(locus = ss0.variant.locus,
#                        alleles = ss0.variant.alleles)
#        ss0 = ss0.key_by(ss0.locus,ss0.alleles)
#
#        ss1 = ss0.annotate(ref = ss0.alleles[0])
#        ss1 = ss1.annotate(alt_AF = (ss1.minor_AF*(ss1.minor_allele!=ss1.ref)+
#                                     (1-ss1.minor_AF)*(ss1.minor_allele==ss1.ref)))
#        ss1 = ss1.annotate(real_raf = ss1.alt_AF*(ss1.real_betahat>0)+(1-ss1.alt_AF)*(ss1.real_betahat<0), #risk allele frequency
#                           real_rbeta = hl.abs(ss1.real_betahat)) #beta of risk allele
#        ss1 = ss1.annotate(real_varexp = 2*ss1.real_raf*(1-ss1.real_raf)*ss1.real_rbeta**2)
#
#        # set only top hits to have beta = 1, otherwise = 0
#        mt_all_loci = mt_top_loci.annotate_rows(sim_truebeta = hl.literal(top_variants).contains(hl.str(mt_top_loci.locus)),
#                                                real_rbeta = ss1[mt_top_loci.locus,mt_top_loci.alleles].real_rbeta,
#                                                real_raf = ss1[mt_top_loci.locus,mt_top_loci.alleles].real_raf,
#                                                real_pval = ss1[mt_top_loci.locus,mt_top_loci.alleles].real_pval,
#                                                minor_AF = ss1[mt_top_loci.locus,mt_top_loci.alleles].minor_AF,
#                                                minor_allele = ss1[mt_top_loci.locus,mt_top_loci.alleles].minor_allele,
#                                                alt_af = ss1[mt_top_loci.locus,mt_top_loci.alleles].alt_AF,
#                                                coding = ss1[mt_top_loci.locus,mt_top_loci.alleles].coding,
#                                                real_varexp = ss1[mt_top_loci.locus,mt_top_loci.alleles].real_varexp)
#        mt_all_loci = mt_all_loci.filter_rows(hl.is_defined(ss1[mt_all_loci.locus,mt_all_loci.alleles]))
#
#        count = mt_all_loci.count()
#        print('\n####################')
#        print(f'Phenotype: {phen}')
#        print(f'Dimensions of matrix table of top {n_top_loci} loci: {count}')
#        print(f'pval_threshold: {pval_threshold}')
#        print(f'ld_window: {ld_window}')
#        print(f'maf: {maf}')
#        print('####################\n')
#
#        mt_all_loci = mt_all_loci.checkpoint(wd+phen+f'.genotypes.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb.mt',overwrite=True)

        mt_all_loci = hl.read_matrix_table(
            wd + phen + f'.genotypes.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb.mt')

# sim_all_loci = mt_all_loci.annotate_cols(y =
# hl.agg.sum(mt_all_loci.sim_beta * mt_all_loci.dosage)) #simulate
# phenotype without noise (no need to normalize genotypes)
        sim_all_loci = calculate_phenotypes(
            mt=mt_all_loci,
            genotype=mt_all_loci.dosage,
            beta=mt_all_loci.sim_truebeta,
            h2=1)

        sim_all_loci = sim_all_loci.checkpoint(
            wd +
            phen +
            f'.simulation.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb.mt',
            overwrite=True)

        return sim_all_loci


def add_real_fields(phen):
    ss0 = hl.import_table(
        wd +
        f'data/{phen}.gwas.imputed_v3.both_sexes.coding.tsv.bgz',
        impute=True)
    ss0 = ss0.annotate(
        variant=hl.parse_variant(
            ss0.variant,
            reference_genome='GRCh37'))
    ss0 = ss0.rename({'beta': 'real_betahat', 'se': 'real_se',
                      'tstat': 'real_tstat', 'pval': 'real_pval'})
    ss0 = ss0.annotate(locus=ss0.variant.locus,
                       alleles=ss0.variant.alleles)
    ss0 = ss0.key_by(ss0.locus, ss0.alleles)

    ss1 = ss0.annotate(ref=ss0.alleles[0])
    ss1 = ss1.annotate(alt_AF=(ss1.minor_AF * (ss1.minor_allele != ss1.ref) +
                               (1 - ss1.minor_AF) * (ss1.minor_allele == ss1.ref)))
    ss1 = ss1.annotate(real_raf=ss1.alt_AF * (ss1.real_betahat > 0) + (1 - ss1.alt_AF) * (ss1.real_betahat < 0),  # risk allele frequency
                       real_rbeta=hl.abs(ss1.real_betahat))  # beta of risk allele
    ss1 = ss1.annotate(real_varexp=2 * ss1.real_raf *
                       (1 - ss1.real_raf) * ss1.real_rbeta**2)


def gwas_all_loci(sim_all_loci_mt, phen, n_top_loci, ld_window):

    print('\n###################')
    print(f'Starting GWAS for phenotype {phen}')
    print('###################\n')

    cov_list = ['isFemale', 'age', 'age_squared', 'age_isFemale',
                'age_squared_isFemale'] + ['PC{:}'.format(i) for i in range(1, 21)]
    cov_list = list(
        map(lambda x: sim_all_loci_mt[x] if isinstance(x, str) else x, cov_list))
    cov_list += [1]

    gwas_all_loci_ht = hl.linear_regression_rows(y=sim_all_loci_mt.y,
                                                 x=sim_all_loci_mt.dosage,
                                                 covariates=cov_list,
                                                 pass_through=['sim_beta', 'real_rbeta',
                                                               'real_raf', 'real_pval',
                                                               'coding', 'real_varexp',
                                                               'minor_AF', 'minor_allele',
                                                               'alt_AF'])
    gwas_all_loci_ht.write(
        wd +
        phen +
        f'.gwas.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb.ht',
        overwrite=True)


def get_loci_ld(phen, n_top_loci, ldwindow):
    mt = hl.read_matrix_table(
        wd +
        phen +
        f'.simulation.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb.mt')
    locus_ld = hl.row_correlation(mt.dosage)
    locus_ld.write(
        wd +
        phen +
        f'.r.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb.bm')


if __name__ == '__main__':
    get_top_loci = True
    n_top_loci = 100
    ld_window = 300e3
    pval_threshold = 1
    maf = 0
    phen = args.phen

    sim_all_loci_mt = get_sim_all_loci_mt(phen=phen,
                                          n_top_loci=n_top_loci,
                                          ld_window=ld_window,
                                          pval_threshold=pval_threshold,
                                          maf=maf)
    gwas_all_loci(sim_all_loci_mt=sim_all_loci_mt,
                  phen=phen,
                  n_top_loci=n_top_loci,
                  ld_window=ld_window)
