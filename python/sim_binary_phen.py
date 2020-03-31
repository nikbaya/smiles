#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 18:35:03 2020

For simulating a binary phenotype upon which to test whether variance explained is monotonic in p-value

@author: nbaya
"""

import argparse
import hail as hl
from hail.experimental import ldscsim
hl.init(log='/tmp/hail.log')


def get_mt(remove_withdrawn=True):
    ## matrix table of 1008898 HM3 variants and white British subset of UKB
    mt = hl.read_matrix_table('gs://nbaya/split/ukb31063.hm3_variants.gwas_samples_repart.mt')    

    if remove_withdrawn:    
        withdrawn = hl.import_table('gs://nbaya/w31063_20200204.csv', missing='', no_header=True)
        withdrawn = withdrawn.rename({'f0': 's'})  # rename field with sample IDs
        withdrawn = withdrawn.key_by('s')
        mt = mt.filter_cols(hl.is_defined(withdrawn[mt.s]), keep=False)
        
    return mt


def get_sim_mt(mt, h2, pi, K):
    r'''
    Simulate phenotype on samples in `mt` with h2=`h2`, probability of being 
    causal `pi`, population prevalence `K`. Uses `gt_field` entry field as the 
    genotypes in the simulation.
    '''
    assert {'GT','dosage'}.intersection(mt.entry)!={}, "mt does not have an entry field named 'dosage' or 'GT' corresponding to genotype data"
    gt_field = mt.dosage if 'dosage' in mt.entry else mt.GT #no need to add .n_alt_alleles(), simulate_phenotypes can take call expressions as genotypes
    sim_mt = ldscsim.simulate_phenotypes(mt=mt, genotype=gt_field, h2=h2, pi=pi)
    sim_mt = ldscsim.binarize(mt=sim_mt, y=sim_mt.y, K=K, exact=False)
    
    return sim_mt


def get_shuffled_ht(ht, phen : str, is_cas : bool, seed=None):
    r'''
    Returns shuffled Table of cases if `is_cas`=True, controls if `is_cas`=False.
    Case status is determined by binary field `phen`.
    '''
    ht = ht.filter(ht[phen]==is_cas)
    
    ht = ht.annotate(tmp_rand = hl.rand_norm(seed=seed))
    ht = ht.order_by('tmp_rand')
    ht = ht.add_index('tmp_idx')
    
    return ht
   
    
def make_subsets(ht, phen : str, n_list : list, n_cas_list : list, seed=None):
    r'''
    Makes `len(n_list)` subsets from `ht`. Subset i has `n_list[i]` samples and
    `n_cas_list[i]` cases. Uses binary field `phen` in `ht` to determine case status. 
    `n_cas_list[i]` is always bounded by `n_list[i]`.
    
    Returns
    -------
    :class:`.Table`
        :class:`.Table` with samples annotated by field "subset_idx", indicating
        which subset they were assigned. (subset_idx=1 is first subset)
    '''
    assert len(n_list)==len(n_cas_list), "lengths of `n` and `n_cas` lists do not match."
    n_total = ht.count()
    assert n_total >= sum(n_list), "not enough samples in `ht` to produce subsets with sample sizes from list `n_list`"
    n_cas_total = ht.filter(ht[phen]==1).count()
    assert n_cas_total >= sum(n_cas_list), "not enough casesin `ht` to produce subsets with case sizes from list `n_cas_list`"
    n_con_list = [n_i - n_cas_i for n_i, n_cas_i in zip(n_list, n_cas_list)]
    
    ht = ht.annotate(subset_idx = 0)
    cas =  get_shuffled_ht(ht=ht, phen=phen, is_cas=True, seed=seed)
    con =  get_shuffled_ht(ht=ht, phen=phen, is_cas=False, seed=seed)
    
    n_cas_prev = 0
    n_con_prev = 0
    for idx in range(len(n_list)):
        n_cas_curr = n_cas_list[idx]
        n_con_curr = n_con_list[idx]
#        cas = cas.annotate(subset_idx = cas.subset_idx+(cas.tmp_idx<n_cas_prev+n_cas_curr)) # reverses order of subset idx's relative to order of n_list and n_cas_list
#        con = con.annotate(subset_idx = con.subset_idx+(con.tmp_idx<n_con_prev+n_con_curr))
        cas = cas.annotate(subset_idx = cas.subset_idx+(cas.tmp_idx<n_cas_prev+n_cas_curr))
        con = con.annotate(subset_idx = con.subset_idx+(con.tmp_idx<n_con_prev+n_con_curr))
        n_cas_prev += n_cas_curr
        n_con_prev += n_con_curr
    
    for idx in range(len(n_list)):
        cas_curr_ct = cas.filter(cas.subset_idx==idx+1).count()
        con_curr_ct = con.filter(con.subset_idx==idx+1).count()
        
        print(f'\n\ncas ct for subset {idx+1}: {cas_curr_ct}')
        print(f'con ct for subset {idx+1}: {con_curr_ct}\n\n')
    
    cas = cas.filter(cas.subset_idx>0)
    con = con.filter(con.subset_idx>0)
    
    ht = cas.union(con)
        
    return ht

def run_gwas(mt, phen : str, sim_name : str,  subset_idx : int, param_suffix : str, wd : str,
             is_logreg=True):
    assert {'GT','dosage'}.intersection(mt.entry)!={}, "mt does not have an entry field named 'dosage' or 'GT' corresponding to genotype data"
    
    mt = mt.filter_cols(mt.subset_idx==subset_idx)
    mt = mt.filter_cols(hl.is_defined(mt[phen]))
    print(f'\n\ngwas sample count (subset {subset_idx}): {mt.count_cols()}\n\n')
    
    if 'dosage' in mt.entry:
        mt = mt.annotate_rows(EAF = hl.agg.mean(mt.dosage)/2)
    elif 'GT' in mt.entry:
        mt = mt.annotate_rows(EAF = hl.agg.mean(mt.GT.n_alt_alleles())/2)
    
    gwas_path = f'{wd}/gwas.{"logreg" if is_logreg else "linreg"}.{sim_name}.subset_{subset_idx}.{param_suffix}.tsv.gz'
    
    if not hl.hadoop_is_file(gwas_path):
        gt_field = mt.dosage if 'dosage' in mt.entry else mt.GT.n_alt_alleles()
        
        if is_logreg:
            gwas_ht = hl.logistic_regression_rows(test='wald', 
                                                  y=mt[phen], 
                                                  x=gt_field,
                                                  covariates=[1],
                                                  pass_through=['EAF'])
        else:
            gwas_ht = hl.linear_regression_rows(y=mt[phen], 
                                                x=gt_field,
                                                covariates=[1],
                                                pass_through=['EAF'])
        gwas_ht.select('EAF','beta','standard_error','p_value').export(gwas_path)
    
    else: 
        print(f'GWAS already run! ({gwas_path})')
        gwas_ht = hl.import_table(gwas_path, impute=True, force=True)
        gwas_ht = gwas_ht.annotate(locus = hl.parse_locus(gwas_ht.locus),
                                   alleles = gwas_ht.alleles.replace('\[\"','').replace('\"\]','').split('\",\"'))
        gwas_ht = gwas_ht.key_by('locus','alleles')

    return gwas_ht
    
def metaanalyze_gwas(subsets, gwas_ht_list, sim_name, param_suffix, wd):

    if len(gwas_ht_list)==1: # if list is single GWAS, don't meta-analyze
        return gwas_ht_list[0]
    
    sample_ct_dict = {}
    
    for subset_idx, tmp_gwas_ht in enumerate(gwas_ht_list,1): 
        sample_ct = subsets.filter(subsets.subset_idx==subset_idx).count()
        sample_ct_dict[subset_idx] = sample_ct
        print(f'\n\nmeta-analysis sample count subset {subset_idx}: {sample_ct}\n\n')
    
    comb_gwas_ht = gwas_ht_list[0].annotate(subset_idx=1,
                                            n = sample_ct_dict[1])
    union_args = [ht.annotate(subset_idx=subset_idx,
                              n = sample_ct_dict[subset_idx]) for subset_idx, ht in enumerate(gwas_ht_list[1:],2)] # list of gwas_ht's to join
    comb_gwas_ht = comb_gwas_ht.union(*union_args)
    
    comb_gwas_ht = comb_gwas_ht.annotate(w = 1/(comb_gwas_ht['standard_error']**2))
    
    agg_expr = {
        'meta_se': hl.sqrt(1/(hl.agg.sum(comb_gwas_ht.w))),
        'meta_beta': hl.agg.sum(comb_gwas_ht['beta'] * comb_gwas_ht.w)/hl.agg.sum(comb_gwas_ht.w),
        'meta_EAF': hl.agg.sum(comb_gwas_ht['EAF']*comb_gwas_ht['n'])/hl.agg.sum(comb_gwas_ht['n'])
    }
    
    comb_gwas_ht = comb_gwas_ht.group_by('locus','alleles').aggregate(**agg_expr)
    
    comb_gwas_ht = comb_gwas_ht.annotate(meta_pval = 2*hl.pnorm(-hl.abs(comb_gwas_ht.meta_beta/
                                                                        comb_gwas_ht.meta_se)))
    
    meta_gwas_path = f'{wd}/gwas.logreg.{sim_name}.{param_suffix}.tsv.gz'
    comb_gwas_ht.export(meta_gwas_path)
    

if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--sim_name', default=None, type=str, help='Which simulation setup to use')
    args = parser.parse_args()
    
    sim_name = args.sim_name
    
    smiles_wd = 'gs://nbaya/smiles'
    
    if sim_name=='meta1':
        seed = 1
        h2=0.6
        pi=0.01
        K=0.05 # population prevalence
        n_cas_list = [3000]*5 
        n_list = [2*x for x in n_cas_list]
    elif sim_name=='meta2':
        seed = 1
        h2=0.6
        pi=0.01
        K=0.05
        n_cas_list = [10000, 2000, 1000, 1000, 500, 500]
        n_list = [2*x for x in n_cas_list]
    elif sim_name=='meta3':
        seed = 1
        h2=0.6
        pi=0.01
        K=0.05
        n_cas_list = [500]*30
        n_list = [2*x for x in n_cas_list]
    elif sim_name=='nonmeta_ascertained':
        seed = 1
        h2=0.6
        pi=0.01
        K=0.05
        n_cas_list = [18000]
        n_list = [2*x for x in n_cas_list]    
    elif sim_name=='nonmeta_ascertained2':
        seed = 1
        h2=0.6
        pi=0.01
        K=0.05
        n_cas_list = [15000]
        n_list = [2*x for x in n_cas_list]  
    elif sim_name=='bn_nonmeta1':
        seed = 1
        h2=0.6
        pi=0.01
        K=0.05
        
        n_pops = 3
        fst = [0.1]*n_pops
        n_vars = int(100e3)

        n_cas_list = [5000]*3
        n_list = [2*x for x in n_cas_list]

        n_sim = int(320e3) # over simulate to be able to have ascertainment
    else:
        raise ValueError(f'sim_name="{sim_name}" does not match any models')
    
    hl.set_global_seed(seed)
    
    gt_sim_suffix = f'bn.npops_{n_pops}.nvars_{n_vars}.nsim_{n_sim}' if sim_name[:3]=='bn_' else '' # suffix for genotype simulation (empty string if using ukb data)
    param_suffix = f'{gt_sim_suffix}.h2_{h2}.pi_{pi}.K_{K}.seed_{seed}'
    betas_path=f'{smiles_wd}/betas.{param_suffix}.tsv.gz'
    phens_path=f'{smiles_wd}/phens.{param_suffix}.tsv.gz'
    
    if sim_name[:3]=='bn_':
        mt = hl.balding_nichols_model(n_populations=n_pops, 
                                      n_samples=n_sim, 
                                      n_variants=n_vars, 
                                      fst=fst)
        
        mt = mt.filter_rows((hl.abs(hl.agg.mean(mt.GT.n_alt_alleles())/2-0.5)<0.5)) # remove invariant SNPs
        mt = mt.annotate_cols(s = hl.str(mt.sample_idx))
    
        
        if hl.hadoop_is_file(betas_path) and hl.hadoop_is_file(phens_path):
#            betas = hl.import_table(betas_path, impute=True, force=True)             
#            betas = betas.annotate(locus = hl.parse_locus(betas.locus),
#                                   alleles = betas.alleles.replace('\[\"','').replace('\"\]','').split('\",\"'))
#            betas = betas.key_by('locus','alleles')

            phens = hl.import_table(phens_path, key=['s'], types={'s':hl.tstr}, impute=True, force=True) 
            
            sim_mt = mt.annotate_cols(y_binarized=phens[mt.s].y_binarized)
            
        else:
            sim_mt = get_sim_mt(mt=mt, 
                                h2=h2, 
                                pi=pi, 
                                K=K)
            
            sim_mt.rows().select('beta').export(betas_path)    
            sim_mt.cols().select('y','y_binarized').export(phens_path)

        # TODO: write out matrix table with sim results? Should always be able to get same exact mt if global seed is set
        
    elif sim_name[:3]!='bn_' and hl.hadoop_is_file(betas_path) and hl.hadoop_is_file(phens_path):
        mt = get_mt(remove_withdrawn=False) # no need to remove withdrawn samples because phenotypes have only been calculated for non-withdrawn samples
        
        betas = hl.import_table(betas_path, impute=True, force=True) 
        phens = hl.import_table(phens_path, key=['s'], types={'s':hl.tstr},impute=True, force=True) 
        
        betas = betas.annotate(locus = hl.parse_locus(betas.locus),
                               alleles = betas.alleles.replace('\[\"','').replace('\"\]','').split('\",\"'))
        betas = betas.key_by('locus','alleles')
        
#        sim_mt = mt.annotate_rows(beta=betas[mt.locus, mt.alleles].beta)
        sim_mt = sim_mt.annotate_cols(y_binarized=phens[sim_mt.s].y_binarized)
        
    else:
        mt = get_mt(remove_withdrawn=True)
        sim_mt = get_sim_mt(mt=mt, 
                            h2=h2, 
                            pi=pi, 
                            K=K)

        sim_mt.rows().select('beta').export(betas_path)    
        sim_mt.cols().select('y','y_binarized').export(phens_path)
    
#    print(f'\n\nphenotype stats: {sim_mt.aggregate_cols(hl.agg.stats(sim_mt.y))}\n\n')
    print(f'\n\nbinarized phenotype stats: {sim_mt.aggregate_cols(hl.agg.stats(sim_mt.y_binarized))}\n\n')
    
    subsets_path = f'{smiles_wd}/subsets.{sim_name}.tsv.gz'

    if hl.hadoop_is_file(subsets_path):
        subsets = hl.import_table(subsets_path, key=['s'], types={'s':hl.tstr}, impute=True, force=True)
        
    else:
        if sim_name[:3]=='bn_':
            # NOTE: This approach does not guarantee equal ascertainment across subsets, it is subject to some stochasticity
            ht = sim_mt.cols()
            cas = get_shuffled_ht(ht=ht, phen='y_binarized', is_cas=True, seed=seed)
            con = get_shuffled_ht(ht=ht, phen='y_binarized', is_cas=False, seed=seed)
            cas = cas.filter(cas.tmp_idx<sum(n_cas_list))
            con = con.filter(con.tmp_idx<(sum(n_list)-sum(n_cas_list)))
            subsets = cas.union(con)
            subsets = subsets.annotate(subset_idx = subsets.pop+1)
            
        else:
            subsets = make_subsets(ht= sim_mt.cols(), 
                                   phen='y_binarized', 
                                   n_list=n_list, 
                                   n_cas_list=n_cas_list)
        
        subsets = subsets.key_by('s')
        subsets.select('subset_idx').export(subsets_path)
        
    sim_mt = sim_mt.annotate_cols(subset_idx = subsets[sim_mt.s].subset_idx)
    
    gwas_ht_list = []
    
    for subset_idx, _ in enumerate(n_cas_list, 1): # subset_idx is 1-indexed

        gwas_ht = run_gwas(mt=sim_mt, 
                           phen='y_binarized', 
                           sim_name=sim_name, 
                           subset_idx=subset_idx, 
                           param_suffix=param_suffix,
                           wd=smiles_wd,
                           is_logreg=True)
        
        gwas_ht_list += [gwas_ht]

    if len(gwas_ht_list)>1:
        
        metagwas_ht = metaanalyze_gwas(subsets=subsets, 
                                       gwas_ht_list=gwas_ht_list, 
                                       sim_name=sim_name,
                                       param_suffix=param_suffix,
                                       wd=smiles_wd)
        
# TODO: Balding-Nichols with two pops, to test allele frequencies in sub-cohorts
   


## get eaf
#mt = mt.annotate_rows(EAF = hl.agg.mean(mt.dosage)/2)
#mt.rows().select('EAF').export(f'{smiles_wd}/ukb.hm3.eaf.n_{n_downsample}.tsv.gz')

# get missingness
#mt = mt.annotate_rows(missing_rate = hl.agg.sum(hl.is_missing(mt.dosage)))
#mt.rows().select('missing_rate').export(f'{smiles_wd}/ukb.hm3.missing_rate.n_{n_downsample}.tsv.gz')

#
#mt_sim = ldscsim.simulate_phenotypes(mt=mt, 
#                                     genotype=mt.dosage, 
#                                     h2=h2, 
#                                     pi=0.01)
#
#
#linreg_ht = hl.linear_regression_rows(y=mt_sim.y, 
#                                    x=mt_sim.dosage,
#                                    covariates=[1])
#
#linreg_path = f'{smiles_wd}/linreg.n_{n_downsample}.h2_{h2}.pi_{pi}.ht'
#linreg_ht.write(linreg_path)
#linreg_ht = hl.read_table(linreg_path)
#linreg_ht.export(f'{smiles_wd}/linreg.n_{n_downsample}.h2_{h2}.pi_{pi}.tsv.gz')

#mt_sim = ldscsim.binarize(mt=mt_sim, 
#                          y=mt_sim.y, 
#                          K=K, 
#                          exact=True)
#
#
#stats = mt_sim.aggregate_cols(hl.agg.stats(mt_sim.y_binarized))
#print(stats)
#
