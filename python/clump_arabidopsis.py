#!/usr/bin/env python3

import pandas as pd
import subprocess
from os import path
import scipy.stats as stats

smiles_dir = '/stanley/genetics/users/nbaya/smiles'
data_dir =f'{smiles_dir}/data'
results_dir = f'{smiles_dir}/results'

def clean_sumstats(pheno, overwrite=False):
    
    ss_fname=f'{data_dir}/clumping.{pheno}.tsv.gz'

    if not path.exists(ss_fname) or overwrite:
        df = pd.read_csv(f'{data_dir}/{pheno}.tsv.gz',compression='gzip',sep='\t')
    
        df = df.rename(columns={'Pval':'pval',
                             'Chr':'chr',
                             'Pos':'pos'})
            
        df['EAF'] = df.AC_1/(df.AC_0+df.AC_1)
        df = df[(df.EAF>0)&(df.EAF<1)]
        df['raf'] = (df.EAF)*(df.beta>0) + (1-df.EAF)*(df.beta<0)
        df['rbeta'] = abs(df.beta)
        
        df['SNP'] = df.SNP.str.replace(' ','')
        
        df = df[['chr','pos','SNP','raf','rbeta','pval']]
        
        df.to_csv(ss_fname, compression='gzip',sep='\t',index=False)
    else:
        df = pd.read_csv(ss_fname, compression='gzip',sep='\t')
    
    print(f'variant ct in {pheno} sumstats: {df.shape[0]}')

    return ss_fname

def clump(pheno, pval_thresh, ss_fname, overwrite=False, suffix='', clump_p1 = 0.0001, clump_p2 = 0.01):
#    bash_script = '/stanley/genetics/users/nbaya/smiles/ld_prune.sh'
#    exit_code = subprocess.call([bash_script,ss_fname,str(pval_thresh)])
    
    bfile = f'{data_dir}/NESP_genotypes.maf_gt_0.nonmissing'
    
    clump_kb = 1000
    clump_r2 = 0.1
    
    print(f'Using {ss_fname}')
    
    out = f'{results_dir}/{pheno}.pval_{pval_thresh}{suffix}'
    
    if not path.exists(f'{out}.clumped') or overwrite:
        exit_code = subprocess.call(['plink',
                                     '--silent',
                                     '--bfile', bfile,
                                     '--chr','1-5',
                                     '--clump',ss_fname,
                                     '--clump-field','pval',
                                     '--clump-p1', str(clump_p1),
                                     '--clump-p2', str(clump_p2),
                                     '--clump-kb', str(clump_kb),
                                     '--clump-r2', str(clump_r2),
                                     '--out', out])
        
        if exit_code != 0:
            print(f'PLINK clumping did not return exit code 0')
        
    try: 
        clumped_gwas_fname = f'{results_dir}/clumped_gwas.{pheno}.pval_{pval_thresh}{suffix}.tsv.gz'
        ss0 = pd.read_csv(ss_fname, compression='gzip', sep='\t')
        
        if not path.exists(clumped_gwas_fname) or overwrite:
            plink_clumped = pd.read_csv(f'{out}.clumped', delim_whitespace=True)
            plink_clumped = plink_clumped.rename(columns={'CHR':'chr',
                                                          'BP':'pos'})
    
            clumped = ss0.merge(plink_clumped[['chr','pos','SNP']], on=['chr','pos','SNP'])
                
    
            clumped.to_csv(clumped_gwas_fname, compression='gzip', sep='\t', index=False)
        else:
            clumped = pd.read_csv(clumped_gwas_fname, compression='gzip', sep='\t')
        print(f'PLINK clumping results for {pheno}, pval={pval_thresh}:')
        print(f'\tInitial variant ct: {ss0.shape[0]}')
        print(f'\tFinal variant ct:   {clumped.shape[0]}')
        
    except FileNotFoundError:
        print(f'PLINK clumping failed for {pheno}, pval={pval_thresh}')

def get_varexp_thresh_sumstats(pheno, pval_thresh, ss_fname):
    print(f'\nThresholding by variance explained for {pheno}, pval={pval_thresh}')
    clumped_gwas_fname = f'{results_dir}/clumped_gwas.{pheno}.pval_0.0001.tsv.gz'
    clumped = pd.read_csv(clumped_gwas_fname, compression='gzip', sep='\t')
    clumped.loc[clumped.pval==0, 'pval'] = 5e-324
    clumped = clumped[(clumped.raf!=0)&(clumped.raf!=1)]
    clumped['var_exp'] = 2*clumped.raf*(1-clumped.raf)*clumped.rbeta**2
    clumped['n_hat'] = stats.chi2.isf(q=clumped.pval,df=1)/clumped.var_exp
    n_hat_mean = clumped.n_hat.mean()
    print(f'estimated n_eff: {n_hat_mean}')
    varexp_thresh = stats.chi2.isf(q=pval_thresh,df=1)/n_hat_mean
    
    ss0 = pd.read_csv(ss_fname, compression='gzip', sep='\t')
    ss0['var_exp'] = 2*ss0.raf*(1-ss0.raf)*ss0.rbeta**2
    ss = ss0[ss0.var_exp>varexp_thresh]
    print(f'Var. exp. thresh. results for {pheno}, pval={pval_thresh}:')
    print(f'\tInitial variant ct: {ss0.shape[0]}')
    print(f'\tFinal variant ct:   {ss.shape[0]}')
    
    if ss.shape[0]>0:
        varexp_thresh_ss_fname = f'{data_dir}/{pheno}.pval_{pval_thresh}.varexp_thresh.tsv.gz'
        ss.to_csv(varexp_thresh_ss_fname, compression='gzip', sep='\t', index=False)
    else:
        print('Error: No SNPs remaining')
        varexp_thresh_ss_fname = None
        
    return varexp_thresh_ss_fname
        

if __name__=='__main__':
    phenos = ["GxE_FS_117ind_SP", "GxE_t50_117ind_SP", "GxE_SL_117ind_SP",
               "GxE_t50_83ind_NE", "GxE_SL_83ind_NE", "GxE_FS_83ind_NE",
               "FSHL_85acc_NE", "FSHL_119acc_SP", "FSLL_85acc_NE",
               "FSLL_119acc_SP", "SLHL_85acc_NE", "SLHL_119acc_SP",
               "SLLL_85acc_NE", "SLLL_119acc_SP", "t50HL_85acc_NE",
               "t50HL_119acc_SP", "t50LL_85acc_NE", "t50LL_119acc_SP"]

#    phenos = ["GxE_SL_117ind_SP"]

    pval_thresh_list = [1e-3, 1e-4, 1e-5, 1e-6, 2e-7, 1e-7] # decided to add 1e-
    
    for pheno in phenos:
        print(f'\nStarting cleaning sumstats for {pheno}')
        fname = clean_sumstats(pheno=pheno,
                               overwrite=True)
        print(f'Finished cleaning sumstats for {pheno}')        
        
        df = pd.read_csv(fname, compression='gzip', sep='\t')
        assert 'pval' in df.columns.values, f"{fname} does not contain 'pval' column"
        
        for pval_thresh in pval_thresh_list:    
            print(f'Starting clumping for {pheno} with pval_thresh={pval_thresh}')
            clump(pheno=pheno, 
                  pval_thresh=pval_thresh,
                  clump_p1 = pval_thresh,
                  ss_fname=fname,
                  overwrite=True)

    for pheno in phenos:
        ss_fname = f'{data_dir}/clumping.{pheno}.tsv.gz'
        
        for pval_thresh in pval_thresh_list:
            fname = get_varexp_thresh_sumstats(pheno=pheno,
                                               pval_thresh=pval_thresh,
                                               ss_fname=ss_fname)
            
            if fname is not None:
                clump(pheno=pheno, 
                      pval_thresh=pval_thresh, 
                      ss_fname=fname, 
                      overwrite=True, 
                      suffix='.varexp_thresh',
                      clump_p1 = 1, # set to be 1 because we already thresholded on variance explained 
                      clump_p2 = 1)