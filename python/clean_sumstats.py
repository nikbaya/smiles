#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 08:20:11 2019

For cleaning up summary statistics before plotting

@author: nbaya
"""

import pandas as pd
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt

wd_data = "/Users/nbaya/Documents/lab/smiles/data"

# Cleaning up T2D (unadjusted for BMI) data and removing unnecessary
# columns to save space
# Neff 
t2d = pd.read_csv(wd_data + 'Mahajan.NatGenet2018b.T2D.European.txt.gz',
                  compression='gzip', delimiter='\t')
t2d = t2d.rename(columns={'Chr': 'chr', 'Pos': 'pos', 'Beta': 'beta',
                          'Pvalue': 'pval', 'Neff': 'n_complete_samples'})
t2d['chr'] = t2d['chr'].astype(str)
t2d = t2d[['chr', 'pos', 'EAF', 'beta', 'pval', 'n_complete_samples']]
t2d.to_csv(wd_data + 'Mahajan.NatGenet2018b.T2D.European.tsv.gz',
           compression='gzip', sep='\t', index=False)


# Cleaning up T2D (adjusted for BMI) data and removing unnecessary columns
# to save space
t2d = pd.read_csv(wd_data + 'Mahajan.NatGenet2018b.T2Dbmiadj.European.txt.gz',
                  compression='gzip', delimiter='\t')
t2d = t2d.rename(columns={'Chr': 'chr', 'Pos': 'pos', 'Beta': 'beta',
                          'Pvalue': 'pval', 'Neff': 'n_complete_samples'})
t2d['chr'] = t2d['chr'].astype(str)
t2d = t2d[['chr', 'pos', 'EAF', 'beta', 'pval', 'n_complete_samples']]
t2d.to_csv(wd_data + 'Mahajan.NatGenet2018b.T2Dbmiadj.European.tsv.gz',
           compression='gzip', sep='\t', index=False)

# Cleaning up SCZ1 data and removing unnecessary columns to save space
scz = pd.read_csv(f'{wd_data}/pgc.scz.full.2012-04.txt.gz',
                  compression='gzip', delim_whitespace=True)
scz = scz.rename(columns={'hg18chr': 'chr'})
# based on abstract (https://www.ncbi.nlm.nih.gov/pubmed/25056061)
scz['n_complete_samples'] = 150064
scz = scz[scz.CEUaf != '.']  # filter out SNPs with missing CEUaf
# convert OR to correct scale for effect size
scz['beta'] = np.log10(scz['or'])
scz.loc[scz.beta > 0, 'EAF'] = scz.loc[scz.beta > 0, 'CEUaf'].astype(float)
scz.loc[scz.beta < 0, 'EAF'] = 1 - scz.loc[scz.beta < 0, 'CEUaf'].astype(float)
scz = scz[['chr', 'EAF', 'beta', 'pval', 'n_complete_samples']]
scz.to_csv(wd_data + 'pgc.scz.full.2012-04.tsv.gz',
           compression='gzip', sep='\t', index=False)

## Liu, van Sommeren data (EUR.*.gwas_info_filtered.assoc)
#CHR - chromosome
#SNP - rsid
#BP - position in b37
#A1 - effect allele 
#A2 - other allele
#FRQ_A_XXXX - effect allele frequency in XXXX cases
#FRQ_U_XXXXX - effect allele frequency in XXXXX controls
#INFO - IMPUTE2 info score
#OR - odds ratio
#SE - standard error
#P - P value
#Direction - dir


# Clean up IBD data
ibd = pd.read_csv(wd_data + 'EUR.IBD.gwas_info03_filtered.assoc.gz',
                  compression='gzip', delim_whitespace=True)
n_cas = int([x for x in ibd.columns.values if 'FRQ_A' in x][0].split('_')[2])
n_con = int([x for x in ibd.columns.values if 'FRQ_U' in x][0].split('_')[2])
ibd = ibd.rename(columns={f'FRQ_U_{n_con}': 'EAF',
                          'CHR': 'chr', 'P': 'pval', 'BP': 'pos'})
ibd['beta'] = np.log10(ibd.OR)
ibd['n'] = n_cas + n_con
ibd = ibd[['chr', 'EAF', 'beta', 'pval', 'n', 'pos']]
ibd.to_csv(wd_data + 'EUR.IBD.gwas_info03_filtered.assoc.tsv.gz',
           compression='gzip', sep='\t', index=False)

# Clean up CD (Crohn's disease) data
df = pd.read_csv(f'{wd_data}/EUR.CD.gwas_info03_filtered.assoc.gz',
                 compression='gzip', delim_whitespace=True)
n_cas = int([x for x in df.columns.values if 'FRQ_A' in x][0].split('_')[2])
n_con = int([x for x in df.columns.values if 'FRQ_U' in x][0].split('_')[2])
df = df.rename(columns={f'FRQ_U_{n_con}': 'EAF',
                        'CHR': 'chr', 'P': 'pval', 'BP': 'pos'})
df['beta'] = np.log10(df.OR)
df['n'] = n_cas + n_con
df = df[['chr', 'EAF', 'beta', 'pval', 'n', 'pos']]
df.to_csv(wd_data + 'EUR.CD.gwas_info03_filtered.assoc.tsv.gz',
          compression='gzip', sep='\t', index=False)

# Clean up UC (ulcerative colitis) data
df = pd.read_csv(wd_data + 'EUR.UC.gwas_info03_filtered.assoc.gz',
                 compression='gzip', delim_whitespace=True)
n_cas = int([x for x in df.columns.values if 'FRQ_A' in x][0].split('_')[2])
n_con = int([x for x in df.columns.values if 'FRQ_U' in x][0].split('_')[2])
df = df.rename(columns={f'FRQ_U_{n_con}': 'EAF',
                        'CHR': 'chr', 'P': 'pval', 'BP': 'pos'})
df['beta'] = np.log10(df.OR)
df['n'] = n_cas + n_con
df = df[['chr', 'EAF', 'beta', 'pval', 'n', 'pos']]
df.to_csv(wd_data + 'EUR.UC.gwas_info03_filtered.assoc.tsv.gz',
          compression='gzip', sep='\t', index=False)

# Clean up daner PGC SCZ data
df = pd.read_csv(f'{wd_data}/daner_PGC_SCZ43_mds9.gz.hq2.gz',
                 compression='gzip', delim_whitespace=True)
n_cas = int([x for x in df.columns.values if 'FRQ_A' in x][0].split('_')[2])
n_con = int([x for x in df.columns.values if 'FRQ_U' in x][0].split('_')[2])
df = df.rename(columns={f'FRQ_U_{n_con}':'EAF','CHR':'chr','P':'pval','BP':'pos'})
# df['beta'] = -np.log10(df.OR) # REVERSE SIGN OF ODDS RATIO because daner
# files use odds ratios that are in reference to A1, meaning: A1*OR = A2
df['beta'] = np.log10(df.OR)
df['n'] = n_cas + n_con
df = df[['chr', 'A1','A2','EAF','beta', 'pval', 'n', 'pos']]
df.to_csv(wd_data + 'daner_PGC_SCZ43_mds9.tsv.gz',
          compression='gzip', sep='\t', index=False)


# PGC SCZ cols
['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A_31335', 'FRQ_U_38765',
 'INFO', 'OR', 'SE', 'P', 'ngt', 'Direction', 'HetISqt', 'HetChiSq',
 'HetDf', 'HetPVa']

df.loc[df.OR > 1, 'raf_U'] = df.loc[df.beta > 0, 'FRQ_U_38765']
df.loc[df.OR < 1, 'raf_U'] = 1 - df.loc[df.beta < 0, 'FRQ_U_38765']

df.loc[df.OR > 1, 'raf_A'] = df.loc[df.beta > 0, 'FRQ_A_31335']
df.loc[df.OR < 1, 'raf_A'] = 1 - df.loc[df.beta < 0, 'FRQ_A_31335']

df['FRQ_A_div_FRQ_U'] = df['FRQ_A_31335'] / df['FRQ_U_38765']

# (df['FRQ_A_31335']/df['FRQ_U_38765']>1)&df.beta>


plt.plot(ss.alt_af, ss.beta, '.', alpha=0.1)
plt.title('Standing Height (UKB)')
plt.xlabel('alt allele frequency')
plt.ylabel('beta')

pval_threshold = 5e-8

df_sig = df[df.P < pval_threshold]


plt.plot(df_sig.FRQ_U_38765, df_sig.OR, '.', alpha=0.1)
plt.title(f'SCZ from daner file (p-value < {pval_threshold})')
plt.xlabel('A1 allele frequency in controls')
plt.ylabel('OR with respect to A1 allele')
plt.savefig('/Users/nbaya/Downloads/tmp_fig.png', dpi=600)

plt.plot(df_sig.FRQ_A_div_FRQ_U, np.exp(df_sig.beta), '.')
plt.xlabel('AF_cases/AF_controls (A1 frequency)')
plt.ylabel('OR ')

plt.plot(df_sig.raf_A / df_sig.raf_U, df_sig.beta, '.')


df_sig.FRQ_U_38765.mean()
df_sig.FRQ_A_31335.mean()

df_sig.raf_U.mean()
df_sig.raf_A.mean()

df.shape

df.OR.max()


plt.hist(df.OR, 50)
plt.xlabel('odds ratio')
plt.ylabel('density')

plt.hist(df.EAF, 50)
plt.xlabel('EAF')
plt.ylabel('density')
plt.title('EAF distribution')

plt.hist(df[df.OR > 1].EAF, 50)
plt.xlabel('EAF')
plt.ylabel('density')
plt.title('EAF distribution for SNPs with OR>1')

plt.hist(df[df.OR < 1].EAF, 50)
plt.xlabel('EAF')
plt.ylabel('density')
plt.title('EAF distribution for SNPs with OR<1')


pval_threshold = 5e-8
plt.hist(df[(df.OR > 1) & (df.pval < pval_threshold)].EAF, 50)
plt.xlabel('EAF')
plt.ylabel('density')
plt.title(f'EAF distribution for SNPs with OR>1 and pval<{pval_threshold}')

pval_threshold = 5e-8
plt.hist(1 - df[(df.OR < 1) & (df.pval < pval_threshold)].EAF, 50)
plt.xlabel('EAF')
plt.ylabel('density')
plt.title(f'EAF distribution for SNPs with OR<1 and pval<{pval_threshold}')


plt.hist(ss0[ss0.alt_a].alt_af)
plt.title('RAF distribution for UKB SCZ data\n(without low confidence variants)')
# plt.hist(ss0[(ss0.pval<1)&(~ss0.low_confidence_variant)].raf.dropna())
plt.title(
    'Odds ratio distribution for UKB SCZ data\n(without low confidence variants)')
plt.hist(ss0[(~ss0.low_confidence_variant)].OR, 50)
plt.xlabel('odds ratio')


plt.title('Odds ratio distribution for UKB SCZ data\n(with low confidence variants)')
plt.hist(ss0.beta.dropna(), 50)
plt.xlabel('odds ratio')


plt.xlabel('RAF')
plt.ylabel('density')
plt.xlim([0.99, 1.01])

ss0['OR'] = np.exp(ss0.beta)

plt.hist(df[df.OR < 1].EAF, 50)
plt.hist(df[df.pval < 1e-8].EAF, 50)


ss0.columns.values


# PGC SCZ (Mary Pat's data)
df_new = pd.read_csv(
    '/Users/nbaya/Documents/lab/smiles/MPReeveDataVizualiationFinalProject/GWASandNaturalSelectionData_SCZ.txt', sep='\t')

df_new.shape

df_new.columns.values
plt.plot(df_new.FRQ_U_94015, df_new['OR-rep'], '.')

plt.plot(df_new.FRQ_U_94015, df_new['Ororig'], '.')

merged = df.merge(df_new, left_on='SNP', right_on='rsid')

merged.columns.values

plt.plot(merged.FRQ_U_94015, merged.FRQ_U_38765, '.')
plt.xlabel("Mary Pat's allele frequencies in controls")
plt.ylabel("Ben's allele frequencies in controls")

# flip alleles with Ben OR < 1
merged.loc[merged.OR_x > 1, 'raf_ben'] = merged.FRQ_U_38765
merged.loc[merged.OR_x < 1, 'raf_ben'] = 1 - merged.FRQ_U_38765
plt.plot(merged.FRQ_U_94015, merged.raf_ben, '.')

# flip MP's alleles with OR < 1
merged.loc[merged.Ororig > 1, 'raf_MP'] = merged.FRQ_U_94015
merged.loc[merged.Ororig < 1, 'raf_MP'] = 1 - merged.FRQ_U_94015
merged['beta_ben'] = np.log10(merged.OR_x)
merged.loc[merged.beta_ben > 0, 'RAF_ben'] = merged.FRQ_U_38765
merged.loc[merged.beta_ben < 0, 'RAF_ben'] = 1 - merged.FRQ_U_38765
plt.plot(merged.RAF, merged.RAF_ben, '.')
plt.xlabel("Mary Pat's risk allele frequency in controls")
plt.ylabel("Ben's risk allele frequency in controls")

plt.plot(merged.RAF_ben, merged.beta_ben, '.')

plt.plot(merged.RAF_ben, np.abs(merged.beta_ben), '.')
plt.plot(merged.RAF, np.log10(merged.OR_y), '.')
plt.xlabel('risk allele frequency')
plt.ylabel('effect size')
plt.legend(["Ben's data", "Mary Pat's data"])

merged_sig = merged[merged.P_x < 5e-8]

plt.plot(merged_sig.RAF_ben, np.abs(merged_sig.beta_ben), '.')
plt.plot(merged_sig.RAF, np.log10(merged_sig.OR_y), '.')
plt.xlabel('risk allele frequency')
plt.ylabel('effect size')
plt.legend(["Ben's data", "Mary Pat's data"])


df['beta_ben'] = np.log10(df.OR)
df.loc[df.beta_ben > 0, 'RAF_ben'] = df.FRQ_U_38765
df.loc[df.beta_ben < 0, 'RAF_ben'] = 1 - df.FRQ_U_38765

df_sig = df[df.P < 5e-8]
df_sig.P

plt.plot(df[df.P < 5e-8].RAF_ben, np.abs(df[df.P < 5e-8].beta_ben), '.')


## Read in Beecham 2013 GWAS of multiple sclerosis
# No effect allele freqency (or other effect allele frequency)
ms =  pd.read_csv('/Users/nbaya/Documents/lab/smiles/data/multiplesclerosis.beecham2013.b37.tsv.gz',
                  compression='gzip',sep='\t')

## Read in MICAD.EUR.ExA.Consortium.PublicRelease.310517.txt.gz from
## Citation: Data on coronary artery disease / myocardial infarction have been 
##           contributed by the Myocardial Infarction Genetics and CARDIoGRAM 
##           Exome investigators and have been downloaded from www.CARDIOGRAMPLUSC4D.ORG
mi = pd.read_csv(f'{wd_data}/MICAD.EUR.ExA.Consortium.PublicRelease.310517.txt.gz',
                 delim_whitespace=True, compression='gzip')
mi.loc[mi.log_OR>0,'beta'] = np.log(mi.loc[mi.log_OR>0,'log_OR']) # all SNPs with log_OR=0 have a p-value>=0.9958, so it's okay to exlude these
mi = mi.rename(columns={'effect_allele_freq':'EAF','Pvalue':'pval','N_samples':'n'})
mi = mi[['chr', 'pos', 'effect_allele','other_allele', 'n','EAF','beta', 'pval']]
mi_tmp = mi.dropna(axis=0)
print(f'Dropped {mi.shape[0]-mi_tmp.shape[0]} rows for having NA values')
for col in mi.columns.values:
    print(f'\tNA values in column {col}: {sum(mi[col].isnull())}')
mi = mi_tmp
mi['chr'] = mi['chr'].astype(int)
mi['pos'] = mi['pos'].astype(int)
mi.to_csv(f'{wd_data}/MICAD.EUR.ExA.Consortium.PublicRelease.310517.tsv.gz',
          sep='\t', index=False, compression='gzip')


## Read in Michailidou 2017 breast cancer GWAS
# From: https://www.ebi.ac.uk/gwas/studies/GCST004988
# Discovery sample description
#   76,192 European ancestry cases, 
#   63,082 European ancestry controls

bc = pd.read_csv(wd_data+'breastcancer.michailidou2017.b37.tsv.gz',
                 sep='\t', compression='gzip', low_memory=False)
bc = bc.rename(columns={'chromosome':'chr','base_pair_location':'pos',
                        'effect_allele_frequency':'EAF', 'p_value':'pval'})
bc = bc[['chr', 'pos', 'effect_allele','other_allele','EAF','beta', 'pval']]
cols_to_filter = ['chr', 'pos', 'EAF','beta', 'pval'] #columns which are necessary for downstream analysis and thus must be filtered for NAs
bc_tmp = bc.dropna(axis=0, subset=cols_to_filter)
print(f'Dropped {bc.shape[0]-bc_tmp.shape[0]} rows for having NA values')
for col in cols_to_filter: 
    print(f'\tNA values in column {col}: {sum(bc[col].isnull())}')
bc = bc_tmp
bc['pos'] = bc['pos'].astype(int)
bc.to_csv(wd_data+'breastcancer.michailidou2017.b37.cleaned.tsv.gz',
          sep='\t', index=False, compression='gzip')




## Read fine-mapped version of BMI
# get columns from UKBB_94traits_release1.cols  
# NOTE: effect allele =  alternative = allele2
def clean_finemapped_sumstats(phen: str):
    print(f'cleaning finemapped sumstats for phen {phen}')
    finemapped_cols = pd.read_csv(wd_data+'/UKBB_94traits_release1.cols', header=None,
                                  sep='\t')
    fm_ss = pd.read_csv(wd_data+f'/UKBB_94traits_release1.{phen}.tsv.gz',
                        names=finemapped_cols[0].values, sep='\t', compression='gzip')
    fm_ss['pos'] = fm_ss['variant'].str.split(':',expand=True)[1].astype(int)
    fm_ss['chr'] = fm_ss['chromosome'].str.replace('chr','').astype(int)
    fm_ss.loc[fm_ss.minorallele==fm_ss.allele2,'EAF'] = fm_ss.maf # NOTE: effect allele =  alternative = allele2
    fm_ss.loc[fm_ss.minorallele!=fm_ss.allele2,'EAF'] = 1-fm_ss.maf
    fm_ss['pval'] = (fm_ss.beta_posterior/fm_ss.sd_posterior).apply(lambda x: 2*stats.norm.cdf(-abs(x)))
    fm_ss = fm_ss.rename(columns={'beta_posterior':'beta', 'allele2':'effect_allele',
                                  'allele1':'other_allele'})
    cols_to_filter = ['chr', 'pos', 'EAF','beta', 'pval'] #columns which are necessary for downstream analysis and thus must be filtered for NAs
    fm_ss_tmp = fm_ss.dropna(axis=0, subset=cols_to_filter)
    print(f'Dropped {fm_ss.shape[0]-fm_ss_tmp.shape[0]} rows for having NA values')
    if fm_ss.shape[0]-fm_ss_tmp.shape[0] > 0:
        for col in cols_to_filter: 
            print(f'\tNA values in column {col}: {sum(fm_ss[col].isnull())}')
    fm_ss = fm_ss_tmp[['chr', 'pos', 'effect_allele','other_allele','EAF','beta', 'pval']]
    
    fm_ss.to_csv(wd_data+f'/UKBB_94traits_release1.{phen}.cleaned.tsv.gz',sep='\t',
              compression='gzip',index=False)
    
    return fm_ss

    
bmi_fm = clean_finemapped_sumstats(phen='bmi')
height_fm = clean_finemapped_sumstats(phen='height')

finemapped_cols = pd.read_csv(wd_data+'UKBB_94traits_release1.cols', header=None,
                              sep='\t')
bmi_fm = pd.read_csv(wd_data+'UKBB_94traits_release1.bmi.tsv.gz',
                     names=finemapped_cols[0].values, sep='\t', compression='gzip')
##### bmi_fm['pos'] = bmi_fm['end'].astype(str) # WARNING: variant position is not always the end position of the variant...
bmi_fm['pos'] = bmi_fm['variant'].str.split(':',expand=True)[1].astype(int)
bmi_fm['chr'] = bmi_fm['chromosome'].str.replace('chr','').astype(int)
bmi_fm.loc[bmi_fm.minorallele==bmi_fm.allele2,'EAF'] = bmi_fm.maf
bmi_fm.loc[bmi_fm.minorallele!=bmi_fm.allele2,'EAF'] = 1-bmi_fm.maf
bmi_fm['pval'] = (bmi_fm.beta_posterior/bmi_fm.sd_posterior).apply(lambda x: 2*stats.norm.cdf(-abs(x)))

## Read fine-mapped version of standing height
finemapped_cols = pd.read_csv(wd_data+'UKBB_94traits_release1.cols', header=None,
                              sep='\t')
height_fm = pd.read_csv(wd_data+'UKBB_94traits_release1.height.tsv.gz',
                        names=finemapped_cols[0].values, sep='\t', compression='gzip')


## Check variance explained for I9_CORATHER, 20160 ("Ever Smoked")
# I9_CORATHER: Coronary atherosclerosis (prevalence: 0.0397)

def plot_ukb_varexp(phen: str, phen_desc: str, maf=0.05):
    import matplotlib.pyplot as plt
    
    ss0 = pd.read_csv(f'{wd_data}/{phen}.gwas.imputed_v3.both_sexes.tsv.bgz',
                      compression='gzip', sep='\t')
    ss = ss0.copy()
    ss['alt_allele'] = ss.variant.str.split(':',expand=True)[3]
    ss.loc[ss.alt_allele==ss.minor_allele,'EAF'] = ss.loc[ss.alt_allele==ss.minor_allele].minor_AF # alt allele is effect allele for Neale lab UKB GWAS
    ss.loc[ss.alt_allele!=ss.minor_allele,'EAF'] = 1-ss.loc[ss.alt_allele!=ss.minor_allele].minor_AF
    ss.loc[ss.beta>=0,'RAF'] = ss.loc[ss.beta>=0].EAF# alt allele is effect allele for Neale lab UKB GWAS
    ss.loc[ss.beta<0,'RAF'] = 1-ss.loc[ss.beta<0].EAF
    ss['var_exp'] = 2*(ss.EAF)*(1-ss.EAF)*ss.beta**2
    
    ss1 = ss.loc[ss.minor_AF>maf]
    pvals = np.logspace(-5,-10,6)
    fig, ax1 = plt.subplots(figsize=(6*1.2, 4*1.2))
    for i, pval in enumerate(pvals):
        if i < len(pvals)-1:
            ss_tmp = ss1[(ss1.pval<pval)&(ss1.pval>pvals[i+1])]
    
        ax1.plot(ss_tmp.RAF, ss_tmp.var_exp,'.')
    plt.xlim([0,1])
    plt.title(f'{phen}: {phen_desc}  (maf>{maf})')
    plt.xlabel('RAF')
    plt.ylabel(r'$v$')
    plt.savefig(f'/Users/nbaya/Downloads/{phen}.raf_varexp.png',dpi=300)

    fig, ax2 = plt.subplots(figsize=(6*1.2, 4*1.2))
    for i, pval in enumerate(pvals):
        if i < len(pvals)-1:
            ss_tmp = ss1[(ss1.pval<pval)&(ss1.pval>pvals[i+1])]
    
        ax2.plot(ss_tmp.RAF, abs(ss_tmp.beta),'.')
    plt.xlim([0,1])
    #plt.title('Coronary atherosclerosis: I9_CORATHER (maf>0.05)')
    plt.title(f'{phen}: {phen_desc} (maf>{maf})')
    plt.xlabel('RAF')
    plt.ylabel(r'$\beta$')
    plt.savefig(f'/Users/nbaya/Downloads/{phen}.raf_rbeta.png',dpi=300)
    
plot_ukb_varexp(phen='20160', phen_desc='Ever Smoked')
plot_ukb_varexp(phen='I9_CORATHER', phen_desc='Coronary atherosclerosis')


def plot_clumped_varexp(phen: str, phen_desc: str, info=None, maf=0.05):
    if info is not None:
        assert len({'chr','pos','info'}.intersection(info.columns.values))==3
        
    import matplotlib.pyplot as plt
    
    ss = pd.read_csv(f'{wd_data}/clumped_gwas.{phen}.ld_wind_cm_0.6.block_mhc.pval_1e-05.tsv.gz',
                      compression='gzip', sep='\t')
    
    if info is None and phen in ['standing_height','bmi','standing_height_male','standing_height_female']:
        info = pd.read_csv(f'{wd_data}/variants.tsv.bgz', compression='gzip',sep='\t')[['chr','pos','info']]
        
    ss = ss.merge(info[['chr','pos','info']], on=['chr','pos'])    
    ss['varexp'] = 2*ss.raf*(1-ss.raf)*ss.rbeta**2
    
    print(ss.columns.values)
    ss.loc[ss.pval==0, 'pval'] = 5e-324
    ss1 = ss.loc[(ss.raf>maf)&(1-ss.raf>maf)]
    ss1['n_hat'] = stats.chi2.isf(q=ss1.pval,df=1)/ss1.varexp
    n_hat_mean = ss1.n_hat.mean()
    ss1['n_hat_resid'] = ss1.n_hat - n_hat_mean
    
#    for field in ['raf','varexp', 'beta','info']:
#        fig, ax = plt.subplots(figsize=(6*1.2, 4*1.2))
#        ax.plot(ss1[field], ss1.n_hat_resid, '.')
#        plt.ylabel(r'$\hat{n_j} - \overline{n}$',fontsize=15)
#        plt.xlabel(field,fontsize=15)
#        plt.title(f'{phen_desc}  (maf>{maf})')
#        plt.tight_layout()
#        plt.savefig(f'/Users/nbaya/Downloads/{phen}.nhatresid_{field}.png',dpi=300)

    fig, ax = plt.subplots(figsize=(6*1.2, 4*1.2))
    ax.plot(ss1['info'], ss1.raf, '.')
    plt.xlabel('info',fontsize=15)
    plt.ylabel('raf',fontsize=15)
    plt.title(f'{phen_desc}  (maf>{maf})')
    plt.tight_layout()
    plt.savefig(f'/Users/nbaya/Downloads/{phen}.info_raf.png',dpi=300)

    
    pvals = [1e-5, 1e-6, 1e-7, 5e-8, 1e-8, 1e-9]
#    pvals = np.logspace(-5,-10,6)*1e-4
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, ax1 = plt.subplots(figsize=(6*1.2, 4*1.2))
    varexp_factor = n_hat_mean # ss.n.values.mean() #
#    if phen=='scz':
#        phi = 31335/(31335+38765)
#        varexp_factor = (31335+38765)*phi*(1-phi)
    for i, pval in enumerate(pvals[:-1]):
        ss_tmp = ss1[(ss1.pval<pval)&(ss1.pval>pvals[i+1])]
        ax1.plot(ss_tmp.raf, ss_tmp.varexp,'.', c=colors[i])
    for i, pval in enumerate(pvals[:-1]):    
        varexp_thresh = stats.chi2.isf(q=pval,df=1)/varexp_factor
        plt.axhline(y=varexp_thresh, c=colors[i])
#        varexp_thresh1 = stats.chi2.isf(q=pval,df=1)/(n-2*std)
#        varexp_thresh2 = stats.chi2.isf(q=pval,df=1)/(n+2*std)
#        plt.fill_between(x=[0,1], y1=varexp_thresh1, y2=varexp_thresh2, 
#                         alpha=0.1, color=colors[i])
#        plt.axhline(y=varexp_thresh, c=colors[i], ls='--', alpha=0.5)
#        plt.axhline(y=varexp_thresh, c=colors[i], ls='--', alpha=0.5)
    plt.xlim([0,1])
    plt.title(f'{phen_desc}  (maf>{maf}, '+r'$\overline n =$'+f'{n_hat_mean} )')
    plt.xlabel('RAF')
    plt.ylabel(r'$v$')
    plt.yscale('log')
    plt.legend(pvals[:-1])
    plt.savefig(f'/Users/nbaya/Downloads/{phen}.raf_varexp.png',dpi=300)

variant_info = pd.read_csv(f'{wd_data}/variants.tsv.bgz', compression='gzip',sep='\t')

plot_clumped_varexp(phen='standing_height', phen_desc='Standing Height', info=variant_info)

plot_clumped_varexp(phen='standing_height_male', phen_desc='Standing Height (Male)')
plot_clumped_varexp(phen='standing_height_female', phen_desc='Standing Height (Female)')

plot_clumped_varexp(phen='bmi', phen_desc='BMI', info=variant_info)

scz_orig = pd.read_csv(f'{wd_data}/daner_PGC_SCZ43_mds9.gz.hq2.gz',
                 compression='gzip', delim_whitespace=True)
scz_info = scz_orig.rename(columns={'CHR':'chr','BP':'pos','INFO':'info'})

plot_clumped_varexp(phen='scz', phen_desc='Schizophrenia', info=scz_info[['chr','pos','info']])


ibd_orig = pd.read_csv(f'{wd_data}/EUR.IBD.gwas_info03_filtered.assoc.gz',
                 compression='gzip', delim_whitespace=True)
ibd_orig = ibd_orig.rename(columns={'CHR':'chr','BP':'pos','INFO':'info'})
plot_clumped_varexp(phen='ibd', phen_desc='IBD', info = ibd_orig[['chr','pos','info']])

cd_orig = pd.read_csv(f'{wd_data}/EUR.CD.gwas_info03_filtered.assoc.gz',
                 compression='gzip', delim_whitespace=True)
cd_orig = cd_orig.rename(columns={'CHR':'chr','BP':'pos','INFO':'info'})
plot_clumped_varexp(phen='cd', phen_desc='CD', info=cd_orig[['chr','pos','info']])   

uc_orig = pd.read_csv(f'{wd_data}/EUR.UC.gwas_info03_filtered.assoc.gz',
                 compression='gzip', delim_whitespace=True)
uc_orig = uc_orig.rename(columns={'CHR':'chr','BP':'pos','INFO':'info'})
plot_clumped_varexp(phen='uc', phen_desc='UC', info = uc_orig[['chr','pos','info']])


## Atrial Fibrillation from BBJ, http://jenger.riken.jp/en/result 
# Build: GRCh37, Effect allele is "ALT" (equivalent to A2?)
afib = pd.read_csv('/Users/nbaya/Downloads/MegaGWAS_summary_Asian.txt',delim_whitespace=True)
n_cas = int(afib.N_CASE.max())
n_con = int(afib.N_CTRL.max())
afib = afib[(afib.N_CASE==afib.N_CASE.max()) & (afib.N_CTRL==afib.N_CTRL.max())] # filter to SNPs with maximum cases and controls such that there is no variation in n_case, n_ctrl
afib['EAF'] = 1- (afib.CTRL_FREQ1*n_con+afib.CASE_FREQ1*n_cas)/(n_cas+n_con) # assume CTRL_FREQ1 is frequency of A1 (non effect allele?) in controls
#afib['EAF'] = 1- (afib.CTRL_FREQ1) #*n_con+afib.CASE_FREQ1*n_cas)/(n_cas+n_con) # assume CTRL_FREQ1 is frequency of A1 (non effect allele?) in controls
afib = afib.rename(columns={'POS':'pos','CHR':'chr','PVALUE':'pval','A1':'other_allele',
                            'A2':'effect_allele','BETA':'beta'})
afib = afib[['chr', 'pos', 'effect_allele','other_allele','EAF','beta', 'pval']]
cols_to_filter = ['chr', 'pos', 'EAF','beta', 'pval'] #columns which are necessary for downstream analysis and thus must be filtered for NAs
afib_tmp = afib.dropna(axis=0, subset=cols_to_filter)
print(f'Dropped {afib.shape[0]-afib_tmp.shape[0]} rows for having NA values')
for col in cols_to_filter: 
    print(f'\tNA values in column {col}: {sum(afib[col].isnull())}')
afib = afib_tmp
afib['pos'] = afib['pos'].astype(int)
afib.loc[afib.beta>=0,'RAF'] = afib.loc[afib.beta>=0].EAF# alt allele is effect allele for Neale lab UKB GWAS
afib.loc[afib.beta<0,'RAF'] = 1-afib.loc[afib.beta<0].EAF
afib['var_exp'] = 2*(afib.EAF)*(1-afib.EAF)*afib.beta**2

phen='atrial fibrillation'; phen_desc=f'(cases: {n_cas}, controls: {n_con})'; maf=0.01
afib1 = afib.loc[(afib.EAF>maf) | (1-afib.EAF>maf)]
pvals = [10**(-x) for x in range(5,11)]

fig, ax1 = plt.subplots(figsize=(6*1.2, 4*1.2))
for i, pval in enumerate(pvals):
    if i < len(pvals)-1: # if not at final p-value in list
        afib_tmp = afib1[(afib1.pval<pval)&(afib1.pval>pvals[i+1])]
    else:
        continue
#        afib_tmp = afib1[(afib1.pval<pval)]

    ax1.plot(afib_tmp.RAF, afib_tmp.var_exp,'.')
plt.xlim([0,1])
plt.legend([f'pval<{pval}' for pval in pvals])
plt.title(f'{phen}: {phen_desc}  (maf>{maf}), AF weighted by cas/con')
plt.xlabel('RAF')
plt.ylabel(r'$v$')
plt.savefig(f'/Users/nbaya/Downloads/tmp1.png',dpi=300)

fig, ax2 = plt.subplots(figsize=(6*1.2, 4*1.2))
for i, pval in enumerate(pvals):
    if i < len(pvals)-1:
        afib_tmp = afib1[(afib1.pval<pval)&(afib1.pval>pvals[i+1])]
    else:
        continue
#        afib_tmp = afib1[(afib1.pval<pval)]

    ax2.plot(afib_tmp.RAF, abs(afib_tmp.beta),'.')
plt.xlim([0,1])
#plt.title('Coronary atherosclerosis: I9_CORATHER (maf>0.05)')
plt.title(f'{phen}: {phen_desc} (maf>{maf})')
plt.xlabel('RAF')
plt.ylabel(r'$beta$')




def plot_saige_varexp(df, phen: str):
    ## SAIGE uses alt allele as effect allele, `af` field is freq of effect allele (alt allele)
    n_cas = df.num_cases.max()
    n_con = df.num_controls.max()
#    df = df[(df.num_cases==n_cas) & (df.num_controls==n_con)]
    df = df[df.Is_Converged==1]
    df = df.rename(columns={'#CHROM':'chr','POS':'pos'})
    df['EAF'] = df.af # WARNING: Not clear which is the effect allele
    df = df[['chr', 'pos','EAF','beta', 'pval']]
    cols_to_filter = ['chr', 'pos', 'EAF','beta', 'pval'] #columns which are necessary for downstream analysis and thus must be filtered for NAs
    df_tmp = df.dropna(axis=0, subset=cols_to_filter)
    print(f'Dropped {df.shape[0]-df_tmp.shape[0]} rows for having NA values')
    for col in cols_to_filter: 
        print(f'\tNA values in column {col}: {sum(df[col].isnull())}')
    df = df_tmp
    df['pos'] = df['pos'].astype(int)
    df.loc[df.beta>=0,'RAF'] = df.loc[df.beta>=0].EAF# alt allele is effect allele for Neale lab UKB GWAS
    df.loc[df.beta<0,'RAF'] = 1-df.loc[df.beta<0].EAF
    df['var_exp'] = 2*(df.EAF)*(1-df.EAF)*df.beta**2
    
    phen_desc=f'(max cases: {n_cas}, max controls: {n_con})'
    maf=0.05
    df = df[(df.EAF>maf) & (df.EAF<1-maf)]
    
    pvals = [10**(-x) for x in range(3,11)]
    
    fig, ax2 = plt.subplots(figsize=(6*1.2, 4*1.2))
    for i, pval in enumerate(pvals):
        if i < len(pvals)-1:
            df_tmp = df[(df.pval<pval)&(df.pval>pvals[i+1])]
        else:
            continue
    #        df_tmp = df1[(df1.pval<pval)]
    
        ax2.plot(df_tmp.RAF, abs(df_tmp.beta),'.')
    plt.xlim([0,1])
    plt.legend([f'pval<{pval}' for pval in pvals])
    plt.title(f'{phen}: {phen_desc} (maf>{maf})')
    plt.xlabel('RAF')
    plt.ylabel(r'$beta$')
    
    fig, ax1 = plt.subplots(figsize=(6*1.2, 4*1.2))
    for i, pval in enumerate(pvals):
        if i < len(pvals)-1: # if not at final p-value in list
            df_tmp = df[(df.pval<pval)&(df.pval>pvals[i+1])]
        else:
            continue
    #        df_tmp = df1[(df1.pval<pval)]
    
        ax1.plot(df_tmp.RAF, df_tmp.var_exp,'.')
    plt.xlim([0,1])
    plt.legend(pvals, loc='upper right')
    plt.title(f'{phen}: {phen_desc}  (maf>{maf})')
    plt.xlabel('RAF')
    plt.ylabel(r'$v$')
    #plt.savefig(f'/Users/nbaya/Downloads/{phen}.raf_varexp.png',dpi=300)
    
    
## Myocardial infarction using SAIGE (PheCode: 411.2)
mi = pd.read_csv('/Users/nbaya/Downloads/PheCode_411.2_SAIGE_MACge20.txt.vcf.gz',
                 delim_whitespace=True, compression='gzip')

## Coronary atherosclerosis using SAIGE (PheCode: 411.4)
ca = pd.read_csv('/Users/nbaya/Downloads/PheCode_411.4_SAIGE_MACge20.txt.vcf.gz',
                 delim_whitespace=True, compression='gzip')

plot_saige_varexp(df=mi, phen='myocardial infarction')
plt.savefig('/Users/nbaya/Downloads/image1.png',dpi=300)

plot_saige_varexp(df=ca, phen='coronoary atherosclerosis')
plt.savefig('/Users/nbaya/Downloads/image2.png',dpi=300)



## Read simulation results

def plot_sim_varexp(df, eaf, phen, phen_desc):
    if 'EAF' not in df.columns.values:
        df= df.merge(eaf, on=['locus','alleles'])

    maf = 0.05
    df = df[(df.EAF>maf) & (df.EAF<1-maf)]
        
    df = df.rename(columns={'p_value':'pval'})
    
    df.loc[df.beta>=0,'RAF'] = df.loc[df.beta>=0].EAF# alt allele is effect allele for Neale lab UKB GWAS
    df.loc[df.beta<0,'RAF'] = 1-df.loc[df.beta<0].EAF
    
#    df['beta_CI_upper'] = df.beta+2*df.standard_error
#    df['beta_CI_lower'] = df.beta-2*df.standard_error
    
    
    df['var_exp'] = 2*(df.EAF)*(1-df.EAF)*df.beta**2
#    df['var_exp_CI_upper'] = 2*(df.EAF)*(1-df.EAF)*df.beta_CI_upper**2
#    df['var_exp_CI_lower'] = 2*(df.EAF)*(1-df.EAF)*df.beta_CI_lower**2

    pvals = [10**(-x) for x in range(3,11)]
    
    fig, ax1 = plt.subplots(figsize=(6*1.2, 4*1.2))
    for i, pval in enumerate(pvals):
        if i < len(pvals)-1: # if not at final p-value in list
            df_tmp = df[(df.pval<pval)&(df.pval>pvals[i+1])]
        else:
            continue
    #        df_tmp = df[(df.pval<pval)]
    
        ax1.plot(df_tmp.RAF, df_tmp.var_exp,'.')
#        plt.errorbar(x=df_tmp.RAF, y=df_tmp.var_exp,yerr=2*(df_tmp.EAF)*(1-df_tmp.EAF)*(df_tmp.standard_error)**2,
#                     elinewidth=0.5,linestyle='None')
#        ax1.hist2d(df_tmp.RAF, df_tmp.var_exp, bins=8, weights=df_tmp.standard_error)
#        sns.kdeplot(df_tmp.RAF, df_tmp.var_exp)
    plt.legend(pvals)
    plt.title(phen_desc)
    plt.xlabel('RAF')
    plt.ylabel(r'$v$')
    plt.xlim([0,1])

m1 = df_tmp.RAF.values
m2 = df_tmp.var_exp.values

xmin = m1.min()
xmax = m1.max()
ymin = m2.min()
ymax = m2.max()
 
X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([m1, m2])
kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)

fig, ax = plt.subplots()
ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
          extent=[xmin, xmax, ymin, ymax])
ax.plot(m1, m2, 'k.', markersize=2)
#ax.set_xlim([xmin, xmax])
#ax.set_ylim([ymin, ymax])
plt.show()


    plt.xlim([0,1])
    plt.legend([f'pval<{pval}' for pval in pvals])
    plt.title(f'{phen}: {phen_desc}  (maf>{maf})')
    plt.xlabel('RAF')
    plt.ylabel(r'$v$')
#    plt.savefig(f'/Users/nbaya/Downloads/{phen}.{phen_desc}.raf_varexp.png',dpi=300)
    
#    fig, ax1 = plt.subplots(figsize=(6*1.2, 4*1.2))
#    for i, pval in enumerate(pvals):
#        if i < len(pvals)-1: # if not at final p-value in list
#            df_tmp = df[(df.pval<pval)&(df.pval>pvals[i+1])]
#        else:
#            continue
#    #        df_tmp = df[(df.pval<pval)]
#    
#        ax1.plot(df_tmp.RAF, df_tmp.var_exp,'.')
##        plt.errorbar(x=df_tmp.RAF, y=df_tmp.var_exp,yerr=2*(df_tmp.EAF)*(1-df_tmp.EAF)*df_tmp.standard_error**2,
##                     elinewidth=0.1,linestyle='None')
#    plt.xlim([0,1])
#    plt.legend([f'pval<{pval}' for pval in pvals])
#    plt.title(f'{phen}: {phen_desc}  (maf>{maf})')
#    plt.xlabel('RAF')
#    plt.ylabel(r'$v$')
#    plt.savefig(f'/Users/nbaya/Downloads/{phen}.{phen_desc}.raf_beta.png',dpi=300)
    
#    fig, ax2 = plt.subplots(figsize=(6*1.2, 4*1.2))
#    for i, pval in enumerate(pvals):
#        if i < len(pvals)-1:
#            df_tmp = df[(df.pval<pval)&(df.pval>pvals[i+1])]
#        else:
#            continue
#    #        df_tmp = df[(df.pval<pval)]
#    
#        plt.errorbar(x=df_tmp.RAF, y=df_tmp.var_exp,yerr=df_tmp.standard_error,
#                     elinewidth=0.1, linestyle='None')
#    plt.xlim([0,1])
#    plt.legend([f'pval<{pval}' for pval in pvals])
#    plt.title(f'{phen}: {phen_desc} (maf>{maf})')
#    plt.xlabel('RAF')
#    plt.ylabel(r'$\beta$')
    

    
linreg = pd.read_csv('/Users/nbaya/Downloads/linreg.n_100000.h2_0.6.pi_0.01.tsv.gz',
                     compression='gzip', sep='\t')

logreg = pd.read_csv('/Users/nbaya/Downloads/logreg.n_100000.h2_0.6.pi_0.01.tsv.gz',
                     compression='gzip', sep='\t')

eaf = pd.read_csv('/Users/nbaya/Downloads/ukb.hm3.eaf.n_100000.tsv.gz',
                  compression='gzip',sep='\t')

missing = pd.read_csv('/Users/nbaya/Downloads/ukb.hm3.missing_rate.n_100000.tsv.gz',
                  compression='gzip',sep='\t')

logreg1 = pd.read_csv(f'{wd_data}/gwas.logreg.nonmeta_ascertained.subset_1.h2_0.6.pi_0.01.K_0.05.seed_1.tsv.gz',
                     compression='gzip', sep='\t')
logreg2 = pd.read_csv(f'{wd_data}/gwas.logreg.nonmeta_ascertained2.subset_1.h2_0.6.pi_0.01.K_0.05.seed_1.tsv.gz',
                     compression='gzip', sep='\t')



meta_logreg1 = pd.read_csv(f'{wd_data}/gwas.logreg.meta1.h2_0.6.pi_0.01.K_0.05.seed_1.tsv.gz',
                     compression='gzip', sep='\t')
meta_logreg1 = meta_logreg1.rename(columns={'meta_pval':'p_value',
                                'meta_beta':'beta',
                                'meta_EAF':'EAF'})
meta_logreg2 = pd.read_csv(f'{wd_data}/gwas.logreg.meta2.h2_0.6.pi_0.01.K_0.05.seed_1.tsv.gz',
                     compression='gzip', sep='\t')
meta_logreg2 = meta_logreg2.rename(columns={'meta_pval':'p_value',
                                'meta_beta':'beta',
                                'meta_EAF':'EAF'})
meta_logreg3 = pd.read_csv(f'{wd_data}/gwas.logreg.meta3.h2_0.6.pi_0.01.K_0.05.seed_1.tsv.gz',
                     compression='gzip', sep='\t')
meta_logreg3 = meta_logreg3.rename(columns={'meta_pval':'p_value',
                                'meta_beta':'beta',
                                'meta_EAF':'EAF'})
    
bn1 = pd.read_csv(f'{wd_data}/smiles_gwas.logreg.bn_nonmeta1.bn.npops_3.nvars_100000.nsim_320000.h2_0.6.pi_0.01.K_0.05.seed_1.tsv.gz',
                     compression='gzip', sep='\t')
bn1 = bn1.rename(columns={'meta_pval':'p_value',
                                'meta_beta':'beta',
                                'meta_EAF':'EAF'})

plot_sim_varexp(df=linreg, eaf=eaf, phen='sim', phen_desc='quant trait linreg, h2=0.6, pi=0.01, n=100k')
plot_sim_varexp(df=logreg, eaf=eaf, phen='sim', phen_desc='logreg, h2=0.6, pi=0.01, K=0.5, P=0.5, n=100k')
plt.savefig('/Users/nbaya/Downloads/image1.png',dpi=300)

plot_sim_varexp(df=logreg1, eaf=None, phen='sim', phen_desc='non-meta logreg: h2=0.6, pi=0.01, K=0.05, P=0.5, n=36k')
plot_sim_varexp(df=logreg2, eaf=None, phen='sim', phen_desc='non-meta logreg: h2=0.6, pi=0.01, K=0.05, P=0.5, n=30k')
plot_sim_varexp(df=meta_logreg1, eaf=None, phen='sim', phen_desc='Uniform cohort size (6 x 5k), logreg: h2=0.6, pi=0.01, K=0.05, P=0.5')
plot_sim_varexp(df=meta_logreg2, eaf=None, phen='sim', phen_desc='Mixed cohort size, logreg: h2=0.6, pi=0.01, K=0.05, P=0.5')
plot_sim_varexp(df=meta_logreg3, eaf=None, phen='sim', phen_desc='Uniform cohort size (30 x 1k), logreg: h2=0.6, pi=0.01, K=0.05, P=0.5')

plot_sim_varexp(df=bn1, eaf=None, phen='sim', phen_desc='BN sim, mixed cohort size (total: 30k), logreg: h2=0.6, pi=0.01, K=0.05')
plt.yscale('log')



## try using 1KG European AFs with cleaned sumstats
#ref = pd.read_csv('/Users/nbaya/Documents/lab/genotype-qc/eur_1kg.frq',
#                      delim_whitespace=True)
#ref = eur_1kg.rename(columns={'MAF':'FRQ_A1'})

# For Liu, van Sommeren data (CD, IBD, UC), EAF is FRQ_A1 because A1 is the effect allele
orig = pd.read_csv(f'{wd_data}/EUR.CD.gwas_info03_filtered.assoc.gz',
                 compression='gzip', delim_whitespace=True)
#df['FRQ_A1'] = 
orig = orig.rename(columns={[x for x in orig.columns.values if 'FRQ_U' in x][0]:'FRQ_A1'}) # use AFs of controls

clumped = pd.read_csv(f'{wd_data}/clumped_gwas.cd.ld_wind_cm_0.6.block_mhc.pval_1e-05.tsv.gz',
                 compression='gzip',sep='\t')
clumped = clumped.rename(columns={'chr':'CHR','pos':'BP'})

clumped_merge = clumped.merge(orig, on=['CHR','BP'])

gwas_ref = ref.merge(clumped_merge, on=['CHR','BP'], suffixes=['_ref','_gwas'])

gwas_ref.loc[(gwas_ref.A1_gwas==gwas_ref.A1_ref),'EAF_A1_ref_flipped'] = gwas_ref.FRQ_A1_ref
gwas_ref.loc[(gwas_ref.A1_gwas==gwas_ref.A2_ref),'FRQ_A1_ref_flipped'] = 1-gwas_ref.FRQ_A1_ref


gwas_ref['rbeta'] = abs(gwas_ref['beta'])
gwas_ref['beta'] = np.log10(gwas_ref['OR'])

gwas_ref['EAF_ref'] = gwas_ref['FRQ_A1_ref_flipped']
gwas_ref.loc[gwas_ref.beta>=0,'RAF_ref'] = gwas_ref.loc[gwas_ref.beta<0,'EAF_ref']
gwas_ref.loc[gwas_ref.beta<0,'RAF_ref'] = 1-gwas_ref.loc[gwas_ref.beta<0,'EAF_ref']

gwas_ref['var_exp_ref'] = 2*(gwas_ref.RAF_ref)*(1-gwas_ref.RAF_ref)*gwas_ref.beta**2


pvals = [10**(-x) for x in range(5,11)]
df = gwas_ref
    
fig, ax1 = plt.subplots(figsize=(6*1.2, 4*1.2))
for i, pval in enumerate(pvals):
    if i < len(pvals)-1: # if not at final p-value in list
        df_tmp = df[(df.P<pval)&(df.P>pvals[i+1])]
    else:
        continue
#        df_tmp = df[(df.pval<pval)]

    ax1.plot(df_tmp.RAF_ref, df_tmp.var_exp_ref,'.')
#        plt.errorbar(x=df_tmp.RAF, y=df_tmp.var_exp,yerr=2*(df_tmp.EAF)*(1-df_tmp.EAF)*(df_tmp.standard_error)**2,
#                     elinewidth=0.5,linestyle='None')
#        ax1.hist2d(df_tmp.RAF, df_tmp.var_exp, bins=8, weights=df_tmp.standard_error)
#        sns.kdeplot(df_tmp.RAF, df_tmp.var_exp)
plt.legend(pvals)
plt.title('Crohns: AF from 1KG European')
plt.xlabel('RAF')
plt.ylabel(r'$v$')
plt.xlim([0,1])
plt.yscale('log')
plt.savefig('/Users/nbaya/Downloads/ref_af.png',dpi=300)


gwas_ref['EAF_gwas'] = gwas_ref['FRQ_A1_gwas']
gwas_ref.loc[gwas_ref.beta>=0,'EAF_gwas'] = gwas_ref.loc[gwas_ref.beta<0,'EAF_gwas']
gwas_ref.loc[gwas_ref.beta<0,'EAF_gwas'] = 1-gwas_ref.loc[gwas_ref.beta<0,'EAF_gwas']
gwas_ref['var_exp_gwas'] = 2*(gwas_ref.RAF_ref)*(1-gwas_ref.RAF_ref)*gwas_ref.beta**2

df = gwas_ref
fig, ax1 = plt.subplots(figsize=(6*1.2, 4*1.2))
for i, pval in enumerate(pvals):
    if i < len(pvals)-1: # if not at final p-value in list
        df_tmp = df[(df.P<pval)&(df.P>pvals[i+1])]
    else:
        continue
#        df_tmp = df[(df.pval<pval)]

    ax1.plot(df_tmp.RAF, df_tmp.var_exp_gwas,'.')
#        plt.errorbar(x=df_tmp.RAF, y=df_tmp.var_exp,yerr=2*(df_tmp.EAF)*(1-df_tmp.EAF)*(df_tmp.standard_error)**2,
#                     elinewidth=0.5,linestyle='None')
#        ax1.hist2d(df_tmp.RAF, df_tmp.var_exp, bins=8, weights=df_tmp.standard_error)
#        sns.kdeplot(df_tmp.RAF, df_tmp.var_exp)
#plt.legend(pvals)
plt.title('Crohns: AF from GWAS controls')
plt.xlabel('RAF')
plt.ylabel(r'$v$')
plt.xlim([0,1])
plt.yscale('log')
plt.savefig('/Users/nbaya/Downloads/gwas_controls_af.png',dpi=300)

plt.plot(clumped.raf, 2*(1-clumped.raf)*clumped.raf*clumped.rbeta**2,'.')
plt.yscale('log')




# Try using UKB standing height AFs
ref = pd.read_csv(f'{wd_data}/50_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz',
                      compression='gzip',sep='\t')
if 'chr' not in ref.columns.values:
    ref['chr'] = ref.variant.str.split(':',n=1,expand=True).iloc[:,0]
ref = ref[ref.chr.isin([str(x) for x in range(1,23)])] # only need to go up to chr 22 because genetic map only covers 22 chr
ref['chr'] = ref['chr'].astype(int) # previous filtering step is necessary for converting to int for cases of chr='X' or 'MT'

if 'pos' not in ref.columns.values:
    ref['pos'] = ref.variant.str.split(':',n=2,expand=True).iloc[:,1]
ref['pos'] = ref['pos'].astype(int)
#maf_thresh = 0.01
#ref = ref[ref.minor_AF>maf_thresh] # MAF filter
ref = ref.rename(columns={'minor_AF':'FRQ_A1',
                          'minor_allele':'A1'})
    
ref = ref[['chr','pos','A1','FRQ_A1']] # only include necessary columns
    

disease_desc_dict = {'uc':'Ulcerative colitis',
                     'ibd':'IBD',
                     'cd':"Crohn's",
                     'scz':'Schizophrenia'}

for disease in ['uc']: #disease = 'uc'
    orig = pd.read_csv(f'{wd_data}/EUR.{disease.upper()}.gwas_info03_filtered.assoc.gz',
                     compression='gzip', delim_whitespace=True)
    orig = orig.rename(columns={[x for x in orig.columns.values if 'FRQ_U' in x][0]:'FRQ_A1',# use AFs of controls
                                'CHR':'chr',
                                'BP':'pos'}) 
    
    clumped = pd.read_csv(f'{wd_data}/clumped_gwas.{disease}.ld_wind_cm_0.6.block_mhc.pval_1e-05.tsv.gz',
                     compression='gzip',sep='\t')
    gwas = clumped.merge(orig, on=['chr','pos'], suffixes=('_clumped','_orig'))
    
    gwas_ref = ref.merge(gwas, on=['chr','pos'], suffixes=('_ref','_gwas'))
    
    gwas_ref.loc[(gwas_ref.A1_gwas==gwas_ref.A1_ref),'FRQ_A1_ref_flipped'] = gwas_ref.FRQ_A1_ref
    gwas_ref.loc[(gwas_ref.A1_gwas!=gwas_ref.A1_ref),'FRQ_A1_ref_flipped'] = 1-gwas_ref.FRQ_A1_ref
    
    for version in ['ref_flipped','gwas']:
        gwas_ref[f'EAF_{version}'] = gwas_ref[f'FRQ_A1_{version}']
        gwas_ref.loc[gwas_ref.beta>=0,f'RAF_{version}'] = gwas_ref.loc[gwas_ref.beta>=0,f'EAF_{version}']
        gwas_ref.loc[gwas_ref.beta<0,f'RAF_{version}'] = 1-gwas_ref.loc[gwas_ref.beta<0,f'EAF_{version}']
        gwas_ref[f'var_exp_{version}'] = 2*(gwas_ref[f'RAF_{version}'])*(1-gwas_ref[f'RAF_{version}'])*gwas_ref.beta**2
    
    pvals = [10**(-x) for x in range(5,11)]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    maf_threshold = 0.05
    for version in ['ref_flipped','gwas']:
        fig, ax1 = plt.subplots(figsize=(6*1.2, 4*1.2))
        df = gwas_ref[(gwas_ref[f'RAF_{version}']>maf_threshold)&(1-gwas_ref[f'RAF_{version}']>maf_threshold)]
        # remove outlier to fix y-axis scale
        df = df[~(df.var_exp_ref_flipped==df.var_exp_ref_flipped.min())]
        for i, pval in enumerate(pvals):
            if i < len(pvals)-1: # if not at final p-value in list
                df_tmp = df[(df.P<pval)&(df.P>pvals[i+1])]
            else:
                continue
        #        df_tmp = df[(df.pval<pval)]
        
            ax1.plot(df_tmp[f'RAF_{version}'], df_tmp[f'var_exp_{version}'],'.')
    #                 'x' if version is 'ref' else '.',
    #                 c=colors[i])
        #        plt.errorbar(x=df_tmp.RAF, y=df_tmp.var_exp,yerr=2*(df_tmp.EAF)*(1-df_tmp.EAF)*(df_tmp.standard_error)**2,
        #                     elinewidth=0.5,linestyle='None')
        #        ax1.hist2d(df_tmp.RAF, df_tmp.var_exp, bins=8, weights=df_tmp.standard_error)
        #        sns.kdeplot(df_tmp.RAF, df_tmp.var_exp)
        #plt.legend(pvals)
    #    plt.title(f'Crohns: AF from {version}')
    #    plt.title(f'IBD: AF from {version}')
        plt.title(f'{disease_desc_dict[disease]}: AF from {version} (maf>{maf_threshold})')
        plt.xlabel(f'RAF_{version}')
        plt.ylabel(r'$v$')
        plt.xlim([0,1])
        plt.yscale('log')
        plt.savefig(f'/Users/nbaya/Downloads/{disease}.{version}_af.png',dpi=300)

        
    
    
    
    
# PGC SCZ

disease = 'scz'
    
orig = pd.read_csv(f'{wd_data}/daner_PGC_SCZ43_mds9.gz.hq2.gz',
                 compression='gzip', delim_whitespace=True)
#orig = orig.rename(columns={[x for x in orig.columns.values if 'FRQ_U' in x][0]:'FRQ_A1',# use AFs of controls
#                            'CHR':'chr',
#                            'BP':'pos'}) 
# "OR reported with respect to A1 allele" -> A1 allele is the effect allele (https://www.cell.com/ajhg/pdf/S0002-9297(12)00044-4.pdf)
# FRQ_A1 field is allele freq in unaffected samples
gwas = pd.read_csv(f'{wd_data}/clumped_gwas.{disease}.ld_wind_cm_0.6.block_mhc.pval_1e-05.tsv.gz',
                      compression='gzip',sep='\t')
gwas['FRQ_A1'] = gwas['EAF']
                      
#gwas = clumped.merge(orig, on=['chr','pos'], suffixes=('_clumped','_orig'))


ref = pd.read_csv(f'{wd_data}/50_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz',
                      compression='gzip',sep='\t')
if 'chr' not in ref.columns.values:
    ref['chr'] = ref.variant.str.split(':',n=1,expand=True).iloc[:,0]
ref = ref[ref.chr.isin([str(x) for x in range(1,23)])] # only need to go up to chr 22 because genetic map only covers 22 chr
ref['chr'] = ref['chr'].astype(int) # previous filtering step is necessary for converting to int for cases of chr='X' or 'MT'

if 'pos' not in ref.columns.values:
    ref['pos'] = ref.variant.str.split(':',n=2,expand=True).iloc[:,1]
ref['pos'] = ref['pos'].astype(int)
#maf_thresh = 0.01
#ref = ref[ref.minor_AF>maf_thresh] # MAF filter
ref = ref.rename(columns={'minor_AF':'FRQ_A1',
                          'minor_allele':'A1'})
    
ref = ref[['chr','pos','A1','FRQ_A1']] # only include necessary columns

gwas_ref = ref.merge(gwas, on=['chr','pos'], suffixes=('_ref','_gwas'))

gwas_ref.loc[(gwas_ref.A1_gwas==gwas_ref.A1_ref),'FRQ_A1_ref_flipped'] = gwas_ref.FRQ_A1_ref
gwas_ref.loc[(gwas_ref.A1_gwas!=gwas_ref.A1_ref),'FRQ_A1_ref_flipped'] = 1-gwas_ref.FRQ_A1_ref

for version in ['ref_flipped','gwas']:
    gwas_ref[f'EAF_{version}'] = gwas_ref[f'FRQ_A1_{version}']
    gwas_ref.loc[gwas_ref.beta>=0,f'RAF_{version}'] = gwas_ref.loc[gwas_ref.beta>=0,f'EAF_{version}']
    gwas_ref.loc[gwas_ref.beta<0,f'RAF_{version}'] = 1-gwas_ref.loc[gwas_ref.beta<0,f'EAF_{version}']
    gwas_ref[f'var_exp_{version}'] = 2*(gwas_ref[f'RAF_{version}'])*(1-gwas_ref[f'RAF_{version}'])*gwas_ref.beta**2

#gwas_ref1 = gwas_ref[~(gwas_ref.var_exp_ref_flipped==gwas_ref.var_exp_ref_flipped.min())]


pvals = [10**(-x) for x in range(5,11)]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
maf_threshold = 0.05
for version in ['ref_flipped','gwas']:
    fig, ax1 = plt.subplots(figsize=(6*1.2, 4*1.2))
    df = gwas_ref[(gwas_ref1[f'RAF_{version}']>maf_threshold)&(1-gwas_ref[f'RAF_{version}']>maf_threshold)]
    # remove outlier to fix y-axis scale
    df = df[~(df.var_exp_ref_flipped==df.var_exp_ref_flipped.min())]
    for i, pval in enumerate(pvals):
        if i < len(pvals)-1: # if not at final p-value in list
            df_tmp = df[(df.pval<pval)&(df.pval>pvals[i+1])]
        else:
            continue
    #        df_tmp = df[(df.pval<pval)]
    
        ax1.plot(df_tmp[f'RAF_{version}'], df_tmp[f'var_exp_{version}'],'.')
    plt.title(f'{disease_desc_dict[disease]}: AF from {version} (maf>{maf_threshold})')
    plt.xlabel(f'RAF_{version}')
    plt.ylabel(r'$v$')
    plt.xlim([0,1])
    plt.yscale('log')
    plt.savefig(f'/Users/nbaya/Downloads/{disease}.{version}_af.png',dpi=300)
