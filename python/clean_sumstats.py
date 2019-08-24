#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 08:20:11 2019

For cleaning up summary statistics before plotting

@author: nbaya
"""

import pandas as pd
import numpy as np

wd_data = "/Users/nbaya/Documents/lab/smiles/data/"

# Cleaning up T2D (unadjusted for BMI) data and removing unnecessary
# columns to save space
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
scz = pd.read_csv(wd_data + 'pgc.scz.full.2012-04.txt.gz',
                  compression='gzip', delim_whitespace=True)
scz = scz.rename(columns={'hg18chr': 'chr', 'or': 'beta'})
# based on abstract (https://www.ncbi.nlm.nih.gov/pubmed/25056061)
scz['n_complete_samples'] = 150064
scz = scz[scz.CEUaf != '.']  # filter out SNPs with missing CEUaf
# convert OR to correct scale for effect size
scz['beta'] = np.log10(scz['beta'])ƒ
scz.loc[scz.beta > 0, 'EAF'] = scz.loc[scz.beta > 0, 'CEUaf'].astype(float)
scz.loc[scz.beta < 0, 'EAF'] = 1 - scz.loc[scz.beta < 0, 'CEUaf'].astype(float)
scz = scz[['chr', 'EAF', 'beta', 'pval', 'n_complete_samples']]
scz.to_csv(wd_data + 'pgc.scz.full.2012-04.tsv.gz',
           compression='gzip', sep='\t', index=False)

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
df = pd.read_csv(wd_data + 'EUR.CD.gwas_info03_filtered.assoc.gz',
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
df = pd.read_csv(wd_data + 'daner_PGC_SCZ43_mds9.gz.hq2.gz',
                 compression='gzip', delim_whitespace=True)
n_cas = int([x for x in df.columns.values if 'FRQ_A' in x][0].split('_')[2])
n_con = int([x for x in df.columns.values if 'FRQ_U' in x][0].split('_')[2])
#df = df.rename(columns={f'FRQ_U_{n_con}':'EAF','CHR':'chr','P':'pval','BP':'pos'})
# df['beta'] = -np.log10(df.OR) # REVERSE SIGN OF ODDS RATIO because daner
# files use odds ratios that are in reference to A1, meaning: A1*OR = A2
df['beta'] = np.log10(df.OR)
df['n'] = n_cas + n_con
df = df[['chr', 'EAF', 'beta', 'pval', 'n', 'pos']]
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
