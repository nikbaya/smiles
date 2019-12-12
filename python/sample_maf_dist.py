#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 09:09:51 2019

Sample 100 SNPs from same MAF distribution as top 100 loci in a phenotype's GWAS

@author: nbaya
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.neighbors import KernelDensity as kd

smiles_wd = '/Users/nbaya/Documents/lab/smiles/data/'


phen = '50_irnt'
n_top_loci = 100
ld_window = int(1000e3)
block_mhc = True
random_betas = False
suffix = f'.top{n_top_loci}loci.ldwindow{int(ld_window/1e3)}kb{".block_mhc" if block_mhc else ""}'

maf = pd.read_csv(smiles_wd+phen+'.maf'+suffix+'.tsv.gz',compression = 'gzip',sep='\t')

toploci_maf = maf[maf.sim_truebeta==1].minor_AF.values

kde = kd(bandwidth=0.06, kernel='gaussian').fit(toploci_maf[:,np.newaxis])
x_plot = np.linspace(0,0.5,1000)[:, np.newaxis]
log_dens = kde.score_samples(x_plot)
pdf = np.exp(log_dens)
plt.plot(x_plot, pdf)
plt.xlim([0,0.5])


sns.kdeplot(maf[maf.sim_truebeta==1].minor_AF,clip=[0,0.5])
plt.xlim([0,0.5])
plt.title('MAF of top loci')

sns.kdeplot(maf.minor_AF,clip=[0,0.5])
plt.xlim([0,0.5])
plt.title(f'MAF of all {len(maf)} SNPs')



cdf = (np.cumsum(pdf)/len(pdf)/2)

plt.plot(x_plot.tolist()+[x_plot[-1]], cdf.tolist()+[1])





