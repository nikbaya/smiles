#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 10:59:21 2019

get binary annotations for whether a variant is in a gene (i.e. coding variant)
Uses Brandon Franco's (or Liam Abbott's?) code to read in gencode data
gencode data: 
    
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/

https://www.gencodegenes.org/human/release_30.html

@author: nbaya
"""

import hail as hl

wd_data = '/Users/nbaya/Documents/lab/smiles/data'

t = hl.import_table(f'{wd_data}/gencode.v31lift37.annotation.gff3.gz',no_header=True,impute=True, comment=('#'),force=True)
#t = hl.import_table('/Users/nbaya/Downloads/gencode.v31lift37.annotation.gtf',no_header=True,impute=True, comment=('#'))

                                                                                                                    
t2 = t.annotate(GFF_Columns = t.f8.split(";").map(lambda x: x.split("=")))
t2 = t2.filter(t2.f2 == "CDS") # only want coding sequences, not entire genes
t2 = t2.filter(hl.is_valid_locus(t2.f0[3:], t2.f3, 'GRCh37'))
t2 = t2.filter(hl.is_valid_locus(t2.f0[3:], t2.f4, 'GRCh37'))
t2 = t2.annotate(interval=hl.interval(hl.locus(t2.f0[3:], t2.f3, 'GRCh37'), hl.locus(t2.f0[3:], t2.f4, 'GRCh37')))
t2 = t2.annotate(GFF_Columns = hl.dict(t2.GFF_Columns.map(lambda x: (x[0], x[1]))))
t2 = t2.annotate(ID=t2.GFF_Columns["ID"], gene_id=t2.GFF_Columns["gene_id"], 
                 gene_name=t2.GFF_Columns["gene_name"], gene_type=t2.GFF_Columns["gene_type"], 
                 level=t2.GFF_Columns["level"])
t2 = t2.annotate(type=t2.f2, gene_score=t2.f5, gene_strand=t2.f6, gene_phase=t2.f7)
t2 = t2.drop(t2.GFF_Columns, t2.f8, t2.f0, t2.f1, t2.f2, t2.f3, t2.f4, t2.f5, t2.f6, t2.f7)
t2 = t2.key_by(t2.interval)

                                                                                                              
fname = '50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz'
fname = '21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz'
fname = '2443.gwas.imputed_v3.both_sexes.tsv.bgz'
fname = 'EUR.IBD.gwas_info03_filtered.assoc.tsv.gz'
fname = 'EUR.CD.gwas_info03_filtered.assoc.tsv.gz'
fname = 'EUR.UC.gwas_info03_filtered.assoc.tsv.gz'
fname = 'daner_PGC_SCZ43_mds9.tsv.gz'
fname = 'Mahajan.NatGenet2018b.T2D.European.tsv.gz'
fname = 'Mahajan.NatGenet2018b.T2Dbmiadj.European.tsv.gz'
fname = 'MICAD.EUR.ExA.Consortium.PublicRelease.310517.tsv.gz'
fname = 'breastcancer.michailidou2017.b37.cleaned.tsv.gz'
fname = 'AD_sumstats_Jansenetal_2019sept.tsv.gz'
#fname = '20544_2.gwas.imputed_v3.both_sexes.tsv.bgz'
fnames = ['EUR.IBD.gwas_info03_filtered.assoc.tsv.gz',
          'EUR.IBD.gwas_info03_filtered.assoc.tsv.gz',
          'EUR.IBD.gwas_info03_filtered.assoc.tsv.gz',
          'Mahajan.NatGenet2018b.T2Dbmiadj.European.tsv.gz',
          'daner_PGC_SCZ43_mds9.tsv.gz',
          'AD_sumstats_Jansenetal_2019sept.tsv.gz']
for fname in fnames:
    ss0 = hl.import_table(f'{wd_data}/{fname}',impute=True,force=True,types={'chr':hl.tstr})
    if 'variant' in list(ss0.row): 
        variant = ss0.variant.split(':')
        ss = ss0.filter(hl.is_valid_locus(variant[0], 
                                          hl.int(variant[1]),
                                          'GRCh37'))
        locus = ss.variant.split(':')
        ss = ss.annotate(locus = hl.parse_locus(locus[0]+':'+locus[1],reference_genome='GRCh37'))
    elif 'chr' in list(ss0.row) and 'pos' in list(ss0.row):
        ss = ss0.annotate(locus = hl.locus(contig=ss0.chr,pos=ss0.pos,reference_genome='GRCh37'))
    ss = ss.annotate(coding=hl.is_defined(t2[ss.locus]))
    ss = ss.drop('locus')
    ss.export(f"{wd_data}/{fname.split('.tsv')[0]}.coding.tsv{fname.split('.tsv')[1]}")






## Code for GRCh38 sumstats
#                                                                                                              
#
#ss = hl.import_table(wd+'50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz',impute=True,force_bgz=True)
#ss = ss.annotate(hl_variant = hl.parse_variant(ss.variant,reference_genome='GRCh38'))
#ss = ss.filter(hl.is_valid_locus(ss.hl_variant.locus))
#
#t2 = t2.key_by(t2.interval)
#ss = ss.annotate(coding=hl.is_defined(t2[ss.hl_variant.locus]))
#
#ss.show()
#
#ss = ss.annotate(hl_variant = hl.parse_variant)
#
#ss.drop('hl_variant').export(wd+'50_irnt.gwas.imputed_v3.both_sexes.coding.tsv.bgz')


## plot Mary Pat's data

mp_path = '/Users/nbaya/Documents/lab/smiles/MPReeveDataVizualiationFinalProject/'

for phen in ['T2D','SCZ','IBD','AD']:
    df = pd.read_csv(mp_path+f'{phen}_GWASandNaturalSelectionData.csv')
    df['rbeta'] = np.log10(df.OR)
    
    df = df.drop(columns='coding')
    df = df.rename(columns={'coding_bool':'coding','P':'pval'})
    
    df['A1'] = df.A1A2.str.split('/',n=2,expand=True).iloc[:,0]
    df['A2'] = df.A1A2.str.split('/',n=2,expand=True).iloc[:,1]
    df[['chr','pos','A1','A2','rsid','risk_allele','RAF','OR','rbeta','FRQ_A_67390',
        'FRQ_U_94015','INFO','pval','coding']].to_csv(mp_path+f'{phen}_GWASandNaturalSelectionData.tsv',
        sep='\t',index=False)
    
    fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))
    ax.plot(df[~(df.coding_bool==1)].RAF, df[~(df.coding_bool==1)].rbeta,'.',ms=4,c='#1f77b4')
    ax.plot(df[(df.coding_bool==1)].RAF, df[(df.coding_bool==1)].rbeta,'o',markerfacecolor='none',ms=10,c='#1f77b4')
    plt.title(f'Disease: {phen}')
    plt.xlabel('Risk allele frequency')
    plt.ylabel('Effect size')
    plt.legend(['non-coding','coding'],loc=9)
    plt.xlim([-0.01,1.01])
    plt.tight_layout()
    plt.savefig(f'/Users/nbaya/Downloads/{phen}.raf_rbeta.png',dpi=300)
    
    fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))
    ax.plot(df[~(df.coding_bool==1)].RAF, 2*(df[~(df.coding_bool==1)].RAF)*(1-df[~(df.coding_bool==1)].RAF)*df[~(df.coding_bool==1)].rbeta**2,'.',ms=4,c='#1f77b4')
    ax.plot(df[(df.coding_bool==1)].RAF, 2*(df[(df.coding_bool==1)].RAF)*(1-df[(df.coding_bool==1)].RAF)*df[(df.coding_bool==1)].rbeta**2,'o',markerfacecolor='none',ms=10,c='#1f77b4')
    plt.title(f'Disease: {phen}')
    plt.title(f'Disease: {phen}')
    plt.xlabel('Risk allele frequency')
    plt.ylabel('Variance explained')
    plt.xlim([-0.01,1.01])
    plt.legend(['non-coding','coding'],loc=9)
    plt.tight_layout()
    plt.savefig(f'/Users/nbaya/Downloads/{phen}.raf_varexp.png',dpi=300)
    
    
    
