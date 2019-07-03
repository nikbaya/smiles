#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 10:59:21 2019

get binary annotations for whether a variant is in a gene (i.e. coding variant)
Uses Brandon Franco's code to read in gencode data
gencode data: 
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/
https://www.gencodegenes.org/human/release_30.html

@author: nbaya
"""

import hail as hl

wd = '/Users/nbaya/Documents/lab/smiles/data/'

t = hl.import_table('/Users/nbaya/Downloads/gencode.v31.annotation.gff3',no_header=True,impute=True, comment=('#'))

t2 = t.annotate(GFF_Columns = t.f8.split(";").map(lambda x: x.split("=")))
t2 = t2.filter(t2.f2 == "gene")
t2 = t2.annotate(interval=hl.interval(hl.locus(t2.f0, t2.f3, 'GRCh38'), hl.locus(t2.f0, t2.f3, 'GRCh38')))
t2 = t2.annotate(GFF_Columns = hl.dict(t2.GFF_Columns.map(lambda x: (x[0], x[1]))))
t2 = t2.annotate(ID=t2.GFF_Columns["ID"], gene_id=t2.GFF_Columns["gene_id"], gene_name=t2.GFF_Columns["gene_name"], gene_type=t2.GFF_Columns["gene_type"], level=t2.GFF_Columns["level"])
t2 = t2.annotate(type=t2.f2, gene_score=t2.f5, gene_strand=t2.f6, gene_phase=t2.f7)
t2 = t2.drop(t2.GFF_Columns, t2.f8, t2.f0, t2.f1, t2.f2, t2.f3, t2.f4, t2.f5, t2.f6, t2.f7)
                                                                                                              

ss = hl.import_table(wd+'50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz',impute=True,force_bgz=True)
ss = ss.annotate(variant = hl.parse_variant(ss.variant,reference_genome='GRCh38'))

