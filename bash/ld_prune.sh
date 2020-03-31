#! /usr/bin/env bash

ss_version=$1

plink=/Users/nbaya/Documents/lab/smiles/smiles/bash/plink/plink

data_wd=/Users/nbaya/Documents/lab/smiles/data
#vcf=${data_wd}/vcf_206_acc_EUROPE.recode.vcf.gz
#bfile=${data_wd}/arabidopsis

# convert from vcf to bfile, while adding variant IDs
#if test ! -f ${bfile}1.bed; then
#	$plink --vcf $vcf --set-missing-var-ids @-# --allow-extra-chr --make-bed --out ${bfile}1
#fi

bfile=${data_wd}/NESP_genotypes.maf_gt_0.nonmissing

# LD clump
# see https://www.cog-genomics.org/plink/1.9/postproc#clump

#ss_version=clumping.t50LL_119acc_SP
ss_path=${data_wd}/${ss_version}.tsv.gz

echo -e "Using file $ss_path as input\n"

clump_p1=0.0001 # p-val threshold for index variants (i.e. sentinel variants) (default: 0.0001)
clump_p2=0.01 # p-val threshold for any variant in a cluster (default: 0.01) 
clump_kb=1000 # clump kb radius (default: 1000)
clump_r2=0.1 # r^2 threshold. SNPs within clump radius (clump_kb) and r^2 w/ index variant > clump_r2 will be assigned to the clump. (default: 0.1) 

if test ! -f ${data_wd}/${ss_version}.clumped; then
	$plink \
		--bfile ${bfile} \
		--chr 1-5 \
		--allow-extra-chr \
		--clump $ss_path \
		--clump-p1 ${clump_p1} \
		--clump-p2 ${clump_p2} \
		--clump-kb ${clump_kb} \
		--clump-r2 ${clump_r2} \
		--out ${data_wd}/${ss_version}
fi

rm ${data_wd}/${ss_version}.nosex
awk '{print $1,$3,$4,$5}' ${data_wd}/${ss_version}.clumped | column -t | gzip > ${data_wd}/${ss_version}.clumped.tsv.gz
rm ${data_wd}/${ss_version}.clumped
