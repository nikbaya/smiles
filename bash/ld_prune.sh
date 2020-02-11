#! /bin/bash

ss_version=$1

PLINK=/Users/nbaya/Documents/lab/smiles/smiles/bash/plink/plink

DATA_WD=/Users/nbaya/Documents/lab/smiles/data
VCF=${DATA_WD}/vcf_206_acc_EUROPE.recode.vcf.gz
BFILE=${DATA_WD}/arabidopsis

# convert from VCF to bfile
if test ! -f ${BFILE}1.bed; then
	$PLINK --vcf $VCF --set-missing-var-ids @-# --allow-extra-chr --make-bed --out ${BFILE}1
fi

# add variant IDs
# if test ! -f ${BFILE}2.bed; then
#	plink --bfile ${BFILE}1 --allow-extra-chr --chr 1-5 --set-missing-var-ids @-# --make-bed --out ${BFILE}2
# fi

# LD clump
# see https://www.cog-genomics.org/plink/1.9/postproc#clump

#ss_version=clumping.t50LL_119acc_SP
ss_path=${DATA_WD}/${ss_version}.tsv.gz

echo -e "Using file $ss_path as input\n"

clump_p1=0.0001 # p-val threshold for index variants (i.e. sentinel variants) (default: 0.0001)
clump_p2=0.01 # p-val threshold for any variant in a cluster (default: 0.01) 
clump_kb=1000 # clump kb radius (default: 1000)
clump_r2=0.1 # r^2 threshold. SNPs within clump radius (clump_kb) and r^2 w/ index variant > clump_r2 will be assigned to the clump. (default: 0.1) 

if test ! -f ${DATA_WD}/${ss_version}.clumped; then
	$PLINK \
		--bfile ${BFILE}1 \
		--chr 1-5 \
		--allow-extra-chr \
		--clump $ss_path \
		--clump-p1 ${clump_p1} \
		--clump-p2 ${clump_p2} \
		--clump-kb ${clump_kb} \
		--clump-r2 ${clump_r2} \
		--out ${DATA_WD}/${ss_version}
fi

rm ${DATA_WD}/${ss_version}.nosex
awk '{print $1,$3,$4,$5}' ${DATA_WD}/${ss_version}.clumped | column -t | gzip > ${DATA_WD}/${ss_version}.clumped.tsv.gz
rm ${DATA_WD}/${ss_version}.clumped
