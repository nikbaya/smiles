## run LD clumping

# fname='Smileplot_GWAS_THP_IB_FT'
ldprune <- function(fname) {
  zz=gzfile(paste0('/Users/nbaya/Documents/lab/smiles/data/',fname, '.tsv.gz'))  
  a=read.csv(zz,header=T,sep='\t')
  
  # a <- readRDS(paste0('/Users/nbaya/Downloads/wetransfer-5821c5/',fname, '.RDS'))
  
  # From Wieters manuscript: SNPs with a significant association with the phenotype were
  # pruned to remove SNPs standing in strong linkage disequilibrium with plink version 1.90 ([53]),
  # following Sohail et al. 2018. The plink -clump function was set to select SNPs below a (GWAS) Pvalue threshold of 0.0001 to start clumps around these index SNPs in windows of 1 Mb and
  # reduces all SNPs with P < 0.01 that are in LD with the index SNPs. The SNP with the lowest p-value
  # in a clump was retained for further analysis.
  
  pval_thresh = 1
  a_tmp <- a[which(a$Pval<pval_thresh),]
  # a_tmp$EAF = a_tmp$AC_1/(a_tmp$AC_0+a_tmp$AC_1)
  # a_tmp$RAF = (a_tmp$EAF)*(a_tmp$beta>0) + (1-a_tmp$EAF)*(a_tmp$beta<0)
  # a_tmp$rbeta = abs(a_tmp$beta)
  # plot(a_tmp$RAF, a_tmp$rbeta, main = fname)
  
  a_tmp$SNP = gsub(' ','',a_tmp$SNP)
  a_tmp$P = a_tmp$Pval
  a_write <- a_tmp[,c('SNP','P')]
  
  b <- gzfile(paste0('/Users/nbaya/Documents/lab/smiles/data/clumping.',fname,'.tsv.gz'))
  write.table(x=a_write,file=b,sep="\t",row.names = FALSE, quote = FALSE)
  
  
  # run bash script
  cmd1 = paste0('/Users/nbaya/Documents/lab/smiles/smiles/bash/ld_prune.sh ','clumping.',fname)
  out <- system(cmd1, intern = TRUE)
  print(fname)
  print(out[28])

}


fnames = c("GxE_FS_117ind_SP", "GxE_t50_117ind_SP", "GxE_SL_117ind_SP",
           "GxE_t50_83ind_NE", "GxE_SL_83ind_NE", "GxE_FS_83ind_NE",
           "FSHL_85acc_NE", "FSHL_119acc_SP", "FSLL_85acc_NE",
           "FSLL_119acc_SP", "SLHL_85acc_NE", "SLHL_119acc_SP",
           "SLLL_85acc_NE", "SLLL_119acc_SP", "t50HL_85acc_NE", 
           "t50HL_119acc_SP", "t50LL_85acc_NE", "t50LL_119acc_SP")

for (fname in fnames) {
  print(paste0("Starting: ",fname))
  ldprune(fname = fname)
}

# -------------------------------
# miscellaneous

fname='t50LL_119acc_SP'
a <- readRDS(paste0('/Users/nbaya/Downloads/',fname,'.RDS'))
b <- gzfile(paste0('/Users/nbaya/Documents/lab/smiles/data/',fname,'.tsv.gz'))
write.table(x=a,file=b,sep="\t")


fname='GxE_FS_83ind_NE'
a = readRDS(paste0('/Users/nbaya/Downloads/GWAS/',fname,'.RDS'))
b <- gzfile(paste0('/Users/nbaya/Documents/lab/smiles/data/',fname,'.tsv.gz'))
write.table(x=a,file=b,sep="\t")

load('/Users/nbaya/Downloads/snpmat_full_beagle_corrected_lowNA_highMAF.R')


# snpmat_beagle... available on Google Drive if not available locally
# Dimensions: 227 x 1448190
# Note: There are no NAs
load('/Users/nbaya/Downloads/snpmat_full_beagle_corrected_lowNA_highMAF.R') 
snpmat <- snpmat_beagle_corrected_filtered_lowNA_highMAF

test = snpmat[,(ncol(snpmat)-1000):ncol(snpmat)]
dim(test)
r_mat = cor(test, method='pearson')
r2_mat = r_mat^2

heatmap(r2_mat,legend='col')

colnames(test)






# write.csv()

