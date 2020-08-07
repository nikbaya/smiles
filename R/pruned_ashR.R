#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

library(optparse)
library(ashr)
# library(ggplot2)

smiles_dir = '/Users/nbaya/Documents/lab/smiles'

if (!dir.exists(smiles_dir)) {
  smiles_dir = '/stanley/genetics/users/nbaya/smiles'
}

data_dir = paste0(smiles_dir, '/data')
plots_dir = paste0(smiles_dir, '/plots')

phens = c(
  'standing_height',
  'bmi',
  't2d_bmiadj',
  'ibd',
  'cd',
  'uc',
  'scz',
  'ad',
  'breast_cancer',
  'cad',
  'ldl',
  'hdl',
  'wbc_count',
  'rbc_count',
  'urate',
  'systolic_bp',
  'diastolic_bp',
  'triglycerides'
  )


get_arguments <- function() {
  option_list = list(
    make_option(c("--phen"), type="character", default=NULL, 
                help="phenotype to run", metavar="phenotype"),
    make_option(c("--maf"), type="double", default=NULL, 
                help="MAF to filter by (keep SNPs with MAF>threshold", metavar="MAF"),
    make_option(c("--betahat"), type='character', default='beta',
                help='ash betahat argument', metavar='betahat'),
    make_option(c("--mixcompdist"), type='character', default='normal',
                help='ash mixcompdist argument', metavar='betahat'),
    make_option(c("--pointmass"), type='int', default=TRUE,
                help='whether to include point mass at zero as a component in mixture distribution (1=TRUE, 0=FALSE, default: TRUE)', metavar='')
  )
  opt_parser = OptionParser(option_list=option_list)
  args = parse_args(opt_parser)
  return(args)
}

read_ss <- function(fname) {
  df = read.csv(file = fname, sep = '\t')
  stopifnot('se' %in% colnames(df))
  stopifnot(('A1' %in% colnames(df))&('A2' %in% colnames(df)))
  if (betahat=='var_exp') {
    stopifnot(('var_exp' %in% colnames(df))&('var_exp_se' %in% colnames(df)))
  }
  return(df)
}

plot_beta_comparison <- function(single.chrom.ash) {
  plt <- ggplot(data = single.chrom.ash, aes(abs(betahat), abs(PosteriorMean))) + 
    ggtitle(paste0('chr',chrom)) + geom_point(aes(col=-log10(sebetahat))) + 
    geom_abline(slope=1, intercept=0)
  
  plt_fname = sprintf('beta_comparison.%s.%s%s%s%s%s%s.pdf', 
                       phen,
                       mixcompdist, 
                       paste0('.chr',chrom),
                       if (block_mhc) '.block_mhc' else '',
                       if (betahat!='beta') paste0('.',betahat) else '',
                       if (!pointmass) '.no_pointmass' else '',
                       if (is.null(maf)) '' else sprintf('.maf_%s', maf))
  ggsave(filename=plt_fname, path=paste0(plots_dir,'/'))
}

main <- function(phen, betahat, mixcompdist, pointmass, maf) {
  cat(paste0('Running ash for ',phen,'\n'))

  pruned_fname = sprintf('%s/pruned_gwas.%s.ld_wind_kb_%s%s%s.tsv.gz', data_dir, phen, 
                         ld_wind_kb, 
                         if (block_mhc) '.block_mhc' else '',
                         if (is.null(maf)) '' else sprintf('.maf_%s',maf))
  ld.pruned = read_ss(fname=pruned_fname)
  
  ash_betahat=if (grepl('abs_',betahat)) abs(ld.pruned[,gsub("abs_","",betahat)])  else if (grepl('log_',betahat)) log10(ld.pruned[,gsub("log_","",betahat)]) else ld.pruned[,betahat]
  ash_sebetahat=if (grepl('beta',betahat)) ld.pruned$se else if (grepl('log_',betahat)) ld.pruned[,paste0(gsub("log_","",betahat),'_se')]/(log(10)*ld.pruned[,gsub("log_","",betahat)]) else ld.pruned[,paste0(betahat,'_se')]
  
  ld.pruned.ash <- ash.workhorse(ash_betahat, ash_sebetahat, 
                       mixcompdist = mixcompdist,
                       pointmass=pointmass)
  ld.pruned.g <- get_fitted_g(ld.pruned.ash)
  
  cat(paste0('...Reading in genome-wide file...\n'))
  genomewide_fname = sprintf('%s/gwas.%s%s%s.tsv.gz', 
                             data_dir, phen, 
                             if (block_mhc) '.block_mhc' else '',
                             if (is.null(maf)) '' else sprintf('.maf_%s',maf))
  genome.wide = read_ss(fname=genomewide_fname)
  
  cat(paste('...Starting ash on genome-wide file...\n'))
  cat(sprintf('block_mhc=%s, mixcompdist=%s, betahat=%s, pointmass=%s\n', block_mhc, mixcompdist, betahat, pointmass))

  genome.wide.ash = NULL
  
  for (chrom in seq(1,22)) {
    df = genome.wide[which(genome.wide$chr==chrom),]
    
    ash_betahat=if (grepl('abs_',betahat)) abs(df[,gsub("abs_","",betahat)])  else if (grepl('log_',betahat)) log10(df[,gsub("log_","",betahat)]) else df[,betahat]
    ash_sebetahat=if (grepl('beta',betahat)) df$se else if (grepl('log_',betahat)) df[,paste0(gsub("log_","",betahat),'_se')]/(log(10)*df[,gsub("log_","",betahat)]) else df[,paste0(betahat,'_se')]
    
    df.ash <- ash(ash_betahat, ash_sebetahat, 
                  g=ld.pruned.g, 
                  fixg=TRUE, 
                  mixcompdist = mixcompdist)
    single.chrom.ash = cbind(df[,c('chr','pos','A1','A2')], df.ash$result[,c('betahat','sebetahat','PosteriorMean','PosteriorSD')])
    
    if (is.null(genome.wide.ash)) {
      # print(df.ash$result)
      genome.wide.ash = single.chrom.ash
    } else {
      genome.wide.ash = rbind(genome.wide.ash, single.chrom.ash) #append rows of latest chrom to the end of existing results
    }
    cat(paste('...chrom',chrom,'complete...\n'))
  }
  cat(paste('All chromosomes complete for',phen,'\n'))
  # ash <- cbind(genome.wide[,c('chr','pos', 'A1', 'A2')], genome.wide.ash[,c('betahat','PosteriorMean','PosteriorSD')])
  fname_out = sprintf('%s/tmp-ash.%s.%s%s%s%s%s.tsv.gz', 
                      data_dir, 
                      phen,
                      mixcompdist, 
                      if (block_mhc) '.block_mhc' else '',
                      if (betahat!='beta') paste0('.',betahat) else '',
                      if (!pointmass) '.no_pointmass' else '',
                      if (is.null(maf)) '' else sprintf('.maf_%s', maf))

  write.table(x = genome.wide.ash, file = gzfile(fname_out), sep='\t', quote=F, row.names = FALSE)
  
}


args = get_arguments()

phen=args$phen
  
ld_wind_kb = 500
block_mhc = TRUE
maf=args$maf # set maf=NULL to use data that is not MAF filtered

betahat = args$betahat # Which "betahat" is used for ash. options: beta, abs_beta, var_exp, log_var_exp (default: beta)
mixcompdist = args$mixcompdist #"+uniform" # options: normal, +uniform (if using var_exp instead of beta), halfnormal, halfuniform
pointmass = args$pointmass # if True, include point mass at zero in prior. (set to FALSE if using betahat=log_var_exp)

stopifnot(!(((betahat=='abs_beta')|(betahat=='var_exp'))&(mixcompdist!='+uniform'))) # assert that mixcompdist must be +uniform if fitting on abs(beta) or variance explained

main(phen=phen,
     betahat=betahat, 
     mixcompdist=mixcompdist, 
     pointmass=pointmass,
     maf=maf)
     
     

# pdf('/Users/nbaya/Downloads/test.pdf')
# for (i in 1:10) {
#   plot(ggplot(mapping = aes(x=c(1,2), y=c(1, 1)))+
#   geom_line())
# }
# dev.off()