#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

library(ashr)

data_wd = '/Users/nbaya/Documents/lab/smiles/data'
if (!dir.exists(data_wd)) {
  data_wd = '/stanley/genetics/users/nbaya/smiles/data'
}

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

ld_wind_kb = 500
block_mhc = TRUE
maf=0.01 # set maf=NULL to use data that is not MAF filtered

betahat = 'log_var_exp' # Which "betahat" is used for ash. options: beta, abs_beta, var_exp, log_var_exp (default: beta)
pointmass=TRUE # if True, include point mass at zero in prior
mixcompdist = '+uniform' #"+uniform" # options: normal, +uniform (if using var_exp instead of beta), halfnormal, halfuniform
stopifnot(!(((betahat=='abs_beta')|(betahat=='var_exp'))&(mixcompdist!='+uniform'))) # assert that mixcompdist must be +uniform if fitting on abs(beta) or variance explained

read_ss <- function(fname) {
  df = read.csv(file = fname, sep = '\t')
  stopifnot('se' %in% colnames(df))
  stopifnot(('A1' %in% colnames(df))&('A2' %in% colnames(df)))
  if (betahat=='var_exp') {
    stopifnot(('var_exp' %in% colnames(df))&('var_exp_se' %in% colnames(df)))
  }
  return(df)
}

main <- function(phen) {
  cat(paste0('Running ash for ',phen,'\n'))

  pruned_fname = sprintf('%s/pruned_gwas.%s.ld_wind_kb_%s%s%s.tsv.gz', data_wd, phen, 
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
                             data_wd, phen, 
                             if (block_mhc) '.block_mhc' else '',
                             if (is.null(maf)) '' else sprintf('.maf_%s',maf))
  genome.wide = read_ss(fname=genomewide_fname)
  
  cat(paste('...Starting ash on genome-wide file...\n'))
  cat(sprintf('block_mhc=%s, mixcompdist=%s, betahat=%s, pointmass=%s\n', block_mhc, mixcompdist, betahat, pointmass))

  genome.wide.ash = NULL
  
  for (chrom in seq(1,22)) {
    df = genome.wide[which(genome.wide$chr==chrom),]
    
    ash_betahat=if (grepl('abs_',betahat)) abs(df[,gsub("abs_","",betahat)])  else df[,betahat]
    ash_sebetahat=if (grepl('beta',betahat)) df$se else df[,paste0(betahat,'_se')]
    
    df.ash <- ash(ash_betahat, ash_sebetahat, 
                  g=ld.pruned.g, 
                  fixg=TRUE, 
                  mixcompdist = mixcompdist)
    if (is.null(genome.wide.ash)) {
      # print(df.ash$result)
      genome.wide.ash = df.ash$result
    }
    else {
      genome.wide.ash = rbind(genome.wide.ash, df.ash$result)
    }
    cat(paste('...chrom',chrom,'complete...\n'))
  }
  cat(paste('All chromosomes complete for',phen,'\n'))
  ash <- cbind(genome.wide[,c('chr','pos', 'A1', 'A2')], genome.wide.ash[,c('PosteriorMean','PosteriorSD')])
  fname_out = sprintf('%s/ash.%s.%s%s%s%s%s.tsv.gz', 
                      data_wd, 
                      phen,
                      mixcompdist, 
                      if (block_mhc) '.block_mhc' else '',
                      if (betahat!='beta') paste0('.',betahat) else '',
                      if (!pointmass) '.no_pointmass' else '',
                      if (is.null(maf)) '' else sprintf('.maf_%s', maf))

  write.table(x = ash, file = gzfile(fname_out), sep='\t', quote=F, row.names = FALSE)
  
}

if (length(args)==0) {
  for ( phen in phens ){
    main(phen = phen)
  }
} else {
  main(phen=args[1])
}

# pdf('/Users/nbaya/Downloads/test.pdf')
# for (i in 1:10) {
#   plot(ggplot(mapping = aes(x=c(1,2), y=c(1, 1)))+
#   geom_line())
# }
# dev.off()