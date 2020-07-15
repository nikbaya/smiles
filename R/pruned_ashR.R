library(ashr)

data_wd = '/Users/nbaya/Documents/lab/smiles/data'
if (!dir.exists(data_wd)) {
  data_wd = '/stanley/genetics/users/nbaya/smiles/data'
}

phenos = c(
  'standing_height',
  # 'bmi',
  't2d_bmiadj',
  # 'ibd',
  # 'cd',
  # 'uc',
  # 'scz',
  # 'ad',
  # 'breast_cancer'
  # 'cad',
  # 'ldl',
  'hdl'
  # 'wbc_count',
  # 'rbc_count',
  # 'urate',
  # 'systolic_bp',
  # 'diastolic_bp',
  # 'triglycerides'
  )

ld_wind_kb = 500
block_mhc = TRUE
betahat = 'abs_beta' # Which "betahat" is used for ash. options: beta, abs_beta, var_exp (default: beta)
pointmass=TRUE # whether to include point mass at zero in prior
mixcompdist = '+uniform' #"+uniform" # options: normal, +uniform (if using var_exp instead of beta), halfnormal, halfuniform
stopifnot(!(((betahat=='abs_beta')|(betahat=='var_exp'))&(mixcompdist!='+uniform'))) # assert that mixcompdist must be +uniform if fitting on abs(beta) or variance explained

read_ss <- function(fname) {
  df = read.csv(file = fname, sep = '\t')
  stopifnot('se' %in% colnames(df))
  if (use_var_exp) {
    stopifnot(('var_exp' %in% colnames(df))&('var_exp_se' %in% colnames(df)))
  }
  return(df)
}

for ( pheno in phenos ) {
  cat(paste0('Running ash for ',pheno,'\n'))

  pruned_fname = sprintf('%s/pruned_gwas.%s.ld_wind_kb_%s%s.tsv.gz', data_wd, pheno, 
                         ld_wind_kb, if (block_mhc) '.block_mhc' else '')
  ld.pruned = read_ss(fname=pruned_fname)
  
  ash_betahat=if (grepl('abs_',betahat)) abs(ld.pruned[,gsub("abs_","",betahat)])  else ld.pruned[,betahat]
  ash_sebetahat=if (grepl('beta',betahat)) ld.pruned$se else ld.pruned[,paste0(betahat,'_se')]
  
  ld.pruned.ash <- ash.workhorse(ash_betahat, ash_sebetahat, 
                       mixcompdist = mixcompdist,
                       pointmass=pointmass)
  ld.pruned.g <- get_fitted_g(ld.pruned.ash)
  
  cat(paste0('...Reading in genome-wide file...\n'))
  genomewide_fname = sprintf('%s/gwas.%s%s.tsv.gz', 
                             data_wd, pheno, 
                             if (block_mhc) '.block_mhc' else '')
  genome.wide = read_ss(fname=genomewide_fname)
  
  cat(paste('...Starting ash on genome-wide file...\n'))
  cat(sprintf('mixcompdist=%s, block_mhc=%s, pointmass=%s\n', mixcompdist, block_mhc, pointmass))

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
      genome.wide.ash = df.ash$result
    }
    else {
      genome.wide.ash = rbind(genome.wide.ash, df.ash$result)
    }
    cat(paste('...chrom',chrom,'complete...\n'))
  }
  cat(paste('All chromosomes complete for',pheno,'\n'))
  ash <- cbind(genome.wide[,c('chr','pos', 'A1', 'A2')], genome.wide.ash[,c('PosteriorMean','PosteriorSD')])
  fname_out = sprintf('%s/ash.%s.%s%s%s.tsv.gz', data_wd, 
                      pheno, mixcompdist, 
                      if (block_mhc) '.block_mhc' else '',
                      if (betahat!='beta') paste0('.',betahat) else '',
                      if (!pointmass) '.no_pointmass' else '')

  write.table(x = ash, file = gzfile(fname_out), sep='\t', quote=F)
  
}