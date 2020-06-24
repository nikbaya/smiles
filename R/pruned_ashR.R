library(ashr)

data_wd = '/Users/nbaya/Documents/lab/smiles/data'

phenos = c(
  # 'standing_height',
  # 'bmi',
  # 't2d_bmiadj',
  # 'ibd',
  # 'cd',
  'uc'
  # 'scz',
  # 'ad'
  )

ld_wind_kb = 500
block_mhc = F
mixcompdist = "halfnormal"

read_ss <- function(fname) {
  df = read.csv(file = paste0(data_wd,'/',fname), sep = '\t')
  stopifnot('se' %in% colnames(df))
  return(df)
}

for ( pheno in phenos ) {
  cat(paste0('Running ash for ',pheno,'\n'))
  cat(paste0('...Reading in genome-wide file...\n'))
  genomewide_fname = paste0('gwas.',pheno,'.tsv.gz')
  genome.wide = read_ss(fname=genomewide_fname)
  cat(paste0('...Reading in genome-wide file complete...\n'))
  if (block_mhc) {
    pruned_fname = paste0('pruned_gwas.',pheno,'.ld_wind_kb_',ld_wind_kb,'.block_mhc.tsv.gz')
  }
  else {
    pruned_fname = paste0('pruned_gwas.',pheno,'.ld_wind_kb_',ld_wind_kb,'.tsv.gz')
  }
  
  ld.pruned = read_ss(fname=pruned_fname)
  
  ld.pruned.ash <- ash(ld.pruned$beta, ld.pruned$se, mixcompdist = mixcompdist)
  ld.pruned.g <- get_fitted_g(ld.pruned.ash)
  
  cat(paste('...Starting ash on genome-wide file...\n'))
  cat(paste('mixcompdist = ',mixcompdist,'\n'))

  genome.wide.ash = NULL
  
  for (chrom in seq(1,22)) {
    df = genome.wide[which(genome.wide$chr==chrom),]
    df.ash <- ash(df$beta, df$se, mixcompdist = "halfuniform",
                           g=ld.pruned.g, fixg=TRUE)
    if (is.null(genome.wide.ash)) {
      genome.wide.ash = df.ash$result
    }
    else {
      genome.wide.ash = rbind(genome.wide.ash, df.ash$result)
    }
    cat(paste('...chrom',chrom,'complete...\n'))
  }
  cat(paste('All chromosomes complete for',pheno,'\n'))
  ash <- cbind(genome.wide[,c('chr','pos')], genome.wide.ash[,c('PosteriorMean','PosteriorSD')])
  if ( block_mhc) {
    fname_out = gzfile(paste0(data_wd,'/ash.',pheno,'.',mixcompdist,'.block_mhc.tsv.gz'))
  }
  else {
    fname_out = gzfile(paste0(data_wd,'/ash.',pheno,'.',mixcompdist,'.tsv.gz'))
  }
  write.table(x = ash, file = fname_out, sep='\t', quote=F)
  
}