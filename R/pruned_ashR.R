library(ashr)

data_wd = '/Users/nbaya/Documents/lab/smiles/data'

phenos = c(
  # 'standing_height',
  'bmi',
  't2d_bmiadj',
  'ibd',
  'cd',
  'uc',
  'scz',
  'ad')
ld_wind_kb = 500

for ( pheno in phenos ) {
  cat(paste0('Running ash for ',pheno,'\n'))
  cat(paste0('...Reading in genome-wide file...\n'))
  genomewide_fname = paste0('gwas.',pheno,'.tsv.gz')
  genome.wide = read.csv(file = paste0(data_wd,'/',genomewide_fname), sep = '\t')
  cat(paste0('...Reading in genome-wide file complete...\n'))
  pruned_fname = paste0('pruned_gwas.',pheno,'.ld_wind_kb_',ld_wind_kb,'.tsv.gz')
  ld.pruned = read.csv(file = paste0(data_wd,'/',pruned_fname), sep = '\t')
  
  ld.pruned.ash <- ash(ld.pruned$beta, ld.pruned$se, mixcompdist = "halfuniform")
  ld.pruned.g <- get_fitted_g(ld.pruned.ash)
  
  cat(paste('...Starting ash on genome-wide file...\n'))

  genome.wide.ash = NULL
  
  for (chrom in seq(22,1)) {
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
  ash <- cbind(genome.wide[,c('chr','pos','eaf')], genome.wide.ash[,c('PosteriorMean','PosteriorSD')])
  fname_out = gzfile(paste0(data_wd,'/ash.',pheno,'.tsv.gz'))
  write.table(x = ash, file = fname_out, sep='\t')
  
}