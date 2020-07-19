library(ashr)
library(truncnorm)
library(ggplot2)

smiles_dir= '/Users/nbaya/Documents/lab/smiles'
 
if (!dir.exists(smiles_dir)) {
  smiles_dir = '/stanley/genetics/users/nbaya/smiles'
}
data_dir=paste0(smiles_dir,'/data')

phenos = c(
  # 'standing_height',
  # 'bmi',
  # 't2d_bmiadj',
  'ibd',
  'cd',
  'uc',
  'scz',
  'ad',
  'breast_cancer',
  'cad',
  'ldl',
  # 'hdl'
  'wbc_count',
  'rbc_count',
  'urate',
  'systolic_bp',
  'diastolic_bp',
  'triglycerides'
)

ld_wind_kb = 500
block_mhc = TRUE
betahat = 'var_exp' # Which "betahat" is used for ash. options: beta, abs_beta, var_exp (default: beta)
pointmass=T # whether to include point mass at zero in prior
mixcompdist = '+uniform' #"+uniform" # options: normal, +uniform (if using var_exp instead of beta), halfnormal, halfuniform
stopifnot(!(((betahat=='abs_beta')|(betahat=='var_exp'))&(mixcompdist!='+uniform'))) # assert that mixcompdist must be +uniform if fitting on abs(beta) or variance explained

read_ss <- function(fname) {
  df = read.csv(file = fname, sep = '\t')
  stopifnot('se' %in% colnames(df))
  if (betahat=='var_exp') {
    stopifnot(('var_exp' %in% colnames(df))&('var_exp_se' %in% colnames(df)))
  }
  return(df)
}

for ( pheno in phenos ) {
  cat(paste0('Running ash for ',pheno,'\n'))
  
  pruned_fname = sprintf('%s/pruned_gwas.%s.ld_wind_kb_%s%s.tsv.gz', data_dir, pheno, 
                         ld_wind_kb, if (block_mhc) '.block_mhc' else '')
  ld.pruned = read_ss(fname=pruned_fname)
  
  ash_betahat=if (grepl('abs_',betahat)) abs(ld.pruned[,gsub("abs_","",betahat)])  else ld.pruned[,betahat]
  ash_sebetahat=if (grepl('beta',betahat)) ld.pruned$se else ld.pruned[,paste0(betahat,'_se')]
  
  ld.pruned.ash <- ash.workhorse(ash_betahat, ash_sebetahat, 
                                 mixcompdist = mixcompdist,
                                 pointmass=pointmass)
  ld.pruned.g <- get_fitted_g(ld.pruned.ash)
  
  if (grepl('normal', mixcompdist)) { 
    if (grepl('half', mixcompdist)) {
      g = data.frame(ld.pruned.g[1:5])   
    } else {
      g = data.frame(ld.pruned.g[1:3])   
    }
    x <- seq(min(g$mean-3*g$sd), max(g$mean+3*g$sd), length=1000)
    if (!(0 %in% x)&(pointmass)) { # in case of point mass at zero, add x=0 entry
      x = append(x, 0)
      x = sort(x)
    }
    y <- rep(0, length(x))
    density <-data.frame('x'=x, 'y'=y)
    for (row_idx in 1:nrow(g)) {
      row = g[row_idx,]
      if (grepl('half', mixcompdist)) {
        density$y = density$y+row$pi*dtruncnorm(density$x, a=row$a, b=row$b, mean=row$mean, sd=row$sd)
      } else {
        density$y = density$y+row$pi*dtruncnorm(density$x, a=row$mean-3*row$sd, b=row$mean+3*row$sd)
      }
    }
  } else if (grepl('uniform', mixcompdist)) {
    g = data.frame(ld.pruned.g[1:3])  
    x <- seq(0, max(g$b), length=1000)
    # if (!(0 %in% x)) { # in case of point mass at zero, add x=0 entry
    #   x = append(x, 0)
    #   x = sort(x)
    # }
    y <- rep(0, length(x))
    density <-data.frame('x'=x, 'y'=y)
    for (row_idx in 1:nrow(g)) {
      row = g[row_idx,]
      if (row$pi>0) {
        if (row$a!=row$b) { # ignore pointmass entries (i.e. distributions with no width)
          # density$y = density$y+row$pi*dunif(density$x, min=row$a, max=row$b)*(row$b-row$a) # max density of uniform function is always 1
          density$y = density$y+row$pi*dunif(density$x, min=row$a, max=row$b) # density integrates to 1
        }
        
      }
    }
    density <- rbind(data.frame(x=0,y=0), density, data.frame(x=max(density$x),y=0))
  }
  
  plot <- ggplot(density, aes(x=x, y=y)) +
    geom_line() +
    xlab(betahat) +
    ylab('density') +
    labs(title=sprintf('%s prior for %s', betahat, pheno),
         subtitle=sprintf('(mixcompdist=%s, pointmass=%s)', mixcompdist, pointmass)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.margin = margin(10, 30, 10, 10))
  
  if (pointmass) {
    row = g[g$a==g$b,]
    if (row$pi>0) {
      plot = plot + 
        geom_segment(aes(x=0, y=0, xend=0, yend=row$pi), color='red') +
        geom_point(aes(x=0, y=row$pi), color='red')
    }
  }
  
  print(plot)
  plots_dir = '/Users/nbaya/Downloads' # paste0(smiles_dir,'/plots')
  plot_fname=sprintf('%s/ashprior.%s.%s.%s%s%s.png',
                     plots_dir,
                     pheno,
                     betahat,
                     mixcompdist,
                     if (block_mhc) '.block_mhc' else '',
                     if (pointmass) '.pointmass' else '')
                     
  ggsave(filename = plot_fname,
         plot=plot,
         dpi=300,
         scale=0.8)

}
