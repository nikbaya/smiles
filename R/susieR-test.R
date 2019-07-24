library(susieR)
set.seed(1)


### TUTORIAL ###

## The data-set ##

data(N3finemapping)
attach(N3finemapping)

dim(data$Y)

b <- data$true_coef[,1]
plot(b, pch=16, ylab='effect size')

which(b != 0)

## Summary statistics from simple regression ##

z_scores <- sumstats[1,,] / sumstats[2,,]
z_scores <- z_scores[,1]
susie_plot(z_scores, y = "z", b=b)

R <- cor(data$X)

## Fine-mapping with `susieR` using summary statistics ##

# Using beta-hat and SE(beta-hat)

# WITH knowledge of variance of y
fitted_bhat <- susie_bhat(bhat = sumstats[1,,1], 
                          shat = sumstats[2,,1], 
                          R = R, n = nrow(data$X), 
                          var_y = var(data$Y[,1]),
                          L = 10, 
                          scaled_prior_variance = 0.1, 
                          estimate_residual_variance = TRUE, 
                          estimate_prior_variance = FALSE, 
                          standardize = TRUE)

# WITHOUT knowledge of variance of y
fitted_bhat_standardize <- susie_bhat(bhat = sumstats[1,,1], 
                                      shat = sumstats[2,,1], 
                                      R = R, n = nrow(data$X), 
                                      L = 10, 
                                      scaled_prior_variance = 0.1, 
                                      estimate_residual_variance = TRUE,
                                      estimate_prior_variance = FALSE)

summary(fitted_bhat)$cs

#-----------------------------------------------------------------------


## Test with real sumstats



