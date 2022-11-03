simulate_poisson <- function(num_cells = 1000,
                             num_genes = 100,
                             freq_cond = 0.5,
                             effsize_normaldist_mean = 0,
                             effsize_normaldist_sd = 1,
                             intercept_normaldist_mean = 0,
                             intercept_normaldist_sd = 1) {


  #Check that input values are reasonable
  #Make sure there is some variance in condition
  stopifnot(freq_cond > 0)

  #Restrict distribution parameters to permitted values
  stopifnot(effsize_normaldist_sd >= 0)
  stopifnot(intercept_normaldist_sd >= 0)

  #Simulate the independent variable

  indep_var = sample(c(0,1),
                     num_cells,
                     replace = TRUE,
                     prob = c(freq_cond, 1 - freq_cond))

  #Simulate gene-level parameters
  effsizes = rnorm(num_genes, effsize_normaldist_mean, effsize_normaldist_sd)
  intercepts = rnorm(num_genes, intercept_normaldist_mean, intercept_normaldist_sd)


  #Simulate count level data
  counts = matrix(data = NA, nrow = num_cells, ncol = num_genes)

  for (gene in c(1:num_genes)) {
    counts[,gene] = rpois(num_cells,
                          exp(effsizes[gene]*indep_var + intercepts[gene] * rep(1, num_cells)))
  }

  #return simulated data

  return(list(param = list(effsize_normaldist_mean = effsize_normaldist_mean,
                           effsize_normaldist_sd = effsize_normaldist_sd,
                           intercept_lognormaldist_mean = intercept_normaldist_mean,
                           intercept_lognormaldist_sd = intercept_normaldist_sd),
              condition = indep_var,
              effsizes = effsizes,
              intercepts = intercepts,
              counts = counts))
}
