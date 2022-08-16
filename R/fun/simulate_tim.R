simulate_tim <- function(num_cells = 1000,
                         num_genes = 100,
                         freq_cond = 0.5,
                         prob_betadist_alpha = 2,
                         prob_betadist_beta = 2,
                         effsize_normaldist_mean = 0,
                         effsize_normaldist_sd = 1,
                         intercept_lognormaldist_mean = 0.2,
                         intercept_lognormaldist_sd = 0.8) {
  
  #Simulate the independent variable
  
  indep_var = sample(c(0,1),
                     num_cells,
                     replace = TRUE,
                     prob = c(freq_cond, 1 - freq_cond))
  
  #Simulate gene-level parameters
  probs = rbeta(num_genes, prob_betadist_alpha, prob_betadist_beta)
  effsizes = rnorm(num_genes, effsize_normaldist_mean, effsize_normaldist_sd)
  intercepts = rlnorm(num_genes, intercept_lognormaldist_mean, intercept_lognormaldist_sd)
  
  
  #Simulate count level data
  counts = matrix(data = NA, nrow = num_cells, ncol = num_genes)
  
  for (gene in c(1:num_genes)) {
    
    #Convert to alternate parameterization for stability
    counts[,gene] = sapply(1:num_cells,
                           function(cell) {
                             rnbinom(1,
                                     2^(effsizes[gene]*indep_var[cell] + intercepts[gene]),
                                     probs[gene])
                           })
  }
  
  #return simulated data
  
  return(list(param = list(prob_betadist_alpha = prob_betadist_alpha,
                           prob_betadist_beta = prob_betadist_beta,
                           effsize_normaldist_mean = effsize_normaldist_mean,
                           effsize_normaldist_sd = effsize_normaldist_sd,
                           intercept_lognormaldist_mean = intercept_lognormaldist_mean,
                           intercept_lognormaldist_sd = intercept_lognormaldist_sd),
              probs = probs,
              effsizes = effsizes,
              intercepts = intercepts,
              counts = counts))
}