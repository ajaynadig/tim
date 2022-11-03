tim_likelihood <- function(data,
                           effsizes,
                           intercepts,
                           probs = NULL,
                           prob_betadist_alpha,
                           prob_betadist_beta,
                           effsize_normaldist_mean,
                           effsize_normaldist_sd,
                           intercept_lognormaldist_mean,
                           intercept_lognormaldist_sd) {

  #Calculate likelihoods of effect sizes and intercepts
  effsize_likelihoods = dnorm(effsizes,
                              effsize_normaldist_mean,
                              effsize_normaldist_sd)

  intercept_likelihoods = dlnorm(intercepts,
                                 intercept_lognormaldist_mean,
                                 intercept_lognormaldist_sd)


  gene_likelihoods <- rep(NA, ncol(data))

  for (gene in 1:ncol(data)){
    effsize_likelihood = effsize_likelihoods[gene]
    intercept_likelihood = intercept_likelihoods[gene]


  }
}
