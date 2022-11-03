infer_poisson <- function(counts,
                          init_effsize_normaldist_mean = 0,
                          init_effsize_normaldist_sd = 3,
                          init_intercept_normaldist_mean = 0,
                          init_intercept_normaldist_sd = 5,
                          iterations = 100,
                          step_size = 0.005) {


  betas_init = get_betas(counts,
                    init_effsize_normaldist_mean,
                    init_effsize_normaldist_sd,
                    init_intercept_normaldist_mean,
                    init_intercept_normaldist_sd)

  iter = 0
  betas_curr = betas_init
  curr_effsize_normaldist_mean = init_effsize_normaldist_mean
  curr_effsize_normaldist_sd = init_effsize_normaldist_sd
  curr_intercept_normaldist_mean = init_intercept_normaldist_mean
  curr_intercept_normaldist_sd = init_intercept_normaldist_sd

  loglikelihood_vals = rep(NA, iterations)

  while (iter < iterations){
    iter = iter + 1
    loglikelihood_current = get_loglikelihood(curr_effsize_normaldist_mean,
                                              curr_effsize_normaldist_sd,
                                              curr_intercept_normaldist_mean,
                                              curr_intercept_normaldist_sd,
                                              counts,
                                              betas_curr)

    loglikelihood_vals[iter] = loglikelihood_current

    print(iter)
    print(loglikelihood_current)
    gradient = (c(get_loglikelihood(curr_effsize_normaldist_mean + 0.001,
                                   curr_effsize_normaldist_sd,
                                   curr_intercept_normaldist_mean,
                                   curr_intercept_normaldist_sd,
                                   counts,
                                   betas_curr),
                 get_loglikelihood(curr_effsize_normaldist_mean,
                                   curr_effsize_normaldist_sd + 0.001,
                                   curr_intercept_normaldist_mean,
                                   curr_intercept_normaldist_sd,
                                   counts,
                                   betas_curr),
              get_loglikelihood(curr_effsize_normaldist_mean,
                                curr_effsize_normaldist_sd,
                                curr_intercept_normaldist_mean + 0.001,
                                curr_intercept_normaldist_sd,
                                counts,
                                betas_curr),
              get_loglikelihood(curr_effsize_normaldist_mean,
                                curr_effsize_normaldist_sd,
                                curr_intercept_normaldist_mean,
                                curr_intercept_normaldist_sd + 0.001,
                                counts,
                                betas_curr)) - rep(loglikelihood_current, 4))/0.001

    curr_effsize_normaldist_mean = curr_effsize_normaldist_mean + gradient[1]*step_size
    curr_effsize_normaldist_sd = curr_effsize_normaldist_sd + gradient[2]*step_size
    curr_intercept_normaldist_mean = curr_intercept_normaldist_mean + gradient[3]*step_size
    curr_intercept_normaldist_sd = curr_intercept_normaldist_sd + gradient[4]*step_size

    betas_curr = get_betas(counts,
                           curr_effsize_normaldist_mean,
                           curr_effsize_normaldist_sd,
                           curr_intercept_normaldist_mean,
                           curr_intercept_normaldist_sd)

  }

  return(list(effsize_mean = curr_effsize_normaldist_mean,
           effsize_sd = curr_effsize_normaldist_sd,
           intercept_mean = curr_intercept_normaldist_mean,
           intercept_sd = curr_intercept_normaldist_sd,
           betas = betas_curr))

}

get_betas <- function(counts,
                      effsize_normaldist_mean,
                      effsize_normaldist_sd,
                      intercept_normaldist_mean,
                      intercept_normaldist_sd) {

  betas = matrix(data = NA, nrow = ncol(counts$counts), ncol = 2)

  for (gene in 1:ncol(counts$counts)){
    y = counts$counts[,gene]
    x = counts$condition

    data = data.frame(indep_var = x,
                      counts = y)

    model = stan_glm(counts ~ indep_var,
                     data,
                     family = poisson(link = "log"),
                     prior = normal(effsize_normaldist_mean,effsize_normaldist_sd),
                     prior_intercept = normal(intercept_normaldist_mean,intercept_normaldist_sd),
                     chains = 1,
                     iter = 50, seed = 84735,refresh = 0)
    betas[gene,] = model$coefficients
  }

  return(betas)

}

get_loglikelihood <- function(curr_effsize_normaldist_mean,
                              curr_effsize_normaldist_sd,
                              curr_intercept_effsize_normaldist_mean,
                              curr_intercept_effsize_normaldist_sd,
                              counts,
                              betas) {
  loglikelihood.df = data.frame(effsize_loglikelihoods = log(dnorm(betas[,2],
                                                                   mean = curr_effsize_normaldist_mean,
                                                                   sd = curr_effsize_normaldist_sd)),
                                intercept_loglikelihoods = log(dnorm(betas[,1],
                                                                     mean = curr_intercept_effsize_normaldist_mean,
                                                                     sd = curr_intercept_effsize_normaldist_sd)),
                                count_loglikelihood = colSums(log(sapply(1:ncol(counts$counts),
                                                                         function(x) {
                                                                           dpois(counts$counts[,x],
                                                                                 exp(betas[x,2]*counts$condition + betas[x,1]*rep(1,length(counts$condition))))
                                                                         }))))

  loglikelihood = sum(loglikelihood.df)
  return(loglikelihood)
}

testfunc <- function(n_iter){
  betas = matrix(data = NA, nrow = ncol(counts$counts), ncol = 2)

  for (gene in 1:ncol(counts$counts)){
    data = data.frame(indep_var = counts$condition,
                      counts = counts$counts[,gene])

    model = stan_glm(counts ~ indep_var,
                     data,
                     family = poisson(link = "log"),
                     prior = normal(0,1),
                     prior_intercept = normal(4,7),
                     chains = 1,
                     iter = n_iter, seed = 84735,refresh = 0)
    betas[gene,] = model$coefficients
  }

  RMSE = sqrt(sum((counts$effsizes- betas[,2])^2)/length(betas[,2]))
  return(RMSE)
}



