
simulation_params = data.frame(effsize_mean = 0,
                               effsize_sd = seq(0.2,2,length.out = 10))

simulate_infer_poisson = function(effsize_sd){
  print(effsize_sd)
  counts = simulate_poisson(effsize_normaldist_sd = effsize_sd)
  inferred_params = infer_poisson(counts)

  return(inferred_params$effsize_sd)

}


inferred_effsize_sd = sapply(seq(0.5,2,length.out = 10),
                             simulate_infer_poisson)
