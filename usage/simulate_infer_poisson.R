
effsize_sd = as.numeric(commandArgs(trailingOnly=TRUE)[1])

simulate_infer_poisson = function(effsize_sd){
  print(effsize_sd)
  counts = simulate_poisson(effsize_normaldist_sd = effsize_sd)
  inferred_params = infer_poisson(counts)

  return(inferred_params$effsize_sd)

}

inferred_effsize_sd = simulate_infer_poisson(as.numeric(commandArgs(trailingOnly=TRUE)[1]))

write.table(inferred_effsize_sd, paste0("/broad/oconnor/ajay/tim/usage/poisson_sim_output/",commandArgs(trailingOnly=TRUE)[1],".txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
