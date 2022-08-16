
test_that("Distributional parameters are within permitted range", {
  expect_error(simulate_tim(num_cells = 1000,
                            num_genes = 500,
                            freq_cond = 0.5,
                            prob_betadist_alpha = -2,
                            prob_betadist_beta = 2,
                            effsize_normaldist_mean = 1,
                            effsize_normaldist_sd = 1,
                            intercept_lognormaldist_mean = 0.2,
                            intercept_lognormaldist_sd = 0.8),
               "prob_betadist_alpha > 0 is not TRUE")
  expect_error(simulate_tim(num_cells = 1000,
                            num_genes = 100,
                            freq_cond = 0.5,
                            prob_betadist_alpha = 2,
                            prob_betadist_beta = -2,
                            effsize_normaldist_mean = 1,
                            effsize_normaldist_sd = 1,
                            intercept_lognormaldist_mean = 0.2,
                            intercept_lognormaldist_sd = 0.8),
               "prob_betadist_beta > 0 is not TRUE")
  expect_error(simulate_tim(num_cells = 1000,
                            num_genes = 100,
                            freq_cond = 0.5,
                            prob_betadist_alpha = 2,
                            prob_betadist_beta = 2,
                            effsize_normaldist_mean = 1,
                            effsize_normaldist_sd = -1,
                            intercept_lognormaldist_mean = 0.2,
                            intercept_lognormaldist_sd = 0.8),
               "effsize_normaldist_sd >= 0 is not TRUE")
  expect_error(simulate_tim(num_cells = 1000,
                            num_genes = 100,
                            freq_cond = 0.5,
                            prob_betadist_alpha = 2,
                            prob_betadist_beta = 2,
                            effsize_normaldist_mean = 1,
                            effsize_normaldist_sd = 1,
                            intercept_lognormaldist_mean = -0.2,
                            intercept_lognormaldist_sd = 0.8),
               "intercept_lognormaldist_mean > 0 is not TRUE")
  expect_error(simulate_tim(num_cells = 1000,
                            num_genes = 100,
                            freq_cond = 0.5,
                            prob_betadist_alpha = 2,
                            prob_betadist_beta = 2,
                            effsize_normaldist_mean = 1,
                            effsize_normaldist_sd = 1,
                            intercept_lognormaldist_mean = 0.2,
                            intercept_lognormaldist_sd = -0.8),
               "intercept_lognormaldist_sd >= 0 is not TRUE")
})

test_run = simulate_tim(num_cells = 1000,
                        num_genes = 500,
                        freq_cond = 0.5,
                        prob_betadist_alpha = 2,
                        prob_betadist_beta = 2,
                        effsize_normaldist_mean = 1,
                        effsize_normaldist_sd = 1,
                        intercept_lognormaldist_mean = 0.2,
                        intercept_lognormaldist_sd = 0.1)

test_that("NB probabilities follow Beta distribution", {
  expect_gte(mean(test_run$probs), (2/(2+2)) - (0.1 * (2/(2+2))))
  expect_lte(mean(test_run$probs), (2/(2+2)) + (0.1 * (2/(2+2))))
  expect_gte(var(test_run$probs), (2*2/(((2+2)^2)*(2 + 2 + 1))) - 0.1*(2*2/(((2+2)^2)*(2 + 2 + 1))))
  expect_lte(var(test_run$probs), (2*2/(((2+2)^2)*(2 + 2 + 1))) + 0.1*(2*2/(((2+2)^2)*(2 + 2 + 1))))
})

test_that("Effect sizes follow Normal distribution", {
  expect_gte(mean(test_run$effsizes), 1 - 1*0.1)
  expect_lte(mean(test_run$effsizes), 1 + 1*0.1)
  expect_gte(var(test_run$effsizes), 1 - 1*0.1)
  expect_lte(var(test_run$effsizes), 1 + 1*0.1)
})

test_that("Intercepts follow Log Normal distribution", {
  expect_gte(mean(log(test_run$intercepts)), 0.2 - 0.2*0.1)
  expect_lte(mean(log(test_run$intercepts)), 0.2 + 0.2*0.1)
  expect_gte(var(log(test_run$intercepts)), 0.1^2 - (0.1^2)*0.1)
  expect_lte(var(log(test_run$intercepts)), 0.1^2 + (0.1^2)*0.1)
})

test_that("Per-gene count means are correlated with distribution expectation", {
  count_means_0 = colMeans(test_run$counts[test_run$condition == 0,])
  expected_means_0 = (test_run$probs * (2^(test_run$effsizes*0 + test_run$intercepts)))/(1 - test_run$probs)

  expect_gte(cor(count_means_0,expected_means_0), 0.99)

  count_means_1 = colMeans(test_run$counts[test_run$condition == 1,])
  expected_means_1 = (test_run$probs * (2^(test_run$effsizes*1 + test_run$intercepts)))/(1 - test_run$probs)

  expect_gte(cor(count_means_1,expected_means_1), 0.99)
})
