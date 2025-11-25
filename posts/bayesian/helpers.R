r_beta_bernoulli = function(n, alpha, beta) {
  probs = rbeta(n, alpha, beta)
  sample = rbinom(n, 1, probs)
  return(sample)
}

r_beta_binom = function(n, size, alpha, beta) {
  probs = rbeta(n, alpha, beta)
  sample = rbinom(n, size, probs)
  return(sample)
}

r_markov_bernoulli = function(n, p0, drift_sd) {
  probs = numeric(n)
  probs[1] = p0
  
  for (i in 2:n) {
    probs[i] = probs[i - 1] + rnorm(1, 0, drift_sd)
    probs[i] = min(max(probs[i], 0), 1)
  }
  sample = rbinom(n, 1, probs)
  return(sample)
}

r_markov_counts = function(n, size, p0, drift_sd) {
  sample = numeric(n)
  for (i in 1:n) {
    sample[i] = sum(r_markov_bernoulli(size, p0[i], drift_sd[i]))
  }
  return(sample)
}

pp_sample_distance = function(pp_sample1, pp_sample2, method = "tvd") {
  
  if (length(pp_sample1) != length(pp_sample2)) {
    print("warning: sample sizes are different")
  }
  
  if (method == "tvd") {
    max_val = max(c(pp_sample1, pp_sample2))
    
    counts1 = tabulate(pp_sample1, nbins = max_val)
    counts2 = tabulate(pp_sample2, nbins = max_val)
    
    counts1 = counts1 / sum(counts1)
    counts2 = counts2 / sum(counts2)
    
    return(0.5 * sum(abs(counts1 - counts2)))
  } else if (method == "binned_tvd") {
    combined = c(pp_sample1, pp_sample2)
    breaks = hist(combined, plot = FALSE)$breaks
  
    counts1 = hist(pp_sample1, breaks = breaks, plot = FALSE)$counts
    counts2 = hist(pp_sample2, breaks = breaks, plot = FALSE)$counts
    
    counts1 = counts1 / sum(counts1)
    counts2 = counts2 / sum(counts2)
    
    return(0.5 * sum(abs(counts1 - counts2)))
  } else {
    print("available methods are 'tvd' and 'binned_tcd'")
  }
}

plot_pp_difference = function(real_sample, miss_sample, prior_or_post, real_model_name, miss_model_name) {
  
  plot_df = tibble(
   D = real_sample,
   M = miss_sample
  ) |> 
  pivot_longer(
    cols = D:M,
    names_to = "source",
    values_to = "values"
  )

  ggplot(plot_df, aes(x = values, fill = source)) +
    geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.5, bins = 30) +
    labs(title = glue("{prior_or_post} Predictive Distribution Comparison"), x = "Count", y = "Density") +
    scale_fill_discrete(
      labels = c(glue("Real Model ({real_model_name})"), glue("Misspecified Model ({miss_model_name})"))) + 
    theme_minimal()
}

sample_from_pp_D1 = function(n, size, mu1, sigma1, mu2, sigma2) {
  alphas = rlnorm(n, mu1, sigma1)
  betas = rlnorm(n, mu2, sigma2)
  probs = rbeta(n, alphas, betas)
  sample = rbinom(n, size, probs)
  return(sample)
}

sample_from_pp_D2 = function(n, size, alpha, beta, rate) {
    p0 = rbeta(n, alpha, beta)
    drift_sd = rexp(n, rate)
    sample = r_markov_counts(n, size, p0, drift_sd)
    return(sample)
}

r_markov_bernoulli_old = function(n, p_initial, p_after_1, p_after_0){
  sample = numeric(n)
  sample[1] = rbinom(1, 1, p_initial)
  for (i in 2:n) {
    if (sample[i - 1] == 1) {sample[i] = rbinom(1, 1, p_after_1)}
    else {sample[i] = rbinom(1, 1, p_after_0)}
  }
  return(sample)
}

r_markov_counts_old = function(n, size, p_initial, p_after_1, p_after_0) {
  sample = numeric(n)
  for (i in 1:n) {
    sample[i] = sum(r_markov_bernoulli_old(size, p_initial, p_after_1, p_after_0))
  }
  return(sample)
}
