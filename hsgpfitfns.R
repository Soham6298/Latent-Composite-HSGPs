library(truncnorm)
# functions
## Sq exponential cov fn
se <- function(x, alpha, rho, delta) {
  K = matrix(nrow = length(x), ncol = length(x))
  N = length(x)
  for(i in 1:N) {
    K[i, i] = alpha^2 + delta
    if(i == N) {
      break
    }
    for(j in (i + 1):N) {
      K[i, j] = alpha^2 * exp(-(x[i] - x[j])^2 / (2*rho^2))
      K[j, i] = K[i, j]
    }
  }
  return(K)
}
# Matern 3/2
m32 <- function(x, alpha, rho, delta) {
  K = matrix(nrow = length(x), ncol = length(x))
  N = length(x)
  sq_alpha = alpha^2
  r = -1/rho
  
  for(i in 1:N) {
    K[i, i] = sq_alpha + delta
    if(i == N) {
      break
    }
    for(j in (i + 1):N) {
      K[i, j] = (1 - (r * sqrt(3) * abs(x[i] - x[j]))) * 
        exp(r * sqrt(3) * abs(x[i] - x[j])) * sq_alpha
      K[j, i] = K[i, j]
    }
  }
  return(K)
}

# Matern 5/2 
m52 <- function(x, alpha, rho, delta) {
  K = matrix(nrow = length(x), ncol = length(x))
  N = length(x)
  sq_rho = rho^2;
  rho4 = rho^4
  sq_alpha = alpha^2
  r = -1/(2 * sq_rho)
  
  for(i in 1:N) {
    K[i, i] = sq_alpha + delta
    if(i == N) {
      break
    }
    for(j in (i + 1):N) {
      K[i, j] = (1 + (sqrt(5) * abs(x[i] - x[j])/rho) + 
                   ((5 * (x[i] - x[j])^2)/ (3 * sq_rho))) * 
        exp(- sqrt(5) * abs(x[i] - x[j])/rho) * sq_alpha
      K[j, i] = K[i, j]
    }
  }
  return(K)
}


## Draw GP means from mvnorm
gp_draw <- function(draws, x, Sigma, ...) {
  mu <- rep(0, length(x))
  mvtnorm::rmvnorm(draws, mu, Sigma)
}
## Function copied from rethinking
rlkjcorr <- function ( n , K , eta = 1 ) {
  
  stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
  stopifnot(eta > 0)
  #if (K == 1) return(matrix(1, 1, 1))
  
  f <- function() {
    alpha <- eta + (K - 2)/2
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R <- matrix(0, K, K) # upper triangular Cholesky factor until return()
    R[1,1] <- 1
    R[1,2] <- r12
    R[2,2] <- sqrt(1 - r12^2)
    if(K > 2) for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m / 2, alpha)
      
      # Draw uniformally on a hypersphere
      z <- rnorm(m, 0, 1)
      z <- z / sqrt(crossprod(z)[1])
      
      R[1:m,m+1] <- sqrt(y) * z
      R[m+1,m+1] <- sqrt(1 - y)
    }
    return(crossprod(R))
  }
  R <- replicate( n , f() )
  if ( dim(R)[3]==1 ) {
    R <- R[,,1]
  } else {
    # need to move 3rd dimension to front, so conforms to array structure that Stan uses
    R <- aperm(R,c(3,1,2))
  }
  return(R)
}
sim_input <- function(n_obs, true_x, s_x){
  x_obs <- rnorm(n_obs, true_x, s_x)
  data.frame(x_true = true_x, x_obs)
}
gp_sim_data <- function(n_obs, dims, true_x, 
                        rho, alpha, sigma, intc, covfn,...){  # dims is no. of multioutput dimensions (>1)
  x_true <- true_x
  delta <- 1e-12
  K <- list()
  for (j in 1:dims){
    if(covfn == 'm32') {
      K[[j]] <- m32(x_true, alpha = alpha[j], rho = rho[j], delta) 
    } else if(covfn=='m52') {
      K[[j]] <- m52(x_true, alpha = alpha[j], rho = rho[j], delta)
    } else {
      K[[j]] <- se(x_true, alpha = alpha[j], rho = rho[j], delta)
    }
  }
  ### Draw from joint GP (using true_t)
  f0 <- matrix(nrow = dims, ncol = length(x_true))
  for( j in 1:dims){
    f0[j,] <- gp_draw(1, x_true, K[[j]])  # Computing GP mean
  } 
  C <- rlkjcorr(n = 1, K = dims, eta = 1)
 #C <- matrix(rep(corr,dims^2), ncol = dims) # Setting within dims output correlation
  #diag(C) <- 1
  f <- t(f0) %*% chol(C)     # GP mean with output dims correlations
  y <- matrix(ncol= ncol(f), nrow = nrow(f))
  for (j in 1:ncol(f)) {
    y[,j] <- rnorm((length(f[,j])), mean = intc[j] + f[,j], sd = sigma[j])
  }
  data.frame(y, f)
}

# For simulating multi-output periodic data
periodic_data <- function(n_obs, dims, true_x, s_x, rho,
                          alpha, sigma, intc) {
  x_true <- true_x   # Rep for obs and grad
  x_obs <- rnorm(n_obs, true_x, s_x)  # obs_t ~ N(true_t, s)
  f0 <- matrix(nrow = n_obs, ncol = dims)
  for(i in 1:n_obs){
    for(j in 1:dims){
      f0[i, j] <- alpha[j] * sin(x_true[i]/rho[j])
    }
  }
  C <- rlkjcorr(n = 1, K = dims, eta = 1)
  #C <- matrix(rep(corr,dims^2), ncol = dims) # Setting within dims output correlation
  #diag(C) <- 1
  f <- f0 %*% chol(C)     # GP mean with output dims correlations
  y <- matrix(nrow = n_obs, ncol = dims)
  for (j in 1:ncol(f)) {
    y[,j] <- rnorm(length(f[,j]), mean = intc[j] + f[,j], sd = sigma[j])        # Original output
  }
  data.frame(y, f, x_true, x_obs)
}

## Bias and RMSE for summary
abs_bias_draws <- function(theta_hat, theta) {
  abs(mean(theta_hat) - theta)
}
bias_draws <- function(theta_hat, theta) {
  (mean(theta_hat) - theta)
}
rmse_draws <- function(theta_hat, theta) {
  sqrt(mean((theta_hat - theta)^2))
}
mae_draws <- function(theta_hat, theta) {
  mean(abs(theta_hat - theta))
}

## For choosing boundary conditions L
choose_L <- function(x, c) {
  if (!length(x)) {
    range <- 1
  } else {
    range <- max(1, max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  c * range
}

## Summary table function
compare_summary <- function(model, variable, dims, true_variable, n_obs, m_approx, 
                            variable_class, model_name, sim_id) {
  out <- as_draws_matrix(model)
  subset <- subset_draws(out, variable = variable)
  summary <- summarise_draws(subset)
  mean <- summary$mean
  rhat <- summary$rhat
  bess <- summary$ess_bulk
  tess <- summary$ess_tail
  sd <- rep(NA, length(variable))
  abs_bias <- rep(NA, length(variable))
  rmse <- rep(NA, length(variable))
  mae <- rep(NA, length(variable))
  for(i in 1:length(variable)) {
    sd[i] <- sd(subset[,i])
    abs_bias[i] <- abs_bias_draws(subset[,i], true_variable[i])
    rmse[i] <- rmse_draws(subset[,i], true_variable[i])
    mae[i] <- mae_draws(subset[,i], true_variable[i])
  }
  true_value <- true_variable
  class <- rep(variable_class, length(variable))
  model_name <- rep(model_name, length(variable))
  n <- rep(n_obs, length(variable))
  m <- rep(m_approx, length(variable))
  d <- rep(dims, length(variable))
  pars <- variable
  runtime <- rep(max(rowSums(get_elapsed_time(model))), length(variable))
  sampler_params <- get_sampler_params(model)
  divergent_by_chain <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
  divergent <- rep(max(divergent_by_chain), length(variable))
  sim_id <- rep(sim_id, length(variable))
  ranks <- colSums(sweep(subset, 2, true_variable, "<"))
  data.frame(sim_id, n, m, d, class, pars, model_name, true_value, mean,
             sd, abs_bias, rmse, mae, ranks, rhat, bess, tess, runtime, divergent)
}

# Summary fucntion for python models
compare_summary_py_test <- function(latent_mean, latent_sd, true_value, dims, n_obs, m_approx, 
                                    variable_class, model_name, sim_no, sample_id, n_draws, 
                                    runtime){
  implied_samples <- rnorm(n_draws, mean = latent_mean, sd = latent_sd)
  rhat <- NA
  bess <- NA
  tess <- NA
  divergent <- NA
  n <- n_obs
  m <- m_approx
  d <- dims
  sd <- sd(implied_samples)
  abs_bias <- abs_bias_draws(implied_samples, true_value)
  rmse <- rmse_draws(implied_samples, true_value)
  mae <- mae_draws(implied_samples, true_value)
  mean <- mean(implied_samples)
  sim_id <- sim_no
  class <- variable_class
  pars <- sprintf('x[%s]', sample_id)
  ranks <- sum(implied_samples < true_value)
  data.frame(sim_id, n, m, d, class, pars, model_name, true_value, mean,
             sd, abs_bias, rmse, mae, ranks, rhat, bess, tess, runtime, divergent)
  
}

## model fit
gp_model <- function(gpmodel, n_obs, dims, outputs, inputs, latent_sd, latent_inputs, 
                     x_min, x_max,covfn, is_vary, is_corr, rho_prior, ls_param, msd_param, 
                     esd_param, mean_intc, sd_intc, iter, warmup, chains, cores, init, adapt_delta){
  ### Exact GP
  gp_data <- list(
    N = n_obs,
    D = dims,
    y = outputs,
    inputs = inputs,
    s = latent_sd,
    latent = latent_inputs,
    x_min = x_min,
    x_max = x_max,
    covfn = covfn,
    is_vary = is_vary,
    is_corr = is_corr,
    rho_prior = rho_prior,
    ls_param = ls_param,
    msd_param = msd_param,
    esd_param = esd_param,
    mean_intc = mean_intc,
    sd_intc = sd_intc
  )
  #### model fitting
  gpfit <- sampling(
    object = gpmodel,
    data   = gp_data,
    iter   = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    init = init,
    control = list(adapt_delta = adapt_delta)
  )
  return(gpfit)
}

hsgp_model <- function(hsgpmodel, 
                       n_obs, 
                       m_obs_f, 
                       m_obs_g,
                       dims, 
                       outputs1,
                       outputs2,
                       inputs,
                       latent_sd,
                       latent_inputs, 
                       x_min, 
                       x_max, 
                       is_vary, 
                       is_corr,
                       rho_prior, 
                       ls_param_f,
                       ls_param_g,
                       msd_param_f, 
                       msd_param_g,
                       esd_param_f,
                       esd_param_g, 
                       intc_y1, 
                       intc_y2,  
                       c, 
                       covfn,
                       iter, 
                       warmup, 
                       chains,
                       cores, 
                       init, 
                       adapt_delta){
  ### HSGP
  L <- choose_L(inputs, c = c) 
  hsgp_data <- list(
    L = L,
    N = n_obs,
    M_f = m_obs_f,
    M_g = m_obs_g,
    D = dims,
    y1 = outputs1,
    y2 = outputs2,
    inputs = inputs,
    s = latent_sd,
    covfn = covfn,
    latent = latent_inputs,
    x_min = x_min,
    x_max = x_max,
    is_vary = is_vary,
    is_corr = is_corr,
    rho_prior = rho_prior,
    ls_param_f = ls_param_f,
    ls_param_g = ls_param_g,
    msd_param_f = msd_param_f,
    msd_param_g = msd_param_g,
    esd_param_f = esd_param_f,
    esd_param_g = esd_param_g,
    intc_y1 = intc_y1,
    intc_y2 = intc_y2
  )
  #### model fitting
  hsgpfit <- sampling(
    object = hsgpmodel,
    data   = hsgp_data,
    iter   = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    init = init,
    control = list(adapt_delta = adapt_delta)
  )
  return(hsgpfit)
}

# Parameters for GP constant across outputs
true_params <- function(n_obs, dims, msd_pars, esd_pars, ls_pars){
  alpha <- rep(rtruncnorm(1, a = 0.1, b = Inf, mean = msd_pars[1], sd = msd_pars[2]), dims)        # GP marginal SD
  sigma <- rep(rtruncnorm(1, a = 0.1, b = Inf, mean = esd_pars[1], sd = esd_pars[2]), dims)        # Error SD 
  rho <- rep(rtruncnorm(1, a = 0.1, b = Inf, mean = ls_pars[1], sd = ls_pars[2]), dims)        # GP length scale
  data.frame(alpha, sigma, rho)
}

# Parameters for GP varying across outputs
true_params_vary <- function(n_obs, dims, msd_pars, esd_pars, ls_pars, intc_pars){
  alpha <- rtruncnorm(dims, a = 0.1, b = Inf, mean = msd_pars[1], sd = msd_pars[2])      # GP marginal SD
  sigma <- rtruncnorm(dims, a = 0.1, b = Inf, mean = esd_pars[1], sd = esd_pars[2])       # Error SD 
  rho <- rtruncnorm(dims, a = 0.1, b = Inf, mean = ls_pars[1], sd = ls_pars[2]) 
  intc <- rnorm(dims, mean = intc_pars[1], sd = intc_pars[2])# GP length scale
  data.frame(alpha, sigma, rho, intc)
}

log_gamma_statistic <- function(ranks, max_rank, ranks_to_check = NULL) {
  rank_counts <- rep(0, max_rank + 1)
  for(i in 0:max_rank) {
    rank_counts[i + 1] <- sum(ranks == i)
  }
  stopifnot(sum(rank_counts) == length(ranks))
  log_gamma_stat_counts(rank_counts, ranks_to_check)
}
log_gamma_stat_counts <- function(rank_counts, ranks_to_check = NULL) {
  max_rank <- length(rank_counts) - 1
  if(is.null(ranks_to_check)) {
    rank_ids_to_check = 1:max_rank
  } else {
    stopifnot(all(ranks_to_check >= 0 & ranks_to_check <= max_rank))
    rank_ids_to_check = setdiff(unique(ranks_to_check), max_rank) + 1
  }
  # Note that the ECDF follows the notation in Säilynoja, Bürkner & Vehtari 2022
  # i.e. it is the ECDF of a uniform distribution over j/(max_rank + 1) where j = 1:max_rank + 1
  scaled_ecdf <- cumsum(rank_counts)[rank_ids_to_check]
  z = rank_ids_to_check / (max_rank + 1)
  N_trials <- sum(rank_counts)
  return(log(2) + min(
    pbinom(scaled_ecdf, N_trials, z, log = TRUE),
    pbinom(scaled_ecdf - 1, N_trials, z, lower.tail = FALSE, log = TRUE)
  )
  )
}

gp_plot <- function (dim_id, output, gpfns, input, msd, esd, plot_type, label) {
  plot_data <- data.frame(output[, dim_id], gpfns[, dim_id], input, 
                          rep(msd[dim_id], N), rep(esd[dim_id], N))
  colnames(plot_data) <- c('output', 'gpfns', 'input', 'msd', 'esd')
  if (plot_type == 'gpfns') {
    plot <- ggplot(data = plot_data, aes(x = input, y = output)) +
      theme_bw(base_size = 20, base_family = 'Times') +
      geom_point() +
      geom_line(aes(y = gpfns), colour = c("#0072B2")) +
      geom_line(aes(y = gpfns - msd), colour = 'red', linetype = 'dashed') +
      geom_line(aes(y = gpfns + msd), colour = 'red', linetype = 'dashed') +
      labs(x = 'input', y = 'gpfns', title = label) 
  } else {
    plot <- ggplot(data = plot_data, aes(x = input, y = output)) +
      theme_bw(base_size = 20, base_family = 'Times') +
      geom_point() +
      geom_line(aes(y = output), colour = c("#0072B2")) +
      geom_line(aes(y = output - esd), colour = 'red', linetype = 'dashed') +
      geom_line(aes(y = output + esd), colour = 'red', linetype = 'dashed') +
      labs(x = 'input', y = 'outputs', title = label)
  }
  return(plot)
}
