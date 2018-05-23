#' Laplace approximation
#' 
#' Function to compute a Laplace type approximation for marginal likelihood of 
#' the mixed model where the random effects are inegrated out.
#' 
#' @param formula formula; formula of the mixed effects BERGM. 
#' @param theta numeric vector; value(s) for the fixed effects parameters.
#' @param phi numeric vector; values for the random effects parameters.
#' @param sigma_phi numeric; value for the variance of the random effects phi.
#' @param mu_phi numeric; value for the mean of the random effectsd phi.
#' @param t_obs count vector; values for the observed network statistics.
#' @param log.scale logical; returned value on log-scale? 
#' Default is \code{TRUE}.
#' @param iter count; number of iterations for the Laplace approximation, 
#' default is 10000.
#' @param aux.iter count; number of iterations for network simulation 
#' during the Laplace approximation, default is 3000.
#' 
#' @return The approximated value, by default on a log-scale (see argument 
#' \code{log.scale}).
#' 
#' @import ergm
#' @import network
#' 
f_laplace <- function(formula,
                      theta, phi, sigma_phi, mu_phi, 
                      t_obs,
                      log.scale = TRUE,
                      iter = 10000,
                      aux.iter = 3000) {
  
  net.start <- network(length(phi), directed = FALSE)
  net.formula <- as.formula(paste("~", 
                                  paste(attr(terms.formula(formula), 
                                             "term.labels"), collapse = " + "), 
                                  collapse = " "))
  sim_stats <- simulate(net.formula, basis = net.start, coef = c(theta, phi), 
                        nsim = iter, statsonly = TRUE, 
                        control = control.simulate(MCMC.burnin = 100000,
                                                   MCMC.interval = aux.iter))
  t_sim <- sim_stats[, -(1:length(theta))]
  
  Sig <- - 1 / sigma_phi * diag(length(phi)) - cov(t_sim)
  
  res <- (- length(phi) / 2) * log(sigma_phi) + t(phi) %*% t_obs - 
    ( t(phi - mu_phi) %*% (phi - mu_phi) ) / (2 * sigma_phi) - 
    0.5 * log(det(Sig))
  
  if(log.scale) {
    return(res)
  } else {
    return(exp(res))
  }
}

#' Path sampling
#' 
#' Function to estimate ratio of normalizing constants via path sampling.
#' @param mixed.formula formula; formula for the mixed model.
#' @param mixed.theta numeric vector; value(s) for the fixed effects parameters 
#' in the mixed model.
#' @param mixed.phi numeric vector; value(s) for the random effects parameters 
#' in the mixed model.
#' @param fixed.formula formula; formula for the fixed model.
#' @param fixed.theta numeric vector; value(s) for the parameters in the fixed 
#' model.
#' @param grid.points count; number of grid points for the path sampling, 
#' default is 1000.
#' @param iter count; number of iterations on each grid point for the path 
#' sampling, default is 10000.
#' @param aux.iter count; number of iterations for network simulation 
#' during the path sampling, default is 3000.
#' @param use.multicore logical; if set to \code{TRUE} (and if package 
#' \code{\link{parallel}} is available) the path sampling is done in parallel, 
#' default is \code{FALSE}.
#' @param mc.cores count; number of cores to use for parallel path sampling if 
#' \code{use.multicore = FALSE}, default is 2.
#' @param log.scale logical; returned value on log-scale? 
#' Default is \code{TRUE}.
#' @param stats.log logical; if set to \code{TRUE} an additional argument is  
#' returned, default is \code{FALSE}.
#' 
#' @return An estimate for kappa(fixed model)/kappa(mixed model) for the 
#' Bayes factor computation, default is on log-scale. \cr
#' Optional argument is \code{path.stats} containing the number of edges of each
#' simulated network, only returned if \code{stats.log = TRUE}.
#' 
#' @import ergm
#' @import network 
#' 
ratio_estimation <- function(mixed.formula, mixed.theta, mixed.phi,
                             fixed.formula, fixed.theta,
                             grid.points = 1000,
                             iter = 10000,
                             aux.iter = 3000,
                             use.multicore = FALSE, mc.cores = 2,
                             log.scale = TRUE,
                             stats.log = FALSE) {
  ## Define grid
  g.grid <- seq(from = 0, to = 1, 
                by = 1/grid.points) # i.e. equal spacing of the points g_i
  
  coef.grid <- list()
  for (i in 1:length(g.grid)) {
    theta_tmp <- (1 - g.grid[i]) * fixed.theta + g.grid[i] * c(0, mixed.theta)
    phi_tmp <- g.grid[i] * mixed.phi
    coef.grid[[i]] <- c(theta_tmp, phi_tmp)
  }
  rm(theta_tmp, phi_tmp)
  
  theta.diff <- fixed.theta - c(0, mixed.theta)
  phi.diff <- - mixed.phi
  
  formula.sim <- as.formula(paste("~ edges + ", 
                                  paste(attr(terms.formula(mixed.formula), 
                                             "term.labels"), collapse = " + "),
                                  collapse = " "))
  
  FUN <- function(coef, formula.sim, nnodes, 
                  theta.diff, phi.diff, 
                  iter, aux.iter, stats.log = FALSE) {
    net.start <- network(nnodes, directed = FALSE)
    sim_stats <- simulate(formula.sim, basis = net.start, coef = coef, 
                          nsim = iter, statsonly = TRUE, 
                          control = control.simulate(MCMC.burnin = 100000,
                                                     MCMC.interval = aux.iter))
    calc_diff <- function(u, theta.diff, phi.diff) {
      t(c(theta.diff, phi.diff)) %*% u
    }
    exp.val <- 1/iter * sum(apply(sim_stats, 1, calc_diff, 
                                  theta.diff = theta.diff, phi.diff = phi.diff))
    if(stats.log) {
      attr(exp.val, "sim_stats") <- sim_stats[, 1] # log no. of edges
    }
    return(exp.val)
  }
  
  if(use.multicore) {
    exp_val <- mclapply(coef.grid,
                        FUN = FUN, 
                        formula.sim = formula.sim, nnodes = length(mixed.phi),
                        theta.diff = theta.diff, phi.diff = phi.diff, 
                        iter = iter, aux.iter = aux.iter,
                        mc.cores = mc.cores, stats.log = stats.log)
  } else {
    exp_val <- lapply(coef.grid, FUN = FUN,
                      formula.sim = formula.sim, nnodes = length(mixed.phi),
                      theta.diff = theta.diff, phi.diff = phi.diff, 
                      iter = iter, aux.iter = aux.iter, stats.log = stats.log)
  }
  
  ## Calculate approximated integral  
  res <- 0
  for (i in 1:(length(exp_val) - 1)) {
    res <- res + 
      (g.grid[i + 1] - g.grid[i]) * (exp_val[[i + 1]] + exp_val[[i]]) / 2
  }
  
  if(stats.log) {
    attr(res, 
         "path.stats") <- lapply(exp_val, 
                                 FUN = function(u) {
                                   return(attr(u, "sim_stats"))
                                 })
  }
  
  if(log.scale) {
    return(res)
  } else {
    return(exp(res))
  }
}