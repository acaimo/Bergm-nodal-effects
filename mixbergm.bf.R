#' Bayes Factor computation for mixbergm model
#' 
#' Function to compute a Bayes Factor to compare two nested models: a Bayesian 
#' exponential random graph model with nodal random effects to one without 
#' random effects.
#' 
#' @param modmixed object; generated using the \code{\link{mixbergm}} function
#' @param modfixed object; generated using the \code{\link{bergm}} function. 
#' @param laplace.iter count; number of iterations for the build-in Laplace 
#' approximation, default is 10000.
#' @param laplace.aux.iter count; number of iterations for network simulation 
#' during the Laplace approximation, default is 3000.
#' @param path.iter count; number of iterations on each grid point for the path 
#' sampling, default is 10000.
#' @param path.aux.iter count; number of iterations for network simulation 
#' during the path sampling, default is 3000.
#' @param grid.points count; number of grid points for the path sampling, 
#' default is 1000.
#' @param use.multicore logical; if set to \code{TRUE} (and if package 
#' \code{\link{parallel}} is available) the path sampling is done in parallel, 
#' default is \code{FALSE}.
#' @param mc.cores count; number of cores to use for parallel path sampling if 
#' \code{use.multicore = FALSE}, if needed and not specified otherwise two cores 
#' are used.
#' @param stats.log logical; if set to \code{TRUE} some additional arguments are 
#' returned, default is \code{FALSE}. See details.
#' 
#' @details The two models need to be nested, in the sense that \code{modmixed} 
#' consists of same model terms, i.e. \code{\link[ergm]{ergm-terms}}, as 
#' \code{modmixed} except for the \code{edges} term, but with nodal random 
#' effects in addition. \cr
#' The fixed model in \code{modfixed}, generated using the \code{\link{bergm}} 
#' function, has to consist of a single chain. Models with multiple chains can 
#' not be handled a the moment. \cr
#' If the argument \code{stats.log} is set to \code{TRUE}, the returned object 
#' contains two additional arguments:
#' \itemize{
#'  \item{\code{components}}{ numeric vector; containing the individual 
#'  components of the log Bayes Factor.}
#'  \item{\code{path.stats}}{ count vector; number of edges of each simulated 
#'  network for the path sampling.}
#' }
#' 
#' @return The computed Bayes factor (on a log-scale) is returned. \cr
#' A value of 0 (or close to 0) means both models are equally good. 
#' If the value is positive the model with random effects is better, 
#' if the value is negative the model without nodal random effects is 
#' preferable. 
#'
#' @references
#' Stephanie Thiemichen, Nial Friel, Alberto Caimo, Göran Kauermann (2014). 
#' Bayesian Exponential Random Graph Models with Nodal Random Effects. 
#' \url{http://arxiv.org/abs/1407.6895}
#' 
#' @seealso \code{\link{mixbergm}}, \code{\link{bergm}}
#' 
#' @examples 
#' \dontrun{
#' # Note: The following example does not contain sufficient iterations, so that
#' # the models do not converge. It is just a rather quick example for 
#' # illustration purposes.
#' 
#' # Load the florentine marriage network
#' data("florentine", package = "ergm")
#' 
#' # Seed is needed for reproducible results
#' set.seed(123)
#' 
#' # Fit standard single-chain BERGM
#' resfixed <- bergm(flomarriage ~ edges + kstar(2), nchains = 1)
#' 
#' # Fit a model with modal random effects in addition
#' resmixed <- mixbergm(flomarriage ~ kstar(2) + sociality(base = 0))
#' 
#' # Compute log Bayes Factor
#' resBF <- mixbergm.bf(modmixed = resmixed, modfixed = resfixed,
#'                      laplace.iter = 100, laplace.aux.iter = 3000,
#'                      path.iter = 100, path.aux.iter = 3000, 
#'                      grid.points = 100)
#' resBF
#' }
#' 
#' @import mvtnorm
#' 
#' @export mixbergm.bf
#' 

mixbergm.bf <- function(modmixed, modfixed,
                        laplace.iter = 10000,
                        laplace.aux.iter = 3000,
                        path.iter = 10000,
                        path.aux.iter = 3000,
                        grid.points = 1000,
                        use.multicore = FALSE,
                        mc.cores = NULL,
                        stats.log = FALSE) {
  
  if(modfixed$nchains > 1) {
    stop(paste("Only fixed models with a single chain can be handled at the 
               moment."))
  }
  
  if(use.multicore && is.null(mc.cores)) {mc.cores <- 2}
  
  ## Extract posterior means
  modmixed.theta <- apply(modmixed$Theta, 2, mean)
  modmixed.phi <- apply(modmixed$Phi, 2, mean)
  modmixed.mu_phi <- mean(modmixed$Mu_phi)
  modmixed.sigma_phi <- exp(mean(log(modmixed$Sigma_phi)))
  
  modfixed.theta <- apply(modfixed$Theta, 2, mean)
  
  ## Extract observed network statistics
  modmixed.soc_index <- (length(modmixed$stats) - 
                           modmixed$nnodes + 1):(length(modmixed$stats))
  modmixed.s_obs <- modmixed$stats[-modmixed.soc_index] 
  modmixed.t_obs <- modmixed$stats[modmixed.soc_index]
  rm(modmixed.soc_index)
  
  modfixed.s_obs <- modfixed$stats 
  
  ## Extract prior densities at posterior means/modes
  modmixed.prior <- dmvnorm(modmixed.theta, 
                            mean = modmixed$m.prior, 
                            sigma = modmixed$sigma.prior, 
                            log = TRUE) + 
    dgamma((1 / modmixed.sigma_phi), 
           shape = modmixed$sigma.ab[1], rate = modmixed$sigma.ab[2], 
           log = TRUE) +
    dnorm(modmixed.mu_phi, mean = 0, sd = modmixed$rho, log = TRUE)
  
  modfixed.prior <- dmvnorm(modfixed.theta, 
                            mean = modfixed$m.prior, 
                            sigma = modfixed$sigma.prior, 
                            log = TRUE)
  
  ## Laplace Approximation for model with random effects #######################
  modmixed.laplace <- f_laplace(formula = modmixed$formula,
                                theta = modmixed.theta, phi = modmixed.phi,
                                sigma_phi = modmixed.sigma_phi, 
                                mu_phi = modmixed.mu_phi,
                                t_obs = modmixed.t_obs, log.scale = TRUE,
                                iter = laplace.iter, 
                                aux.iter = laplace.aux.iter)

  ## Path Sampling for ratio of normalizing constants ##########################
  ratio.est <- ratio_estimation(mixed.formula = modmixed$formula, 
                                mixed.theta = modmixed.theta,
                                mixed.phi = modmixed.phi,
                                fixed.formula = modfixed$formula, 
                                fixed.theta = modfixed.theta,
                                grid.points = grid.points,
                                iter = path.iter,
                                aux.iter = path.aux.iter,
                                use.multicore = use.multicore,
                                mc.cores = mc.cores,
                                log.scale = TRUE,
                                stats.log = stats.log)
  
  ## Estimate posterior densities
  modmixed.post <- dmvnorm(c(modmixed.theta, modmixed.mu_phi, 
                             log(modmixed.sigma_phi)), 
                           mean = c(modmixed.theta, modmixed.mu_phi, 
                                    log(modmixed.sigma_phi)), 
                           sigma = cov(cbind(modmixed$Theta, 
                                             modmixed$Mu_phi, 
                                             log(modmixed$Sigma_phi))), 
                           log = TRUE) - 
    log(modmixed.sigma_phi)
  modfixed.post <- dmvnorm(modfixed.theta, mean = modfixed.theta, 
                           sigma = cov(modfixed$Theta), log = TRUE)
  
  ## Calculate BF on log-scale
  bf.log <- t(modmixed.theta) %*% modmixed.s_obs + modmixed.laplace - 
    t(modfixed.theta) %*% modfixed.s_obs + 
    ratio.est + 
    modfixed.post - modmixed.post +
    modmixed.prior - modfixed.prior
  
  ## Optional output and return arguments
  if(stats.log) {
    attr(bf.log, "components") <- list(thetaS.modmixed = t(modmixed.theta) %*% 
                                         modmixed.s_obs,
                                       modmixed.laplace = modmixed.laplace,
                                       thetaS.modfixed = t(modfixed.theta) %*% 
                                         modfixed.s_obs,
                                       ratio.est = ratio.est,
                                       modmixed.post = modmixed.post,
                                       modfixed.post = modfixed.post,
                                       modmixed.prior = modmixed.prior,
                                       modfixed.prior = modfixed.prior)
    attr(bf.log, "path.stats") <- attr(ratio.est, "path.stats")
  }
  
  return(bf.log)
}
