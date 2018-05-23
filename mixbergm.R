#' Bayesian exponential random graph models with nodal random effects
#' 
#' Function to fit Bayesian exponential random graphs models with nodal random 
#' effects using the exchange algorithm. 
#' 
#' @param formula formula; an \code{R} formula object, of the form 
#' <network> ~ <model terms> where <network> is a \code{\link[network]{network}} 
#' object and <model terms> are \code{\link[ergm]{ergm-terms}}. 
#' The term \code{sociality(base = 0)} is added automatically to the right hand 
#' side of the formula argument.
#' 
#' @param burn.in count; number of burn-in iterations at the beginning of an 
#' MCMC run, default is 100. 
#' @param main.iters count; number of iterations for the MCMC chain excluding 
#' burn-in iterations, default is 1000. 
#' @param aux.iters count; number of auxiliary iterations used for network 
#' simulation in the exchange algorithm, default is 1000.
#' 
#' @param m.prior numeric vector; mean of the multivariate Normal prior for the 
#' fixed effects terms. By default set to a vector of 0's.
#' @param sigma.prior variance / covariance matrix for the multivariate Normal 
#' prior for the fixed effects terms. By default set to a diagonal matrix with 
#' every diagonal entry equal to 100.
#' 
#' @param gamma numeric; variance of the Normal proposal distribution for the 
#' fixed effects in the model, the default value is 0.5.
#' @param sigma.epsilon variance / covariance matrix for the multivariate Normal 
#' proposal. By default set to a diagonal matrix with every diagonal entry equal 
#' to 0.0025. 
#' If only one fixed effect term is included in the model, \code{sigma.epsilon} 
#' is set equal to \code{gamma}.
#' @param phi.gamma numeric; variance parameter for the Normal proposal density 
#' for the nodal random effects parameters, the default value is 0.5.
#' @param mu_phi.gamma numeric; variance parameter for the Normal proposal 
#' density for the mean parameter of the nodal random effects, the default value 
#' is 0.5.
#' @param sigma_phi.gamma numeric; defines the width of the proposal density for 
#' the variance parameter of the random effects, see details. The dafult value
#' is 0.5.
#' 
#' @param rho numeric; variance parameter of the prior Normal distribution of 
#' the mean parameter \code{mu_phi} of the nodal random effects, the default 
#' value is 100.
#' @param sigma.ab numeric vector; defining the shape and the rate parameter for 
#' the Gamma prior distribution of the variance parameter, the dafault is 
#' \code{c(shape = 0.001, rate = 0.001)}.  
#' 
#' @details  In its current form the implemented algorithm can only handle 
#' undirected networks. \cr
#' The fixed effects parameters are updated using a block-update. For the random 
#' effects a single-site update is used to achieve better acceptance rates.
#' 
#' @seealso \code{\link{mixbergm.bf}}, \code{\link{mixbergm.output}}
#' 
#' @references
#' Stephanie Thiemichen, Nial Friel, Alberto Caimo, Göran Kauermann (2014). 
#' Bayesian Exponential Random Graph Models with Nodal Random Effects. 
#' \url{http://arxiv.org/abs/1407.6895}
#' 
#' @examples 
#' \dontrun{
#' # Load the florentine marriage network
#' data("florentine", package = "ergm")
#' 
#' # Seed is needed for reproducible results
#' set.seed(123)
#' 
#' # Fit a model to measure the propensity to form 2-stars and capture nodal 
#' # heterogeneity
#' res <- mixbergm(flomarriage ~ kstar(2) + sociality(base = 0))
#' 
#' # Summarize results
#' mixbergm.output(res)
#' }
#' 
#' @import network
#' @import ergm
#' @import mvtnorm
#' 
#' @export
#' 

mixbergm <- function (formula, 

                      burn.in = 100,
                      main.iters = 1000,
                      aux.iters = 1000,
                      
                      m.prior = NULL, 
                      sigma.prior = NULL, 
                      gamma = 0.5, 
                      sigma.epsilon = NULL,
                      phi.gamma = 0.5,
                      mu_phi.gamma = 0.5,
                      sigma_phi.gamma = 0.5,
                      rho = 100,
                      sigma.ab = c(shape = 0.001, rate = 0.001)) {
  
  ## 1) Preparation and parameter initialisation ###############################
  y <- ergm.getnetwork(formula)
  if(is.directed(y)) {
    stop(paste("Only undirected networks can be handled at the moment."))
  }
  nnodes <- get.network.attribute(y, "n") ## number of nodes in the network
  
  mixed.only <- FALSE
  ## rearange formula
  formula <- formula_mixed(formula, form.env = environment(formula))
  if (length(attr(terms.formula(formula), "term.labels")) == 1) {
    mixed.only <- TRUE
  } 
  
  model <- ergm.getmodel(formula, y)
  Clist <- ergm.Cprepare(y, model)  
  stats0 <- summary(formula)
  
  control <- control.simulate.formula(MCMC.burnin = aux.iters, 
                                      MCMC.interval = 0)
  control$MCMC.samplesize <- 1
  
  MHproposal <- MHproposal.ergm(object = model, constraints = ~., 
                                arguments = control$MCMC.prop.args, 
                                nw = y, weights = control$MCMC.prop.weights, 
                                class = "c", reference = ~Bernoulli, 
                                response = NULL)     
  
  nstats <- Clist$nstats - nnodes
  
  ## results and current values
  Phi <- matrix(NA, main.iters, nnodes)
  phi <- runif(nnodes, min = -.1, max = .1)
  Sigma_phi <- vector(mode = "numeric", length = main.iters)
  sigma_phi <- 1   
  Mu_phi <- vector(mode = "numeric", length = main.iters)
  mu_phi <- runif(1, min = -.1, max = .1)
  
  ## acceptance counts
  phi.acc.counts <- rep(0, nnodes)
  mu_phi.acc.counts <- 0
  sigma_phi.acc.counts <- 0 

  if (is.null(m.prior))  m.prior <- rep(0, nstats) 
  if (is.null(sigma.prior)) sigma.prior <- diag(100, nstats) 
  if (is.null(sigma.epsilon)) sigma.epsilon <- diag(0.0025, nstats)
  if (nstats == 1) { 
    sigma.epsilon <- diag(gamma, nstats)  
  }
  
  Theta <- matrix(NA, main.iters, nstats)
  theta <- runif(nstats, min = -.1, max = .1)
  theta1 <- rep(0, nstats)
  
  if (mixed.only) {
    Theta <- NULL
    theta <- NULL
    theta1 <- NULL
  }
  
  tot.iters <- burn.in + main.iters
  
  acc.counts <- 0 
  
  
  ## 2) Actual Iterations ######################################################
  for (k in 1:tot.iters) {
    
    if (!(mixed.only)) {
      # (a) Update of theta ....................................................
      
      theta1 <- theta + rmvnorm(1, sigma = sigma.epsilon)[1,]
      pr <- dmvnorm(rbind(theta1, theta), mean = m.prior, 
                    sigma = sigma.prior)
      
      delta <- ergm.mcmcslave(Clist, MHproposal, eta0 = c(theta1, phi), 
                              control, verbose = FALSE)$s 
      
      beta <- t(c(theta, phi) - c(theta1, phi)) %*% t(delta) + log(pr[1] / pr[2])
      
      if (beta >= log(runif(1))) {
        theta <- theta1
        if (k > burn.in) acc.counts <- acc.counts + 1
      }     
      if (k > burn.in) Theta[(k - burn.in), ] <- theta
    }
    
    
    # (b) Update of phi ......................................................
    res <- update_phi(phi, mu_phi, sigma_phi, theta,
                      phi.gamma,
                      Clist, MHproposal, control,
                      phi.acc.counts, after.burn.in = (k > burn.in), 
                      mixed.only)
    phi <- res$phi
    phi.acc.counts <- res$phi.acc.counts
    
    if (k > burn.in) Phi[(k - burn.in), ] <- phi
    # ........................................................................         
    
    # (c) Update of mu_phi ...................................................
    res <- update_mu_phi(mu_phi, mu_phi.gamma,
                         phi, sigma_phi, nnodes, rho,
                         mu_phi.acc.counts,
                         after.burn.in = (k > burn.in))
    mu_phi <- res$mu_phi
    mu_phi.acc.counts <- res$mu_phi.acc.counts
    
    if (k > burn.in) Mu_phi[(k - burn.in)] <- mu_phi
    # ........................................................................            
    
    # (d) Update of sigma_phi ................................................
    res <- update_sigma_phi(sigma_phi, sigma_phi.gamma,
                            phi, mu_phi, nnodes,
                            sigma.ab, sigma_phi.acc.counts,
                            after.burn.in = (k > burn.in))
    sigma_phi <- res$sigma_phi
    sigma_phi.acc.counts <- res$sigma_phi.acc.counts
    
    if (k > burn.in) Sigma_phi[(k - burn.in)] <- sigma_phi
    # ........................................................................
  }
  
  ## 3) Return arguments ######################################################
  out = list(Clist = Clist,
             MHproposal = MHproposal,
             control = control,
             formula = formula,
             model = model,
             nnodes = nnodes,
             specs = model$coef.names,
             dim = Clist$nstats,
             stats = stats0,
             Theta = Theta,
             Phi = Phi,
             Sigma_phi = Sigma_phi,
             Mu_phi = Mu_phi,
             
             acc.rates = acc.counts / main.iters,
             phi.acc.rates = phi.acc.counts / main.iters,
             mu_phi.acc.rates = mu_phi.acc.counts / main.iters,
             sigma_phi.acc.rates = sigma_phi.acc.counts / main.iters,
             
             m.prior = m.prior,
             sigma.prior = sigma.prior,
             rho = rho, 
             sigma.ab = sigma.ab,
             aux.iters = aux.iters,
             time = time) 
  return(out)
}
