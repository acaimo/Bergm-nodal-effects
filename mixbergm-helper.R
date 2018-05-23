#' Mixed Bergm formula
#' 
#' Rearranges formula argument for mixed effects model.
#' 
#' @param formula formula; an \code{R} formula object, of the form 
#' <network> ~ <model terms> where <network> is a \code{\link{network}} object 
#' and <model terms> are \code{\link{ergm-terms}}. 
#' @param form.env \code{R} environment; default is \code{environment(formula)}.
#' 
#' @details The function restructers the right hand side of the \code{formula} 
#' object, so that \code{sociality(base = 0)} is the last term and additional 
#' \code{sociality} terms are removed. \cr
#' A warning is produced if the \code{edges} term is present in the 
#' \code{formula}. This term should not be included for the mixed Bergm.
#' 
#' @return An object of class \code{"formula"}.
#' 
formula_mixed <- function(formula, 
                          form.env = environment(formula)) {
  ## extract term.labels included in formula
  terms.vec <- attr(terms.formula(formula), "term.labels")
  
  if (any(grepl("edges", terms.vec))) {
    warning("The edges term should not be included when using mixbergm().")
  }
  
  if (length(terms.vec) == 1) {
    ## just random effects
    if (any(grepl("sociality", terms.vec))) { # sociality in formula?
      formula.string <- paste(attr(terms(formula, keep.order = TRUE), 
                                   "variables")[[2]], "~", 
                              "sociality(base = 0)")
      
    } else {# add sociality(base = 0) at the end of the formula argument
      formula.string <- paste(attr(terms(formula, keep.order = TRUE), 
                                   "variables")[[2]], "~", 
                              paste(terms.vec, collapse = " + "), 
                              "+ sociality(base = 0)")
    }
  } else {
    index <- 1:length(terms.vec)
    
    ## sociality already somewhere in the formula? - change index accordingly
    if (any(grepl("sociality", terms.vec))) {
      sociality.where <- which(grepl("sociality", terms.vec))
      index <- 1:length(terms.vec)
      index <- index[-sociality.where] # all terms except for sociality
    } 
    
    formula.string <- paste(attr(terms(formula, keep.order = TRUE), 
                                 "variables")[[2]], "~", 
                            paste(terms.vec[index], collapse = " + "), 
                            "+ sociality(base = 0)")
  }
  
  return(as.formula(formula.string, env = form.env))
}

#' Update nodal random effects
#' 
#' Function for the update of the nodal random effects.
#' 
#' @param phi numeric vector; current values of the nodal random effects 
#' parameters.
#' @param mu_phi numeric; current value of the mean parameter for the nodal 
#' random effects.
#' @param sigma_phi numeric; current value of the variance parameter for the 
#' nodal random effects.
#' @param theta numeric (vector); current value(s) of the remaining fixed 
#' effects parameters in the model, if there are any.
#' @param phi.gamma numeric; variance parameter for the Normal proposal density.
#' @param Clist
#' @param MHproposal
#' @param control
#' @param phi.acc.counts count vector; current acceptance counts for the nodal 
#' random effects parameters.
#' @param after.burn.in logical; indicator if the current iteration a main 
#' iteration (\code{TRUE}) or part of the burn-in (\code{FALSE}, default value).
#' @param mixed.only logical; indicator if the model consists of nodal random 
#' effects only, default value is \code{FALSE}.
#' 
#' @details
#' Each element nodal random effect parameter, i.e. each element of \code{phi} 
#' is updated individually in random order using the exchange algorithm. A new 
#' value is proposed using a Normal distribution centered at the current value 
#' with variance \code{phi.gamma}. The acceptance probability is then computed 
#' using the exchange algorithm.
#' 
#' @return A named list with two elements:
#' \itemize{
#'  \item{\code{phi}}{ numeric vector; containing the values of the nodal random 
#'               effects parameters after the update}
#'  \item{\code{phi.acc.counts}}{ count vector; acceptance counts after the update.}
#' }
#' 
#' @import mvtnorm
#' @import ergm
#' 
update_phi <- function(phi,
                       mu_phi,
                       sigma_phi,
                       theta = NULL,
                       phi.gamma,
                       Clist, MHproposal, control,
                       phi.acc.counts,
                       after.burn.in = FALSE,
                       mixed.only = FALSE) {
  nnodes <- length(phi)
  for (i in sample(seq(1:nnodes), nnodes)) { # random ordering of updates
    phi1 <- phi ## set working vector to current values
    
    ## propose new value for phi_i
    phi1[i] <- rnorm(n = 1, mean = phi[i], sd = sqrt(phi.gamma))             
    
    pr <- dmvnorm(rbind(phi1, phi), mean = rep(mu_phi, nnodes), 
                  sigma = diag(sigma_phi, nnodes))
    
    if (mixed.only) {
      eta0 <- phi1
      eta <- phi
    } else {
      eta0 <- c(theta, phi1)
      eta <- c(theta, phi)
    }
    
    delta <- ergm.mcmcslave(Clist, MHproposal, eta0 = eta0, 
                            control, verbose = FALSE)$s 
    
    beta <- (eta - eta0) %*% t(delta) + log(pr[1] / pr[2])
    
    if (beta >= log(runif(1))) {
      phi <- phi1
      if (after.burn.in) phi.acc.counts[i] <- phi.acc.counts[i] + 1
    }
  }
  return(list(phi = phi,
              phi.acc.counts = phi.acc.counts))
}

#' Update the mean parameter for the nodal random effects
#' 
#' Function for the update of the mean parameter of the nodal random effects.
#' 
#' @param mu_phi numeric; current value of the mean parameter for the nodal 
#' random effects.
#' @param mu_phi.gamma numeric; variance parameter for the Normal proposal 
#' density.
#' @param phi numeric vector; current values of the nodal random effects 
#' parameters.
#' @param sigma_phi numeric; current value of the variance parameter for the 
#' nodal random effects.
#' @param nnodes count; number of nodes in the network, i.e. number of random 
#' effects terms. 
#' @param rho numeric; variance parameter of the prior Normal distribution of 
#' the parameter \code{mu_phi}.
#' @param mu_phi.acc.counts count; current acceptance count.
#' @param after.burn.in logical; indicator if the current iteration a main 
#' iteration (\code{TRUE}) or part of the burn-in (\code{FALSE}, default value).
#' 
#' @details
#' The update for the mean mean parameter of the nodal random effects 
#' \code{mu_phi} is done using a classical Metropolis-Hastings update. 
#' A new value is proposed using a Normal distribution centered at the current 
#' value of \code{mu_phi} with variance \code{mu_phi.gamma}. If the new proposed 
#' value is rejected, the current value \code{mu_phi} is the new value.
#' 
#' @return A named list with two elements:
#' \itemize{
#'  \item{\code{mu_phi}}{ numeric; value of the mean parameter for the nodal 
#'                        random effects after the update. }
#'  \item{\code{mu_phi.acc.counts}}{ count; acceptance count after the update. }
#' }
#' 
#' @import mvtnorm
#' 
update_mu_phi <- function(mu_phi,
                          mu_phi.gamma,
                          phi, 
                          sigma_phi, 
                          nnodes = length(phi), 
                          rho,
                          mu_phi.acc.counts,
                          after.burn.in = FALSE) {
  ## propose new value for mu_phi
  mu_phi1 <- rnorm(n = 1, mean = mu_phi, sd = sqrt(mu_phi.gamma))
  
  beta <- dmvnorm(phi, mean = rep(mu_phi1, nnodes), 
                  sigma = diag(sigma_phi, nnodes), log = TRUE) + 
    dnorm(mu_phi1, mean = 0, sd = sqrt(rho), log = TRUE) -
    dmvnorm(phi, mean = rep(mu_phi, nnodes), 
            sigma = diag(sigma_phi, nnodes), log = TRUE) - 
    dnorm(mu_phi, mean = 0, sd = sqrt(rho), log = TRUE)
  
  if (beta >= log(runif(1))) {
    mu_phi <- mu_phi1 # i.e. proposal is accepted
    if (after.burn.in) mu_phi.acc.counts <- mu_phi.acc.counts + 1
  }     
  
  return(list(mu_phi = mu_phi,
              mu_phi.acc.counts = mu_phi.acc.counts))
}

#' Update the variance parameter for the nodal random effects
#' 
#' Function for the update of the variance parameter of the nodal random 
#' effects.
#' 
#' @param sigma_phi numeric; current value of the variance parameter for the 
#' nodal random effects.
#' @param sigma_phi.gamma numeric; defines the width of the proposal density, 
#' see details.
#' @param phi numeric vector; current values of the nodal random effects 
#' parameters.
#' @param mu_phi numeric; current value of the mean parameter for the nodal 
#' random effects.
#' @param nnodes count; number of nodes in the network, i.e. number of random 
#' effects terms. 
#' @param sigma.ab numeric vector; defining the shape and the rate parameter for 
#' the Gamma prior distribution of the variance parameter.  
#' @param sigma_phi.acc.counts count; current acceptance count.
#' @param after.burn.in logical; indicator if the current iteration a main 
#' iteration (\code{TRUE}) or part of the burn-in (\code{FALSE}, default value).
#' 
#' @details
#' The update for the mean variance parameter of the nodal random effects 
#' \code{sigma_phi} is done using a classical Metropolis-Hastings update. 
#' A new value is proposed using a uniform distribution centered at the current 
#' value of \code{sigma_phi} ranging from \code{sigma_phi - sigma_phi.gamma} to 
#' \code{sigma_phi + sigma_phi.gamma}. The distribution is cut a 0, to avoid 
#' negative values.
#' If the new proposed value is rejected, the current value \code{sigma_phi} is 
#' the new value.
#' 
#' @return A named list with two elements:
#' \itemize{
#'  \item{\code{sigma_phi}}{ numeric; value of the variance parameter for the 
#'                           nodal random effects after the update. }
#'  \item{\code{sigma_phi.acc.counts}}{ count; acceptance count after the 
#'                                      update. }
#' }
#' 
#' @import mvtnorm
#' 

update_sigma_phi <- function(sigma_phi,
                             sigma_phi.gamma,
                             phi,
                             mu_phi,
                             nnodes = length(phi),
                             sigma.ab,
                             sigma_phi.acc.counts,
                             after.burn.in = FALSE) {
  ## propose new value for sigma_phi
  if (sigma_phi > sigma_phi.gamma) {
    sigma_phi1 <- runif(n = 1, min = sigma_phi - sigma_phi.gamma, 
                        max = sigma_phi + sigma_phi.gamma) 
    h_0to1 <- 1 / (2 * sigma_phi.gamma)
    if (sigma_phi1 > sigma_phi.gamma) {
      h_1to0 <- 1 / (2 * sigma_phi.gamma)
    } else { ## sigma_phi1 < sigma_phi.gamma, i.e. sigma_phi1 too close to 0
      h_1to0 <- 1 / (sigma_phi.gamma + sigma_phi1)
    }
  } else { ## sigma_phi < sigma_phi.gamma, i.e. sigma_phi too close to 0
    sigma_phi1 <- runif(n = 1, min = 0, max = sigma_phi + sigma_phi.gamma)
    h_0to1 <- 1 / (sigma_phi.gamma + sigma_phi)
    if (sigma_phi1 > sigma_phi.gamma) {
      h_1to0 <- 1 / (2 * sigma_phi.gamma)
    } else {
      h_1to0 <- 1 / (sigma_phi.gamma + sigma_phi1)
    }
  }
  
  beta <- dmvnorm(phi, mean = rep(mu_phi, nnodes), 
                  sigma = diag(sigma_phi1, nnodes), log = TRUE) + 
    dgamma((1 / sigma_phi1), shape = sigma.ab[1], rate = sigma.ab[2], 
           log = TRUE) + log(h_1to0) -
    dmvnorm(phi, mean = rep(mu_phi, nnodes), 
            sigma = diag(sigma_phi, nnodes), log = TRUE) - 
    dgamma((1 / sigma_phi), shape = sigma.ab[1], rate = sigma.ab[2], 
           log = TRUE) - log(h_0to1)
  
  if (beta >= log(runif(1))) {
    sigma_phi <- sigma_phi1 # i.e. proposal is accepted
    if (after.burn.in) sigma_phi.acc.counts <- sigma_phi.acc.counts + 1
  }     
  return(list(sigma_phi = sigma_phi,
              sigma_phi.acc.counts = sigma_phi.acc.counts))
}