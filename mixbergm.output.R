#' MCMC output diagnostics
#' 
#' Function to generate simple MCMC output diagnsotics for a Bayesian 
#' Exponential Random Graph Model with nodal random effects.
#'  
#' @param x object; generated using \code{\link{mixbergm}}.
#' @param ... additional arguments, to be passed to lower level functions.
#' 
#' @details This function returns the posterior parameter density estimates and 
#' creates simple diagnostic plots for the MCMC produced from a 
#' \code{\link{mixbergm}} fit.
#' 
#' @seealso mixbergm
#' 
#' @examples
#' \dontrun{
#' # Load the florentine marriage network
#' data("florentine", package = "ergm")
#' 
#' set.seed(123)
#' res <- mixbergm(flomarriage ~ kstar(2) + sociality(base = 0))
#' 
#' # Summarize results
#' mixbergm.output(res)
#' }
#' 
#' @export mixbergm.output
#' 
#' @import coda
#' 
mixbergm.output <- function(x,
                            ...){
  
    cat("\n", "MCMC results for Model: y ~", paste(x$formula[3]), "\n")
    
    cat("\n", "Nodal random effects posterior density estimate:", "\n")
    
    RFtab <- rbind(c(mean(x$Mu_phi), mean(x$Sigma_phi)), 
                   c(sd(x$Mu_phi), sd(x$Sigma_phi)),
                   c(x$mu_phi.acc.rates, x$sigma_phi.acc.rates))
    rownames(RFtab) <- c("Post. mean:", "Post. sd:", "Acceptance rate:")
    colnames(RFtab) <- c("Mean mu_phi", "Variance sigma_phi")

    all <- as.table(RFtab)
    print(RFtab)
    
    FF <- x$Theta
    FFdim <- (x$dim - x$nnodes)
    
    if (!(is.null(FF))) {
      cat("\n",' Fixed effects posterior density estimate:',"\n")
      overall <- rbind(apply(FF, 2, mean), apply(FF, 2, sd))
      rownames(overall) <- c("Post. mean:", "Post. sd:")
      colnames(overall) <- paste("theta", seq(1, FFdim), " (",
                                 x$specs[seq(1, FFdim)], ")", 
                                 sep = "")
      all <- as.table(overall)
      print(overall)
    
      rates <- matrix(x$acc.rate, 1, FFdim)	
      cat("\n","Acceptance rate for fixed effects:", 
          format(x$acc.rates, digits = 2), "\n", "\n", "\n")
  }		
  
  # Plot results
  
  dev.new()
  
  oma <- c(0, 0, 3, 0) # outer margins of plot
  mar <- c(4, 3, 1.5, 1) # margins of plot
  
  # First plot random effects estimates
  par(mfrow = c(2, 3), oma = oma, mar = mar)
  
  plot(density(x$Mu_phi), main = "", axes = FALSE, 
       xlab = bquote(paste(mu[phi]," (mean)")),
       ylab = "", lwd = 2)
  axis(1); axis(2)
  plot(x$Mu_phi, type = "l", xlab = "Iterations", ylab = "")
  autocorr.plot(mcmc(x$Mu_phi), auto.layout = FALSE, ...)
  
  plot(density(x$Sigma_phi), main = "", axes = FALSE, 
       xlab = bquote(paste(sigma[phi]," (variance)")),
       ylab = "", lwd = 2)
  axis(1); axis(2)
  plot(x$Sigma_phi, type = "l", xlab = "Iterations", ylab = "")
  autocorr.plot(mcmc(x$Sigma_phi), auto.layout = FALSE, ...)
  
  title(paste("MCMC output for Model: y ~", x$formula[3]), outer = TRUE)
  
  # Plot fixed effects estimates
  dev.new()
  
  seqq <- 4 # max. no. of parameter diagnostics plots to display in each 
            # graphics window
   
  K <- mcmc(data = FF)
  par(mfrow = c(min(4, FFdim), 3), oma = oma, mar = mar)
  
  for(i in 1:FFdim){
    if(i%in%c(5, 9, 13)) {
      dev.new()
      par(mfrow = c(min(4, FFdim - (i-1)), 3), oma = oma, mar = mar)
    }
    plot(density(FF[,i]), main = "", axes = FALSE, 
         xlab = bquote(paste(theta[.(i)], " (", .(x$specs[i]), ")")),
         ylab = "", lwd = 2)
    axis(1); axis(2)
    plot(FF[,i], type = "l", xlab = "Iterations", ylab = "")
    autocorr.plot(K[,i], auto.layout = FALSE, ...)
    
    if(FFdim > 4) 
      seqq <- seq(4, FFdim, 4)
    if(i%in%union(FFdim, seqq)) 
      title(paste("MCMC output for Model: y ~", x$formula[3]), outer = TRUE)
  }
}