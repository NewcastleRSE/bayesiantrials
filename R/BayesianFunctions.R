################################################################################
#' Function for performing Bayesian inference using rjags
#'
#' @param tot Total sample size
#' @param y Response
#' @param X Treatment labels
#' @param C Number of clusters
#' @param cl Cluster label for each individual
#' @param K Number of MCMC samples
#' @param eps Analysis Gamma prior parameters
#'
#' @return MCMC samples
#' @import rjags
#' @export
#'
#' @examples
#'
#'
#'
Bayes.inf = function(tot, y, X, C, cl, K, eps = 0.001){

  dat <- list("N" = tot, "Y" = y, "X" = X, "Nclust" = C, "Cl"= cl,eps = eps)

  # mod = jags.model(file = "files/mixed.R", data = dat, n.adapt = 1000, quiet = TRUE)

  mod = jags.model(file = system.file("jags","mixed.R", package = "bayesiantrials"), data = dat, n.adapt = 1000, quiet = TRUE)

  update(mod, n.iter = K,progress.bar="none")
  params <- c("delta")

  samps <- coda.samples(mod, params, n.iter = K, progress.bar = "none")

  delta.out = samps[[1]]

  return(delta.out)

}

################################################################################

#' Function for calculating the assurance
#'
#' @param ss Sample size
#' @param design Design priors
#' @param L Monte Carlo Samples
#' @param K MCMC Samples
#' @param C Number of clusters
#' @param a Dirichlet parameter
#' @param eps Analysis Gamma prior parameters
#' @param alpha Posterior probability
#'
#' @return Assurance
#' @import DirichletReg
#' @export
#'
#' @examples
getAssurance = function(ss, design, L, K, C, a = 100, eps = 0.001, alpha = 0.05){

  I = 0

  for (i in 1:L){

    p = rdirichlet(1,rep(a,C))
    n = t(rmultinom(1,ss*C,p))

    cl = unlist(mapply(rep,seq_along(n), n))

    tot = sum(n)

    tot1 = sum(n[1:(C/2)])

    X = c(rep(1,tot1),rep(0,(tot-tot1)))

    epsilon = rnorm(tot,0,design[i,5])
    ce = rnorm(C,0,design[i,4])

    c.reg = ce[cl]

    y = design[i,1] + X*design[i,2] + c.reg + epsilon

    delta.out = Bayes.inf(tot, y, X, C, cl, K, eps)

    I = I + ifelse(quantile(delta.out,alpha)>0,1,0)

  }

    A = I/L

  return(A)
}


################################################################################

#' Function for finding the required sample size
#'
#' @param ss Vector of sample sizes to check
#' @param design Design priors
#' @param L Monte Carlo samples
#' @param K MCMC Samples
#' @param C Cluster Size
#' @param target Target assurance
#' @param a Dirichlet parameter
#' @param eps Analysis Gamma prior parameters
#' @param alpha Posterior probability
#'
#' @return Minimum sample size for target assurance. Returns target assurance not met for any sample size in ss.
#' @export
#'
#' @examples
getSampleSize = function(ss, design, L, K, C, target = 0.8, a = 100, eps = 0.001, alpha = 0.05){

  for (i in 1:length(ss)){
    assure = getAssurance(ss[i], design, L, K, C, a, eps, alpha)
    if (assure >= target){
      ss = ss[i]
      break
    }
  }
  if(length(ss) > 1){
    ss = NA
  }
  return(ss)
}

################################################################################




################################################################################
