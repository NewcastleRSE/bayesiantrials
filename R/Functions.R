library(mvtnorm)

################################################################################
#' Calculate mode of vector
#'
#' @param v A numeric vector
#'
#' @return A number
#' @export
#'
#' @examples
#' getMode(c(1,2,3,2,3,4,2))
getMode = function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


################################################################################
#' Function to simulate from a bivariate Gaussian copula
#'
#' @param N Number of samples
#' @param rho Samples from design prior distribution for rho
#' @param m_sigma Mean of sigma
#' @param s_sigma Stdev of sigma
#' @param gamma Correlation parameter
#'
#' @return Joint prior distribution between rho and sigma
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats pnorm qnorm quantile rmultinom rnorm sigma update
#' @export
#'
#' @examples
simulateGaussianCopula <- function(N, rho, m_sigma, s_sigma, gamma){

  qmarg1 = function(p) quantile(rho,p)
  qmarg2 = function(p) qnorm(p,m_sigma,s_sigma)

  R <- rbind(c(1,gamma),c(gamma,1))

  dat <- rmvnorm(N, mean = c(0,0), sigma = R)

  dat[,1] <- qmarg1(pnorm(dat[,1]))
  dat[,2] <- qmarg2(pnorm(dat[,2]))

  return(dat)
}
