#' Support function for \code{CDS.HN.lik}
#'
#'  Not to be called by user; helper function to CDS.HN.lik
#'
#' @param x distances
#' @param sigma1 half normal parameter
#'
#' @return vector of detection probabilities
#'
CDS.HN.fx.denom = function(x, sigma1) {
  denom <- detfn(x, pars = sigma1, key = "HN", adjn = 0)
  return(denom)
}

#' Support function for \code{CDS.HN.lik}
#'
#'  Not to be called by user; helper function to CDS.HN.lik
#'
#' @param x distances
#' @param sigma1 half normal parameter
#' @param w trunction
#'
#' @return component of likelihood
#'
CDS.HN.fx <- function(x, sigma1, w) {
  CDS.HN.fx.denK <- integrate(CDS.HN.fx.denom, lower = 0, upper = w, sigma1 = sigma1)
  CDS.HN.fx.ret <- CDS.HN.fx.denom(x, sigma1 = sigma1) / CDS.HN.fx.denK$value
  return(CDS.HN.fx.ret)
}


#' Likelihood of half-normal CDS detection function
#'
#'  To be used in conjunction with \code{optim()} to find MLEs of sigma parameter of half normal detection function
#'
#' @param par parameter vector being optimised
#' @param distances detection distances
#' @param truncation truncation distance
#'
#' @return evaluated likelihood for HN CDS
#'
CDS.HN.lik <- function(par, distances, truncation) {
  sigma1 <- par[1]
  likHN <- sum(log(CDS.HN.fx(distances, sigma1 = sigma1, truncation)))
  return(likHN)
}
