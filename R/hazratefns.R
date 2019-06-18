#' Support function for \code{CDS.HR.lik}
#'
#'  Not to be called by user; helper function to CDS.HR.lik
#'
#' @param x detection distances
#' @param sigma1 shape parameter of hazard rate
#' @param b scale parameter of hazard rate
#'
#' @return component of likelihood
#'
CDS.HR.fx.denom <- function(x, sigma1, b) {
  denom <- detfn(x, pars = c(sigma1, b), key = "HR", adjn = 0)
}

#' Support function for \code{CDS.HR.lik}
#'
#'  Not to be called by user; helper function to CDS.HR.lik
#'
#' @param x distances
#' @param sigma1 shape parameter of hazard rate
#' @param b scale parameter of hazard rate
#' @param w truncation distance
#'
#' @return component of likelihood
#'
CDS.HR.fx <- function(x, sigma1, b, w) {
  CDS.HR.fx.denK <- integrate(CDS.HR.fx.denom, lower = 0, upper = w, sigma1 = sigma1, b = b)
  CDS.HR.fx.ret <- CDS.HR.fx.denom(x, sigma1 = sigma1, b = b) / CDS.HR.fx.denK$value
  return(CDS.HR.fx.ret)
}

#' Likelihood of hazard rate CDS detection function
#'
#' @param par parameters of the hazard rate detection function
#' @param distances detection distances
#' @param truncation truncation distance
#'
#' @return evaluated likelihood for HR CDS
#'
CDS.HR.lik <- function(par, distances, truncation) {
  sigma1 <- par[1]
  b <- par[2]
  lik2 <- sum(log(CDS.HR.fx(distances, sigma1 = sigma1, b = b, w = truncation )))
  return(lik2)
}
