#' Component to modelling the density gradient
#'
#' #' Not to be called by user; support function for \code{likdensity}
#'
#' @param x the perpendicular (or radial) distance
#' @param sigma a gaussian standard deviation
#' @param beta a multiplicative constant, helping to improve the fit notice that if beta>0, animals avoid linear structure, if beta<0, they prefer it
#'
#' @return component of likelihood
#'
pi.x.den <- function(x, sigma, beta) {
  denom <- (1 - beta * dnorm(x, mean=0, sd=sigma))
  return(denom)
}

#' One form of modelling the density gradient
#'
#' Not to be called by user; support function for \code{likdensity}
#'
#' @param x the perpendicular distance
#' @param sigma a gaussian standard deviation
#' @param beta a multiplicative constant, helping to improve the fit notice that if beta>0, animals avoid linear structure, if beta<0, they prefer it
#' @param w truncation distance
#'
#' @return component of likelihood
#'
pi.x <- function(x, sigma, beta, w) {
  pi.x.denK <- integrate(pi.x.den, lower = 0, w, sigma = sigma, beta = beta)
  pi.x.ret <- (1 - beta * dnorm(x, mean=0, sd=sigma)) / pi.x.denK$value
  return(pi.x.ret)
}

# #just a test
# pi.x(2, 20, 1.2, 100)
# #just a few plots
# xs = 0:150
# par(mfrow = c(2, 2), mar = c(4, 4, 0.5, 0.5))
# plot(xs, pi.x(xs, 20, 1.2, 100), type = "l")
# plot(xs, pi.x(xs, 20, -1.2, 100), type = "l")
# plot(xs, pi.x(xs, 50, 1.2, 100), type = "l")
# plot(xs, pi.x(xs, 50, -1.2, 100), type = "l")


#'  Likelihood function to estimate parameters of the density gradient
#'
#' @param par parameters of density gradient
#' @param distances distances of animal distribution from transect
#' @param truncation truncation distance
#'
#' @return evaluated likelihood of density gradient function
#' @export
#'
likdensity <- function(par, distances, truncation) {
  sigma <- par[1]
  beta <- par[2]
  likhood <- sum(log(pi.x(x = distances, sigma = sigma, beta = beta, w = truncation)))
  return(likhood)
}
