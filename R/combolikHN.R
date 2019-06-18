#' Component of combined likelihood with HN detection
#'
#' Not to be called by user; support function for \code{full.lik.fx.HN}
#'
#' @param x detection distances
#' @param sigma1 parameter of half normal detection fn
#' @param sigma2 parameter of density gradient
#' @param beta attraction/repulsion parameter of density gradient
#' @param w truncation distance
#'
#' @return component of likelihood
#'
fx.denom <- function(x,sigma1,sigma2,beta,w){
  denom <- detfn(x,pars=sigma1,key="HN",adjn=0,w=w)*pi.x(x,sigma=sigma2,beta=beta,w=w)
  return(denom)
}

#' Support function for full.lik.fx.HN
#'
#' Not to be called by user; support function for \code{full.lik.fx.HN}
#'
#' @param x detection distances
#' @param sigma1 parameter of half normal detection fn
#' @param sigma2 parameter of density gradient
#' @param beta attraction/repulsion parameter of density gradient
#' @param w truncation distance
#'
#' @return component of likelihood
#' @export
#'
fx <- function(x,sigma1,sigma2,beta,w){
  fx.denK <- integrate(fx.denom,lower=0,upper=w,sigma1=sigma1,sigma2=sigma2,beta=beta,w=w)
  fx.ret <- fx.denom(x,sigma1=sigma1,sigma2=sigma2,beta=beta,w=w)/fx.denK$value
  return(fx.ret)
}

#' Combined likelihood with half normal detection and half normal gradient
#'
#' @param par parameters of combined likelihood (sigma of detection, sigma of gradient, beta gradient)
#' @param detectdists detection distances
#' @param gpsdists distances from feature from collared animals
#' @param truncdet truncation distance of distance sampling
#' @param truncden truncation distance of collaring study
#'
#' @return value of combined likelihood with HN detection function
#' @export
#'
full.lik.fx.HN <- function(par,detectdists,gpsdists,truncdet,truncden){
  sigmadet <- par[1]
  sigmaden <- par[2]
  beta <- par[3]
  lik2 <- sum(log(fx(detectdists,sigma1=sigmadet,sigma2=sigmaden,beta=beta,truncdet)))+
          sum(log(pi.x(gpsdists,sigma=sigmaden,beta=beta,w=truncden)))
  return(lik2)
}
