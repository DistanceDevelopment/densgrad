#' Support function for full.lik.fx.HR
#'
#' Not to be called by user; support function for \code{full.lik.fx.HR}
#'
#' @param x detection distances
#' @param sigma1 shape parameter of hazard rate detection fn
#' @param b scale parameter of hazard rate detection fn
#' @param sigma2 parameter of density gradient
#' @param beta attraction/repulsion parameter of density gradient
#' @param w truncation distance (distance sampling)
#'
#' @return component of likelihood
#'
fx.denomHR <- function(x,sigma1,b,sigma2,beta,w){
  denom <- detfn(x,pars=c(sigma1,b),key="HR",adjn=0,w=w)*
           pi.x(x,sigma=sigma2,beta=beta,w=w)
  return(denom)
}

#' Support function for full.lik.fx.HR
#'
#' Not to be called by user; support function for \code{full.lik.fx.HR}
#'
#' @param x detection distances
#' @param sigma1 shape parameter of hazard rate detection fn
#' @param b scale parameter of hazard rate detection fn
#' @param sigma2 parameter of density gradient
#' @param beta attraction/repulsion parameter of density gradient
#' @param w truncation distance (distance sampling)
#'
#' @return component of likelihood
#' @export
#'
fxHR <- function(x,sigma1,b,sigma2,beta,w){
  fx.denK <- integrate(fx.denomHR,0,upper=w,sigma1=sigma1,b=b,sigma2=sigma2,beta=beta,w=w)
  fx.ret <- fx.denomHR(x,sigma1=sigma1,b=b,sigma2=sigma2,beta=beta,w=w)/fx.denK$value
  return(fx.ret)
}

# fxHR(x=20,sigma1=50,b=1.6,sigma2=35,beta=35,w=150)

#' Combined likelihood with hazard rate detection and half normal gradient
#'
#' @param par four parameters of combined likelihood (2 for hazard rate, 2 for density gradient)
#' @param detectdists detection distances
#' @param gpsdists distances from feature from collared animals
#' @param truncdet truncation distance of distance sampling
#' @param truncden truncation distance of collaring study
#'
#' @return combined likelihood evaluated with HR detection function
#' @export
#'
full.lik.fx.HR <- function(par, detectdists, gpsdists, truncdet, truncden){
  sigma1 <- par[1]
  b <- par[2]
  sigma2 <- par[3]
  beta <- par[4]
  lik2 <- sum(log(fxHR(detectdists,sigma1=sigma1,b=b,sigma2=sigma2,beta=beta,truncdet)))+
          sum(log(pi.x(gpsdists,sigma=sigma2,beta=beta,w=truncden)))
  return(lik2)
}
