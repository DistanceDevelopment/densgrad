#' Variance estimation via bootstrap over transects
#'
#' @param keyfn key function of detection function (HN or HR)
#' @param survey distance sampling survey data frame
#' @param gpsdata animal location data frame from telemetry
#' @param disttrunc truncation distance for distance sampling component
#' @param gradientmaxdist maximum distance from transect of telemetry data
#' @param Nboot number of bootstrap replicates
#' @param alpha Type I error rate for computing CI bounds
#' @param plotit logical indicating whether to produce histogram of bootstrap estimates
#'
#' @return upper and lower confidence interval bounds on individual density estimates
#'
#' @section To Do
#'    incorporate alpha into computations
#'    might it be wise to send back the entire set of simulation results (par.estimates)?
#'
varbootstrap <- function(keyfn, survey, gpsdata, disttrunc=150, gradientmaxdist=100, Nboot, alpha=0.05, plotit=TRUE) {
  library(MIfuns)
  set.seed(123)
  surv.resamp <- resample.data.frame(survey,names=1:Nboot,key="Transect",rekey=TRUE)
  real.resamp <- resample.data.frame(gpsdata,names=1:Nboot,key="animal",rekey=TRUE)
  par.estimates <- data.frame(ID=1:Nboot,L=NA,n=NA,Es=NA,sigma1=NA,b=NA,sigma2=NA,beta=NA,P=NA,D=NA)
  for (i in 1:Nboot) {
    surv.data.b <- surv.resamp[[i]]
    real.resamp.b <- real.resamp[[i]]
    surv.data.b.dist <- surv.data.b$Distance[!is.na(surv.data.b$Distance) & surv.data.b$Distance<disttrunc]
    real.data.b.dist <- real.resamp.b$dist
    par.estimates$L[i] <- sum(tapply(X=surv.data.b$Length, INDEX=surv.data.b$Transect, FUN=mean))
    par.estimates$n[i] <- length(surv.data.b.dist)
    par.estimates$Es[i] <- mean(surv.data.b$group[!is.na(surv.data.b$group)])
    #quite dodgy, assumes correction is constant...
    #but any way, this is not used in the end
#    par.estimates$EsSb[i]=par.estimates$Es[i]*size.bias.corr
    par.estimates$EsSb[i] <- size.bias(keyfn, surv.data.b.dist, surv.resamp[[i]], grad.data=real.data.b.dist,
                                       disttrunc, gradientmaxdist, plotit=FALSE)
    if (keyfn =="HN") {
      mle2.b <- optim(par=c(50,35,35),fn=full.lik.fx.HN,detectdists=surv.data.b.dist,
                   gpsdists=real.data.b.dist,truncdet=disttrunc,truncden=gradientmaxdist,control=list(fnscale=-1))#,method="L-BFGS-B",lower=c(20,-100),upper=c(200,100))
      par.estimates$sigma1[i] <- mle2.b$par[1]
      par.estimates$sigma2[i] <- mle2.b$par[2]
      par.estimates$beta[i] <- mle2.b$par[3]
    } else {
      mle2.b <- optim(par=c(50,35,35),fn=full.lik.fx.HR,detectdists=surv.data.b.dist,
                   gpsdists=real.data.b.dist,truncdet=disttrunc,runcden=gradientmaxdist,control=list(fnscale=-1))#,method="L-BFGS-B",lower=c(20,-100),upper=c(200,100))
      par.estimates$sigma1[i] <- mle2.b$par[1]
      par.estimates$b[i] <- mle2.b$par[2]
      par.estimates$sigma2[i] <- mle2.b$par[3]
      par.estimates$beta[i] <- mle2.b$par[4]
    }
    if (keyfn =="HN") {
      par.estimates$P[i] <- integrate(fx.denom,0,disttrunc,sigma1=mle2.b$par[1],
                                   sigma2=mle2.b$par[2],beta=mle2.b$par[3],w=disttrunc)$value
    } else {
      par.estimates$P[i] <- integrate(fx.denomHR,0,disttrunc,sigma1=mle3.b$par[1],
                                   b=mle3.b$par[2],sigma2=mle3.b$par[3],beta=mle3.b$par[4],w=disttrunc)$value
    }
    par.estimates$D[i] <- 10000 * par.estimates$n[i]*par.estimates$EsSb[i]/
                                  (2*par.estimates$L[i]*disttrunc*par.estimates$P[i])
    if (i/10==round(i/10)) print(i)
  }
#  par.estimates$D <- 10000 * par.estimates$n*par.estimates$ECsSb2/
#                             (2*par.estimates$L*disttrunc*par.estimates$P)
  CIbounds <- quantile(par.estimates$D, probs = c(0.025, 0.975))
  if(plotit) {
    hist(par.estimates$D, main="Distribution of bootstrap density estimates",
         xlab="Estimated density (numbers per hectare)")
    abline(v=CIbounds, lwd=3, lty=2)
  }
  return(CIbounds)
}