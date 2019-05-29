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
#' @param conversion conversion factor to get density units reported
#'
#' @return upper and lower confidence interval bounds on individual density estimates
#'
#' @section To Do
#'    incorporate alpha into computations (done)
#'    might it be wise to send back the entire set of simulation results (par.estimates)?
#'
varbootstrap <- function(keyfn, survey, gpsdata, disttrunc=150, gradientmaxdist=100,
                         Nboot, alpha=0.05, plotit=TRUE, conversion=10000) {
#  library(MIfuns)
  set.seed(123)
  surv.resamp <- MIfuns::resample.data.frame(survey,names=1:Nboot,key="Transect",rekey=TRUE)
  real.resamp <- MIfuns::resample.data.frame(gpsdata,names=1:Nboot,key="animal",rekey=TRUE)
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
      par.estimates$P[i] <- integrate(fx.denomHR,0,disttrunc,sigma1=mle2.b$par[1],
                                   b=mle2.b$par[2],sigma2=mle2.b$par[3],beta=mle2.b$par[4],w=disttrunc)$value
    }
    par.estimates$D[i] <- conversion * par.estimates$n[i]*par.estimates$EsSb[i]/
                                  (2*par.estimates$L[i]*disttrunc*par.estimates$P[i])
  }
  alphaover2 <- alpha/2
  CIlimits <- c(alphaover2, 1-alphaover2)
  CIbounds <- quantile(par.estimates$D, probs = CIlimits)
  if(plotit) {
    hist(par.estimates$D, main="Bootstrap density estimates",
         xlab="Estimated density", sub=paste("Successes",sum(!is.na(par.estimates$D))))
    abline(v=CIbounds, lwd=2, lty=2)
  }
  return(CIbounds)
}
