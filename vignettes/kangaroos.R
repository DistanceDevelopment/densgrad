## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(densgrad)

## ------------------------------------------------------------------------
real <- read.csv(file="../data/GPSpoints.txt", header=TRUE, sep="\t")

## ------------------------------------------------------------------------
names(real) <- c("infid","nearfid","dist","reserve","veg","animal")

## ------------------------------------------------------------------------
survey <- read.csv(file="../data/nonstratified_farrar.txt", header=TRUE, sep="\t")

## ------------------------------------------------------------------------
survey.truncation <- 150
gradient.max.distance <- 100
survey.data <- survey$Distance[!is.na(survey$Distance) & survey$Distance<survey.truncation]

## ------------------------------------------------------------------------
gradient <- optim(par=c(35,1.2),
                  fn=likdensity,
                  distances=real$dist,
                  truncation=gradient.max.distance,
                  control=list(fnscale=-1))#,method="L-BFGS-B",lower=c(20,-100),upper=c(200,100))
if (gradient$convergence == 0) {
  sigma.est <- gradient$par[1]
  beta.est <- gradient$par[2]
  likelihood <- gradient$value
} else {
  print("Convergence not achieved, no parameter estimates produced")
} 

## ---- fig.cap="Fitted animal density gradient", fig.height=3-------------
xs <- seq(from=0, to=gradient.max.distance)
par(mfrow=c(1,1),mar=c(4,4,0.5,0.5))
hist(real$dist,prob=T,main="",xlab="Distance (m)")
lines(xs,pi.x(xs,sigma.est,beta.est,gradient.max.distance),type="l",lty=2,lwd=2)

## ------------------------------------------------------------------------
hn.det <- optim(par=c(50,35,35),
                fn=full.lik.fx.HN,
                detectdists=survey.data,
                gpsdists=real$dist,
                truncdet=survey.truncation,
                truncden=gradient.max.distance,
                control=list(fnscale=-1))#,method="L-BFGS-B",lower=c(20,-100),upper=c(200,100))
if (hn.det$convergence == 0) {
  sigma.detect.est <- hn.det$par[1]
  sigma.gradient.est <- hn.det$par[2]
  beta.est <- hn.det$par[3]
  likelihood.hn <- hn.det$value
  aic.hn <- -2 * likelihood.hn + 2 * length(hn.det$par)
} else {
  print("Convergence not achieved, no parameter estimates produced")
} 

## ------------------------------------------------------------------------
hz.det <- optim(par=c(50,1.2,35,35),
                fn=full.lik.fx.HR,
                detectdists=survey.data,
                gpsdists=real$dist,
                truncdet=survey.truncation,
                truncden=gradient.max.distance,
                control=list(fnscale=-1))#,method="L-BFGS-B",lower=c(20,-100),upper=c(200,100))
if (hz.det$convergence == 0) {
  sigma.detect.est <- hz.det$par[1]
  beta.detect.est <- hz.det$par[2]
  sigma.gradient.est <- hz.det$par[3]
  beta.gradient.est <- hz.det$par[4]
  likelihood.hz <- hz.det$value
  aic.hz <- -2 * likelihood.hz + 2 * length(hz.det$par)
} else {
  print("Convergence not achieved, no parameter estimates produced")
} 

## ------------------------------------------------------------------------
hn.estimates <- t(c(aic.hn, hn.det$par[1], NA, hn.det$par[2:3]))
hz.estimates <- t(c(aic.hz, hz.det$par))
param.table <- as.data.frame(rbind(hn.estimates, hz.estimates))
rownames(param.table) <- c("Half normal", "Hazard rate")
colnames(param.table) <- c("AIC","$\\widehat{\\sigma_D}$","$\\widehat{\\beta_D}$",
                           "$\\widehat{\\tau_G}$","$\\widehat{\\beta_G}$")
knitr::kable(param.table, digits=2, caption="Parameter estimates for two joint models.")

## ---- fig.cap="Fitted half normal detection function accounting for non-uniform distribution.", fig.height=3----
hist(survey$Distance[survey$Distance<survey.truncation],
     main="",yaxt="n",ylab="g(y)",xlab="Distance (m)")
xs <- seq(from=0, to=survey.truncation)
lines(xs,60*detfn(z=xs,pars=hn.det$par[1],key="HN",adjn=0,w=survey.truncation),
      ylab="g(y)",lwd=2,type="l",lty=3)
axis(2,at=seq(0,60,length=5),labels=seq(0,1,length=5))

## ------------------------------------------------------------------------
Phat <- estimate.Phat("HN", hn.det$par, survey.truncation)

## ------------------------------------------------------------------------
survey.design <- survey.design.summary(survey, survey.truncation)
n <- survey.design$n
L <- survey.design$L
Es <- survey.design$Es

## ------------------------------------------------------------------------
Dhat <- (n * Es)/(2*survey.truncation*L*Phat)*10000

## ------------------------------------------------------------------------
Esadj <- size.bias("HN", survey.data, survey, real$dist, survey.truncation, gradient.max.distance)
Dadj <- (n * Esadj)/(2*survey.truncation*L*Phat)*10000

## ------------------------------------------------------------------------
out.table <- data.frame(n=n, Es=Es, Esadj=Esadj, w=survey.truncation, L=L, 
                        P=Phat, D=Dhat, Dadj=Dadj)
colnames(out.table) <- c("n","E(s)","$E(s)_{sb}$","w","L","$\\hat{P}$","$\\hat{D}$","$\\hat{D_{sb}}$")
knitr::kable(out.table, caption="Summary statistics and estimate of density from selected model", 
             digits = 3)

## ---- warning=FALSE, fig.height=3----------------------------------------
bounds <- varbootstrap("HN", survey, real, survey.truncation, gradient.max.distance, 
                       Nboot=999, plotit=TRUE, alpha=0.05)
knitr::kable(data.frame(LCI=signif(unname(bounds$CIbounds[1]),4), UCI=signif(unname(bounds$CIbounds[2]),4)))

## ------------------------------------------------------------------------
alpha <- 0.05
alphaover2 <- alpha/2
CIlimits <- c(alphaover2, 1-alphaover2)
par(mfrow=c(2,2))
hist(bounds$par.estimates$sigma1, main=expression(widehat(sigma[d])), xlab="Estimated value")
abline(v=quantile(bounds$par.estimates$sigma1, probs = CIlimits), lwd=2, lty=3)
hist(bounds$par.estimates$sigma2, main=expression(widehat(sigma[g])), xlab="Estimated value")
abline(v=quantile(bounds$par.estimates$sigma2, probs = CIlimits), lwd=2, lty=3)
hist(bounds$par.estimates$beta, main=expression(widehat(beta)), xlab="Estimated value")
abline(v=quantile(bounds$par.estimates$beta, probs = CIlimits), lwd=2, lty=3)
hist(bounds$par.estimates$P, main=expression(widehat(P)), xlab="Estimated value")
abline(v=quantile(bounds$par.estimates$P, probs = CIlimits), lwd=2, lty=3)

