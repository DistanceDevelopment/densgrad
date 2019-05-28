#' Compute size-bias adjusted group size estimates
#'
#' @param key key function (HN or HR)
#' @param survey.data perpendicular detection distances from distance sampling survey
#' @param survey complete data frame of distance sampling survey data
#' @param grad.data animal telemetry data
#' @param survey.truncation truncation distance for distance sampling survey
#' @param gradient.max.distance maximum distance from transects for telemetry data
#' @param plotit flag determining plot of regression scatter plot
#'
#' @return adjusted group size
#'
size.bias <- function(key, survey.data, survey, grad.data=real$dist, survey.truncation,
                      gradient.max.distance, plotit=FALSE) {
  if (key=="HN") {
    g.of.x <- optim(par=c(50,35,35),fn=full.lik.fx.HN,
                    detectdists=survey.data,gpsdists=grad.data,
                    truncdet=survey.truncation,truncden=gradient.max.distance,
                    control=list(fnscale=-1))
    gxs <- detfn(z=survey$Distance[survey$Distance<survey.truncation],
                 pars=g.of.x$par[1],key="HN",adjn=0,w=survey.truncation)
  } else {
    g.of.x <- optim(par=c(50,1.2,35,35),fn=full.lik.fx.HR,
                    detectdists=survey.data,gpsdists=grad.data,
                    truncdet=survey.truncation,truncden=gradient.max.distance,
                    control=list(fnscale=-1))
    gxs <- detfn(z=survey$Distance[survey$Distance<survey.truncation],
                 pars=g.of.x$par[1:2],key="HR",adjn=0,w=survey.truncation)
  }
  group.size <- survey$group[survey$Distance<survey.truncation]
  regr.sizeongx <- lm(group.size~gxs)
  if(plotit) {
    plot(gxs,group.size)
    abline(lm(group.size~gxs))
  }
  adj.Es <- predict(regr.sizeongx, newdata=list(gxs=1))
  return(adj.Es)
}
