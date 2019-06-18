#' Utility function to estimate detection probability given joint likelihood parameter estimates
#'
#' @param key key function (HN or HR)
#' @param params estimated parameters from the combined likelihood with either HN or HR key
#' @param survey.truncation truncation distance for detection function
#'
#' @return Phat estimated probability of detection, integrated over detection distances out to w
#'
estimate.Phat <- function(key, params, survey.truncation) {
  if (key == "HN") {
    Phat <- integrate(fx.denom,0,survey.truncation,
                      sigma1=params[1],sigma2=params[2],beta=params[3],
                      w=survey.truncation)$value
  } else {
    Phat <- integrate(fx.denomHR,0,survey.truncation,
                      sigma1=params[1],b=params[2],sigma2=params[3],beta=params[4],
                      w=survey.truncation)$value
  }
  return(Phat)
}
