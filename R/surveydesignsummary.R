
#' Extract information about a distance sampling survey from survey data frame
#'
#' @param survey data frame containing 'usual' distance sampling data
#' @param survey.truncation truncation distance for detection function fitting
#'
#' @return named list containing number of detections, survey effort and average group size
#'
survey.design.summary <- function(survey, survey.truncation) {
  num.detects <- sum(!is.na(survey$Distance) & survey$Distance<survey.truncation)
  total.effort <- sum(aggregate(Length~Transect, survey, min)$Length)
  average.group.size <- mean(survey$group[survey$Distance<survey.truncation], na.rm=TRUE)
  return(list(n=num.detects, L=total.effort, Es=average.group.size))
}
