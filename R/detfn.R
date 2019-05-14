#' Distance sampling detection function
#'
#' Computes probability of detection at a given distance when provided with key function, adjustments and parameters
#'
#' @param z generic distance (radial or perpendicular) - can be scalar or vector
#' @param pars the detection function parameter values (HN:s,adj1,adj2,...)(HR:s,b,adj1,adj2,...)(UNI:adj1,adj2,...)
#' @param key the detection function key, default is Half Normal (HR and UNI might be alternatives)
#' @param adjn if 0, no adjustment terms are used, if > 0, number of adjustments to use
#' @param adjt type of adjustments to use, defaults to cosine series
#' @param w truncation distance, only fundamental for adjustments, by default the max of the distances
#'
#' @return probability of detection
#' @export
#'
#'
detfn <- function(z,
                  pars,
                  key = "HN",
                  adjn = 0,
                  adjt = "cos",
                  w = max(z)) {
  #the standard HN and HR
  if (adjn == 0 & key == "HN") {
    p <- exp(-z ^ 2 / (2 * pars[1] ^ 2))
  }
  if (adjn == 0 & key == "HR") {
    p <- 1 - exp(-(z / pars[1]) ^ -pars[2])
  }
  #including cosines in the HN key
  if (adjn > 0 & key == "HN") {
    hn.bit <- exp(-z ^ 2 / (2 * pars[1] ^ 2))
    cos.bits <- matrix(nrow = length(z), ncol = adjn) #storage for the cosine terms
    for (j in 1:adjn) {
      #cos.bits[,j]<-pars[j+1]*cos(2*j*pi*z/w) #como estava inicialmente INCORRECTO
      cos.bits[, j] <- pars[j + 1] * cos((j + 1) * pi * z / w)
    }
    p <- (hn.bit * (1 + apply(cos.bits, 1, sum))) /
         (1 + sum(pars[(length(pars) - adjn + 1):length(pars)]))
  }
  #including cosines in the Uniform key
  if (adjn > 0 & key == "UNI") {
    cos.bits <- matrix(nrow = length(z), ncol = adjn) #storage for the cosine terms
    for (j in 1:adjn) {
      cos.bits[, j] <- pars[j + 1] * cos(j * pi * z / w)
    }
    p <- (1 * (1 + apply(cos.bits, 1, sum))) /
         (1 + sum(pars[(length(pars) - adjn + 1):length(pars)]))
  }
  return(p)
}
