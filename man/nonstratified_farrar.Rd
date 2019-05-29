\name{nonstratified_farrar}
\alias{nonstratified_farrar}
\docType{data}
\title{
Data from line transect survey of kangaroos
}
\description{
Traditional line transect distance sampling data set.  Required fields are \code{Transect}, \code{Length}, \code{Distance} and \code{group}
}
\usage{data("nonstratified_farrar")}
\format{
  A data frame with 373 observations on the following 6 variables.
  \describe{
    \item{\code{Stratum}}{a factor with levels \code{Farrar}}
    \item{\code{Area}}{size of stratum}
    \item{\code{Transect}}{a factor with levels \code{transect0} \code{transect1} \code{transect10} \code{transect12} \code{transect14} \code{transect15} \code{transect16} \code{transect18} \code{transect19} \code{transect20} \code{transect21} \code{transect23} \code{transect4} \code{transect7} \code{transect8}}
    \item{\code{Length}}{transect length}
    \item{\code{Distance}}{perpendicular distance from transect}
    \item{\code{group}}{number of individuals in detected group}
  }
}
\details{
Placeholder
}
\source{Marques, T. A., Buckland, S. T., Bispo, R., & Howland, B. (2013). Accounting for animal density gradients using independent information in distance sampling surveys. Statistical Methods & Applications, 22(1), 67-80.
\url{https://doi.org/10.1007/s10260-012-0223-2}
}
\references{
Placeholder
}
\examples{
data(nonstratified_farrar)
## maybe str(nonstratified_farrar) ; plot(nonstratified_farrar) ...
}
\keyword{datasets}
