\name{GPSpoints}
\alias{GPSpoints}
\docType{data}
\title{
Locations of telemetered animals
}
\description{
Data derived from telemetry study; required fields are \code{dist} and \code{animal}
}
\usage{data("GPSpoints")}
\format{
  A data frame with 5393 observations on the following 6 variables.
  \describe{
    \item{\code{infid}}{a numeric vector}
    \item{\code{nearfid}}{a numeric vector}
    \item{\code{dist}}{a numeric vector}
    \item{\code{reserve}}{a factor with levels \code{Black Mt} \code{Crace Grassland} \code{Farrer Ridge} \code{Googong} \code{Gossan Hill} \code{Gungaderra} \code{Kaleen Horse Paddocks} \code{Lawson Territory} \code{Mt Ainslie} \code{Mt Majura} \code{Mt Taylor} \code{Mulanggari} \code{Oakey Hill} \code{Red Hill} \code{Wanniassa Hills}}
    \item{\code{veg}}{a factor with levels \code{} \code{Forest} \code{Grassland} \code{Open Forest} \code{Open Woodland} \code{Woodland}}
    \item{\code{animal}}{a numeric vector of animal (collar) ID}
  }
}
\details{
Placeholder
}
\source{
Marques, T. A., S. T. Buckland, R. Bispo, and B. Howland. 2013. Accounting for animal density gradients using independent information in distance sampling surveys. Statistical Methods and Applications 22:67–80. \url{https://doi.org/10.1007/s10260-012-0223-2}
}
\references{
Placeholder
}
\examples{
data(GPSpoints)
## maybe str(GPSpoints) ; plot(GPSpoints) ...
}
\keyword{datasets}
