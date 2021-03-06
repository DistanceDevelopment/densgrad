% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combolikHR.R
\name{full.lik.fx.HR}
\alias{full.lik.fx.HR}
\title{Combined likelihood with hazard rate detection and half normal gradient}
\usage{
full.lik.fx.HR(par, detectdists, gpsdists, truncdet, truncden)
}
\arguments{
\item{par}{four parameters of combined likelihood (2 for hazard rate, 2 for density gradient)}

\item{detectdists}{detection distances}

\item{gpsdists}{distances from feature from collared animals}

\item{truncdet}{truncation distance of distance sampling}

\item{truncden}{truncation distance of collaring study}
}
\value{
combined likelihood evaluated with HR detection function
}
\description{
Combined likelihood with hazard rate detection and half normal gradient
}
