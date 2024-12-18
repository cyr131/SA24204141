% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{statCpp}
\alias{statCpp}
\title{ST correlation matrix function}
\usage{
statCpp(t, a)
}
\arguments{
\item{t}{Specify matrix dimensions}

\item{a}{Specify non-diagonal element}
}
\value{
ST correlation matrix with t_dimensions and non-diagonal element a
}
\description{
generate stationary working correlation matrix
}
\examples{
\dontrun{
t <- 6
a <- c(0.5,0.4,0.3,-0.2,-0.1)
statCpp(t,a)
}
}
