% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SA24204141_R.R
\name{logit}
\alias{logit}
\title{logistic function}
\usage{
logit(x)
}
\arguments{
\item{x}{independent variable (Matrix)}
}
\value{
the value of the logistic function corresponding to x
}
\description{
Convert real numbers to numbers between (0,1) using logistic functions
}
\examples{
\dontrun{
 tMu <- matrix(c(0.3,-3,1,-0.5), 2, 2)
 logit(tMu)
}
}
