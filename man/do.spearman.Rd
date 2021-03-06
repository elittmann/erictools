% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do.spearman.R
\name{do.spearman}
\alias{do.spearman}
\title{Do Spearman Function
This function lets you perform spearman correlations in dplyr syntax.}
\usage{
do.spearman(x, y)
}
\arguments{
\item{x}{Numeric Predictor}

\item{y}{Numeric Outcome}
}
\description{
Do Spearman Function
This function lets you perform spearman correlations in dplyr syntax.
}
\examples{
mtcars \%>\%
group_by(cyl) \%>\%
 summarize(pvalue=as.numeric(do.spearman(wt,mpg)[3]),
           rho=as.numeric(do.spearman(wt,mpg)[4]))
}
\keyword{spearman}
