\name{Fst}
\alias{Fst}
\title{F-Statistics}
\description{
  This function computes the \eqn{F_{IT}}{Fit}, \eqn{F_{ST}}{Fst} and
  \eqn{F_{IS}}{Fis} for each locus in the data.
}
\usage{
Fst(x, pop = NULL)
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
  \item{pop}{a vector or factor giving the population assignment of each
    row of \code{x}, or a single numeric value specifying which column
    of \code{x} to use as population indicator. By default, the column
    labelled \code{"population"} is used.}
}
\details{
  The formulae in Weir and Cockerham are used for each allele, and then
  averaged within each locus over the different alleles as suggested by
  these authors.
}
\value{
  a matrix with genes (loci) as rows and the three \emph{F}-statistics
  as columns.
}
\note{
  This function has not been tested thoroughly, so the results need to
  be taken with some caution.
}
\references{
  Weir, B. S. and Cockerham, C. C. (1984) Estimating \emph{F}-statistics
  for the analysis of population structure. \emph{Evolution}, \bold{38},
  1358--1370.
}
\author{Emmanuel Paradis}
\examples{
}
\keyword{htest}
