\name{haplotype}
\alias{haplotype}
\alias{plot.haplotype}
\title{Haplotype Extraction and Frequencies}
\description{
  \code{haplotype} extracts the haplotypes from a set of DNA
  sequences. The result can be plotted with the appropriate function.
}
\usage{
haplotype(x, labels = NULL)
\method{plot}{haplotype}(x, ...)
}
\arguments{
  \item{x}{a set of DNA sequences (as an object of class
    \code{"DNAbin"}, or an object of class \code{"haplotype"}.}
  \item{labels}{a vector of character strings used as names for the rows
    of the returned object. By default, Roman numbers are given.}
  \item{\dots}{further arguments passed to \code{\link[graphics]{barplot}}.}
}
\value{
  \code{haplotype} returns an object of class \code{c("haplotype",
    "DNAbin")} which is an object of class \code{"DNAbin"} with an
  attribute \code{"index"} identifying the index of each observation
  that share the same haplotype.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link[ape]{DNAbin}} for manipulation of DNA sequences in R.
}
\keyword{manip}
\keyword{hplot}
