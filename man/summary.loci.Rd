\name{summary.loci}
\alias{summary.loci}
\alias{print.summary.loci}
\alias{plot.summary.loci}
\alias{print.loci}
\title{Print and Summaries of Loci Objects}
\description{
  These functions print and summarize table of alleles and loci (objects
  of class \code{"loci"}.
}
\usage{
\method{print}{loci}(x, details = FALSE, ...)
\method{summary}{loci}(object, ...)
\method{print}{summary.loci}(x, ...)
\method{plot}{summary.loci}(x, loci, what = "both", layout = 1, col = c("blue", "red"), ...)
}
\arguments{
  \item{x, object}{an object of class \code{"loci"} or \code{"genind"},
    a data frame, a factor, or a vector.}
  \item{details}{a logical value: if \code{TRUE} the data are printed as
    a data frame; the default is \code{FALSE}.}
  \item{loci}{the loci (genes) to be plotted. By default, all loci are
    plotted}
  \item{what}{the frequencies to be plotted. Three choices are possible:
    \code{"alleles"}, \code{"genotypes"}, and \code{"both"} (the default),
    or any unambiguous abbreviations.}
  \item{layout}{the number of graphs to be plotted simultaneously.}
  \item{col}{the colours used the barplots.}
  \item{\dots}{further arguments to be passed to or from other methods.}
}
\details{
  Genotypes not observed in the data frame are not counted.
}
\value{
  \code{summary.loci} returns a list with the genes as names and each
  element made a list with two vectors \code{"genotype"} and
  \code{"allele"} with the frequencies (numbers) of genotypes and
  alleles, respectively. The names of these two vectors are the observed
  genotypes and alleles.

  \code{print} and \code{plot} methods return NULL.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.loci}}, \code{\link{getAlleles}}
}
\keyword{manip}
\keyword{hplot}
