\name{as.loci}
\alias{as.loci}
\alias{as.loci.genind}
\alias{as.loci.data.frame}
\alias{as.loci.factor}
\alias{as.loci.vector}
\alias{genind2loci}
\alias{loci2genind}
\title{Conversion Among Allelic Data Classes}
\description{
  These functions do conversion among different allelic data classes.
}
\usage{
as.loci(x, ...)
\method{as.loci}{genind}(x, ...)
genind2loci(x)
\method{as.loci}{data.frame}(x, allele.sep = "/", col.pop = "none", col.loci = NULL, ...)
loci2genind(x)
\method{as.loci}{factor}(x, allele.sep = "/", ...)
\method{as.loci}{vector}(x, allele.sep = "/", ...)
}
\arguments{
  \item{x}{an object of class \code{"loci"} or \code{"genind"}, a data
    frame, a factor, or a vector.}
    \item{allele.sep}{the character(s) separating the alleles for each locus
    in the data file (a slash by default).}
  \item{col.pop}{specifies whether one of the column of the data file
    identifies the population; default \code{"none"}, otherwise an
    integer giving the number of the column.}
  \item{col.loci}{a vector of integers specifying the indices of the
    columns that are loci. By default, all columns are taken as loci
    except one labelled "population", if present.}
  \item{...}{further arguments to be passed to or from other methods.}
}
\details{
  \code{genind2loci(x)} and \code{as.loci(x)} are the same if \code{x}
  is of class \code{"genind"}.
}
\value{
  a data frame with class \code{c("loci", "data.frame")} for
  \code{as.loci} and \code{genind2loci}; an object of class
  \code{"genind"} for \code{loci2genind}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.loci}}
}
\keyword{IO}
