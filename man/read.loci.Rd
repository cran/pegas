\name{read.loci}
\alias{read.loci}
\title{Read Allelic Data Files}
\description{
  This function reads allelic data from a text file: rows are
  individuals, and each column is an allele. By default, the header of
  the file gives the locus names. If one column is labelled
  ``population'', it is taken as a population variable.
}
\usage{
read.loci(file, header = TRUE, loci.sep = " ", allele.sep = "/",
          col.pop = "none", col.loci = NULL, skip = 0)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a quoted string.}
  \item{header}{a logical specifying whether the first line of the data
    file gives the names of the loci (\code{TRUE} by default).}
  \item{loci.sep}{the character(s) separating the loci (columns) in the
    data file (a space by default).}
  \item{allele.sep}{the character(s) separating the alleles for each locus
    in the data file (a slash by default).}
  \item{col.pop}{specifies whether one of the column of the data file
    identifies the population; default \code{"none"}, otherwise an
    integer giving the number of the column.}
  \item{col.loci}{a vector of integers specifying the indices of the
    columns that are loci. By default, all columns are taken as loci
    except one labelled "population", if present.}
  \item{skip}{an integer specifying the number of lines in the data file
    to be skipped before reading the data.}
}
\details{
  The rownames of the returned object identify the individual genotypes;
  they are either taken from the data file if present, or given the
  values \code{"1"}, \code{"2"}, \dots

  In the returned object, alleles are separated by \code{"/"}, even if
  it is not the case in the data file.
}
\value{
  a data frame with class \code{c("loci", "data.frame")}. Details on the
  structure can be found in
  \url{http://ape.mpl.ird.fr/pegas/DefinitionDataClassesPegas.pdf}
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.gtx}}, \code{\link{write.loci}}
}
\keyword{IO}
