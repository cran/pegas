\name{read.loci}
\alias{read.loci}
\title{Read Allelic Data Files}
\description{
  This function reads allelic data from a text file: rows are
  individuals, and columns are loci and optional variables. By default,
  the first line of the file gives the locus names. If one column is
  labelled ``population'', it is taken as a population variable.
}
\usage{
read.loci(file, header = TRUE, loci.sep = "", allele.sep = "/",
          col.pop = "none", col.loci = NULL, ...)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a quoted string.}
  \item{header}{a logical specifying whether the first line of the data
    file gives the names of the loci (\code{TRUE} by default).}
  \item{loci.sep}{the character(s) separating the loci (columns) in the
    data file (a white space by default).}
  \item{allele.sep}{the character(s) separating the alleles for each locus
    in the data file (a forward slash by default).}
  \item{col.pop}{specifies whether one of the column of the data file
    identifies the population; default \code{"none"}, otherwise an
    integer giving the number of the column.}
  \item{col.loci}{a vector of integers specifying the indices of the
    columns that are loci. By default, all columns are taken as loci
    except one labelled "population", if present.}
  \item{\dots}{further arguments passed to \code{read.table}.}
}
\details{
  The rownames of the returned object identify the individual genotypes;
  they are either taken from the data file if present, or given the
  values \code{"1"}, \code{"2"}, \dots Similarly for the colnames: if
  absent in the file (in which case \code{header = FALSE} must be set),
  they are given the values \code{"V1"}, \code{"V2"}, \dots

  In the returned genotypes, alleles are separated by \code{"/"}, even
  if it is not the case in the data file.

  The vignette ``Reading Genetic Data Files Into R with \pkg{adegenet}
  and \pkg{pegas}'' explains how to read various file formats including
  Excel files (type \code{vignette("ReadingFiles")} in R).

}
\value{
  A data frame with class \code{c("loci", "data.frame")}. It is a data
  frame with an attribute \code{"locicol"} specifying the columns that
  must be treated as loci. The latter are factors. The other columns can
  be of any type.

  Details on the structure can be found in
  \url{http://ape.mpl.ird.fr/pegas/DefinitionDataClassesPegas.pdf}
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.gtx}}, \code{\link{write.loci}},
  \code{\link{summary.loci}}
}
\keyword{IO}
