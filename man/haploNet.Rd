\name{haploNet}
\alias{haploNet}
\alias{plot.haploNet}
\title{Haplotype Networks}
\description{
  \code{haploNet} computes an haplotype network. There is a plot method.
}
\usage{
haploNet(h)
\method{plot}{haploNet}(x, size = 1, col = "black", bg = "white",
             col.link = "black", lwd = 1, lty = 1, pie = NULL,
             labels = TRUE, scale.ratio = 1, legend = FALSE,
             fast = FALSE, ...)
}
\arguments{
  \item{h}{an object of class \code{"haplotype"}.}
  \item{x}{an object of class \code{"haploNet"}.}
  \item{size}{a numeric vector giving the diameter of the circles
    representing the haplotypes: this is in the same unit than the links
    and eventually recycled.}
  \item{col}{a character vector specifying the colours of the circles;
    eventually recycled.}
  \item{bg}{a character vector specifying the colours of the background
    of the circles; eventually recycled.}
  \item{col.link}{a character vector specifying the colours of the links;
    eventually recycled.}
  \item{lwd}{a numeric vector giving the width of the links; eventually
    recycled.}
  \item{lty}{idem for the line types.}
  \item{pie}{a matrix used to draw pie charts for each haplotype; it
    must have as many rows than there is haplotypes.}
  \item{labels}{a logical specifying whether to identify the haplotypes
    with Roman numerals (the default).}
  \item{scale.ratio}{the ratio of the scale of the links representing
    the number of steps on the scale of the circles representing the
    haplotypes. It may be needed to give a value greater than one to
    avoid overlapping circles.}
  \item{legend}{a logical specifying whether to draw the legend, or a
    vector of length two giving the coordinates where to draw the
    legend; \code{FALSE} by default. If \code{TRUE}, the user is asked
    to click where to draw the legend.}
  \item{fast}{a logical specifying whether to optimize the spacing of
    the circles; \code{FALSE} by default.}
  \item{\dots}{further arguments passed to \code{plot}.}
}
\value{
  \code{haploNet} returns an object of class \code{"haploNet"} which is
  a matrix where each row represents a link in the network, the first
  and second columns give the number of linked haplotypes, the third
  column, named \code{"step"}, gives the number of steps in this link,
  and the fourth column, named \code{"Prob"}, gives the probability of a
  parsimonious link as given by Templeton et al. (1992).
}
\author{Emmanuel Paradis}
\references{
  Templeton, A. R., Crandall, K. A. and Sing, C. F. (1992) A cladistic
  analysis of phenotypic association with haplotypes inferred from
  restriction endonuclease mapping and DNA sequence data. III. Cladogram
  estimation. \emph{Genetics}, \bold{132}, 619--635.
}
\seealso{
  \code{\link{haplotype}}
}
\keyword{models}
\keyword{hplot}
