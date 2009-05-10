## IO.R (2009-05-10)

##   Input/Ouput

## Copyright 2009 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../COPYING for licensing issues.

read.loci <-
    function(file, header = TRUE, loci.sep = " ", allele.sep = "/",
             col.pop = "none", col.loci = NULL, skip = 0)
{
    res <- read.table(file = file, header = header, sep = loci.sep,
                      skip = skip, colClasses = "factor")
    as.loci.data.frame(res, allele.sep = allele.sep, col.pop = col.pop,
                       col.loci = col.loci)
}

read.gtx <- function(file)
{
    x <- scan(file, what = "", sep = "\n", quiet = TRUE)
    last <- x[length(x)]
    nloci <- length(strsplit(substr(last, 12, nchar(last)), " +")[[1]])
    npop <- as.integer(strsplit(x[2], " +")[[1]][1])
    ## the number of individuals is found easily:
    n <- length(x) - 2L * (1L + nloci + npop)
    loci.nms <- x[seq(from = 3, by = 2, length.out = nloci)]
    j <- 2L * nloci + 3L
    k <- 1L
    pop.nms <- character(npop)
    keep <- logical(length(x))
    pop <- integer(n)
    for (i in 1:npop) {
        pop.nms[i] <- x[j]
        m <- as.integer(x[j + 1L])
        pop[k:(k + m - 1L)] <- i
        k <- k + m
        j <- j + 2L
        keep[j:(j + m - 1L)] <- TRUE
        j <- j + m
    }
    x <- gsub("^ +", "", x[keep])
    x <- matrix(unlist(strsplit(x, " +")), n, nloci + 1L, byrow = TRUE)
    obj <- as.data.frame(x[, -1])
    for (i in 1:ncol(obj)) {
        levels(obj[, i]) <-
            paste(substr(levels(obj[, i]), 1, 3),
                  substr(levels(obj[, i]), 4, 6), sep = "/")
    }
    dimnames(obj) <- list(x[, 1], loci.nms)
    pop.nms <- gsub("^ +", "", pop.nms)
    pop.nms <- gsub(" +$", "", pop.nms)
    class(pop) <- "factor"
    levels(pop) <- pop.nms
    obj$population <- pop
    attr(obj, "locicol") <- 1:nloci
    obj <- .check.order.alleles(obj)
    class(obj) <- c("loci", "data.frame")
    obj
}

write.loci <- function(x, file = "", loci.sep = " ", allele.sep = "/")
{
    if (allele.sep != "/") {
        for (i in 1:attr(x, "locicol"))
            levels(x[[i]]) <- gsub("/", allele.sep, levels(x[[i]]))
    }
    write.table(x, file = file, sep = loci.sep)
}
