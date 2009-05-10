## summary.loci.R (2009-05-10)

##   Print and Summaries of Loci Objects

## Copyright 2009 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../COPYING for licensing issues.

getPloidy <- function(x)
    unlist(lapply(x[, attr(x, "locicol"), drop = FALSE],
                  function(x) sum(charToRaw(levels(x)[1]) == 47) + 1))

getAlleles <- function(x)
    lapply(x[, attr(x, "locicol"), drop = FALSE],
           function(x) unique(unlist(strsplit(levels(x), "/"))))

getGenotypes <- function(x)
    lapply(x[, attr(x, "locicol"), drop = FALSE], levels)

print.loci <- function(x, details = FALSE, ...)
{
    if (details) print.data.frame(x) else {
        n <- dim(x)
        nloci <- length(attr(x, "locicol"))
        cat("Allelic data frame:", n[1], "individuals\n")
        cat("                   ", nloci, "loci\n")
        nav <- n[2] - nloci
        if (nav) cat("                   ", nav, "additional variables\n")
    }
}

summary.loci <- function(object, ...)
{
    L <- attr(object, "locicol")
    ans <- vector("list", length(L))
    names(ans) <- names(object[L])
    for (i in L) {
        geno <- levels(object[, i])
        alle <- strsplit(geno, "/")
        unialle <- sort(unique(unlist(alle)))
        l <- tabulate(object[, i], length(geno))
        names(l) <- geno
        tab <- matrix(0, length(unialle), length(geno),
                      dimnames = list(unialle, geno))
        for (j in seq_along(alle))
            for (k in alle[[j]]) tab[k, j] <- tab[k, j] + 1
        ans[[i]] <- list(genotype = l, allele = drop(tab %*% l))
    }
    class(ans) <- "summary.loci"
    ans
}

print.summary.loci <- function(x, ...)
{
    nms <- names(x)
    for (i in 1:length(x)) {
        cat("Locus", nms[i], ":\n")
        cat("-- Genotype frequencies:\n")
        print(x[[i]][[1]])
        cat("-- Allele frequencies:\n")
        print(x[[i]][[2]])
        cat("\n")
    }
}
