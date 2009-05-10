## Fst.R (2009-05-12)

##   F-Statistics

## Copyright 2009 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../COPYING for licensing issues.

Fst <- function(x, pop = NULL)
{
    if (is.null(pop)) {
        pop <- x$population
        if (is.null(pop)) stop("no 'population' column in x")
    } else {
        pop <- if (is.numeric(pop)) x[, pop] else factor(pop)
    }

    ALLELES <- getAlleles(x)

    npop <- length(levels(pop)) # 'r' in Weir & Cockerham's notation
    ngene <- length(attr(x, "locicol"))
    nBYpop <- table(pop)
    N <- nrow(x)
    nbar <- N/npop

    ## get the frequencies in each population for each locus
    p <- vector("list", ngene)
    for (j in 1:ngene)
        p[[j]] <- matrix(0, length(ALLELES[[j]]), npop)
    h <- p
    for (i in 1:npop) {
        s <- summary(x[as.integer(pop) == i, ]) # levels are preserved
        for (j in 1:ngene) {
            tmp <- s[[j]]
            p[[j]][, i] <- tmp$allele
            allel <- names(tmp$allele)
            genot <- names(tmp$genotype)
            for (k in seq_along(allel)) {
                for (l in seq_along(genot)) {
                    ag <- unlist(strsplit(genot[l], "/"))
                    if (sum(ag %in% allel[k]) == 1) # assume diploidy I guess
                        h[[j]][k, i] <- h[[j]][k, i] + tmp$genotype[l]
                }
            }
        }
    }
    ## 'p' is a list with, for each locus, a matrix with alleles as rows and populations ar columns, and its entries are the counts
    ## 'h' is the same but with the number of heterozygotes

    nC <- (N - sum(nBYpop^2)/N)/(npop - 1)

    obj <- matrix(0, ngene, 3)
    for (j in 1:ngene) {
        ptild <- rowSums(p[[j]])
        #browser()
        pbar <- ptild/N # for each allele in the locus
        s2 <- (rowSums(p[[j]] - pbar))^2/((npop - 1) * nbar)
        ## s2 <- nBYpop * (ptild - pbar)^2/((npop - 1) * nbar)
        hbar <- rowSums(h[[j]])/N # id.
        ## hbar <- nBYpop * rowSums(h[[j]])/N # id.
        A <- pbar * (1 - pbar) - (npop - 1) * s2/npop
        a <- (nbar/nC) * (s2 - (A - hbar/4)/(nbar - 1))
        b <- nbar * (A - 2 * (nbar - 1) * hbar/(4 * nbar))/(nbar - 1)
        c <- hbar/2
        obj[j, 1] <- 1 - sum(c)/sum(a + b + c)
        obj[j, 2] <- sum(a)/sum(a + b + c)
        obj[j, 3] <- 1 - sum(c)/sum(b + c)
    }
    dimnames(obj) <- list(names(x)[attr(x, "locicol")],
                          c("Fit", "Fst", "Fis"))
    obj
}
