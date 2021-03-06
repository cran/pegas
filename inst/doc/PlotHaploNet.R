### R code from vignette source 'PlotHaploNet.Rnw'

###################################################
### code chunk number 1: PlotHaploNet.Rnw:23-24
###################################################
options(width = 80, prompt = "> ")


###################################################
### code chunk number 2: PlotHaploNet.Rnw:98-100
###################################################
library(pegas) # loads also ape
data(woodmouse)


###################################################
### code chunk number 3: PlotHaploNet.Rnw:109-114
###################################################
set.seed(10)
x <- woodmouse[sample.int(nrow(woodmouse), 80, TRUE), ]
region <- rep(c("regA", "regB"), each = 40)
pop <- rep(paste0("pop", 1:4), each = 20)
table(region, pop)


###################################################
### code chunk number 4: PlotHaploNet.Rnw:120-125
###################################################
h <- haplotype(x)
h
d <- dist.dna(h, "N")
nt <- rmst(d, quiet = TRUE)
nt


###################################################
### code chunk number 5: PlotHaploNet.Rnw:130-131
###################################################
plot(nt)


###################################################
### code chunk number 6: PlotHaploNet.Rnw:136-137
###################################################
plot(nt, fast = TRUE)


###################################################
### code chunk number 7: PlotHaploNet.Rnw:144-145
###################################################
plot(nt, threshold = c(1, 14))


###################################################
### code chunk number 8: PlotHaploNet.Rnw:156-157
###################################################
args(pegas:::plot.haploNet)


###################################################
### code chunk number 9: PlotHaploNet.Rnw:180-181
###################################################
(sz <- summary(h))


###################################################
### code chunk number 10: PlotHaploNet.Rnw:187-188
###################################################
(nt.labs <- attr(nt, "labels"))


###################################################
### code chunk number 11: PlotHaploNet.Rnw:194-196
###################################################
sz <- sz[nt.labs]
plot(nt, size = sz)


###################################################
### code chunk number 12: PlotHaploNet.Rnw:203-205
###################################################
(R <- haploFreq(x, fac = region, haplo = h))
(P <- haploFreq(x, fac = pop, haplo = h))


###################################################
### code chunk number 13: PlotHaploNet.Rnw:211-213
###################################################
R <- R[nt.labs, ]
P <- P[nt.labs, ]


###################################################
### code chunk number 14: PlotHaploNet.Rnw:219-220
###################################################
plot(nt, size = sz, pie = R, legend = c(-25, 30))


###################################################
### code chunk number 15: PlotHaploNet.Rnw:223-224
###################################################
plot(nt, size = sz, pie = P, legend = c(-25, 30))


###################################################
### code chunk number 16: PlotHaploNet.Rnw:289-304
###################################################
par(xpd = TRUE)
size <- c(1, 3, 5, 10)
x <- c(0, 5, 10, 20)

plot(0, 0, type="n", xlim=c(-2, 30), asp=1, bty="n", ann=FALSE)
other.args <- list(y = -5, inches = FALSE, add = TRUE,
                   bg = rgb(1, 1, 0, .3))
o <- mapply(symbols, x = x, circles = sqrt(size / pi),
            MoreArgs = other.args)
other.args$y <- 5
o <- mapply(symbols, x = x, circles = size / 2,
            MoreArgs = other.args)
text(x, -1, paste("size =", size), font = 2, col = "blue")
text(30, -5, expression("circles = "*sqrt(size / pi)))
text(30, 5, "circles = size / 2")


###################################################
### code chunk number 17: PlotHaploNet.Rnw:313-322
###################################################
x <- c(0, 6, 13, 25)
plot(0, 0, type="n", xlim=c(-2, 30), asp=1, bty="n", ann=FALSE)
other.args$y <- 0
o <- mapply(symbols, x = x, circles = size/2, MoreArgs = other.args)
other.args$col <- "black"
other.args$add <- other.args$inches <- NULL
o <- mapply(pegas:::square, x = x, size = size, MoreArgs = other.args)
o <- mapply(pegas:::diamond, x = x, size = size, MoreArgs = other.args)
text(x, -7, paste("size =", size), font = 2, col = "blue")


###################################################
### code chunk number 18: PlotHaploNet.Rnw:366-368
###################################################
plot(nt)
mutations(nt, 18, x = -8.9, y = 16.3, data = h)


###################################################
### code chunk number 19: PlotHaploNet.Rnw:377-380
###################################################
plot(nt)
mutations(nt, 18, x = -8.9, y = 16.3, data = h)
mutations(nt, 18, x = 10, y = 17, data = h, style = "s")


###################################################
### code chunk number 20: PlotHaploNet.Rnw:399-400
###################################################
names(getHaploNetOptions())


###################################################
### code chunk number 21: PlotHaploNet.Rnw:406-407
###################################################
plot(nt, size = 2)


###################################################
### code chunk number 22: PlotHaploNet.Rnw:410-414
###################################################
setHaploNetOptions(haplotype.inner.color = "#CCCC4D",
                   haplotype.outer.color = "#CCCC4D",
                   show.mutation = 3, labels = FALSE)
plot(nt, size = 2)


###################################################
### code chunk number 23: PlotHaploNet.Rnw:417-422
###################################################
setHaploNetOptions(haplotype.inner.color = "blue",
                   haplotype.outer.color = "blue",
                   show.mutation = 1)
par(bg = "yellow3")
plot(nt, size = 2)


###################################################
### code chunk number 24: PlotHaploNet.Rnw:425-429
###################################################
setHaploNetOptions(haplotype.inner.color = "navy",
                   haplotype.outer.color = "navy")
par(bg = "lightblue")
plot(nt, size = 2)


