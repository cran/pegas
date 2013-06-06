### R code from vignette source 'ReadingFiles.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ReadingFiles.Rnw:29-30
###################################################
options(width=60)


###################################################
### code chunk number 2: ReadingFiles.Rnw:97-102
###################################################
library(pegas)
x <- read.loci("toto", header = FALSE)
x
print(x, details = TRUE)
class(x)


###################################################
### code chunk number 3: ReadingFiles.Rnw:118-121
###################################################
y <- read.loci("titi", header = FALSE, allele.sep = "-")
print(y, details = TRUE)
identical(x, y)


###################################################
### code chunk number 4: ReadingFiles.Rnw:126-127
###################################################
args(read.loci)


###################################################
### code chunk number 5: ReadingFiles.Rnw:152-153
###################################################
print(read.loci("tutu", FALSE), TRUE)


###################################################
### code chunk number 6: ReadingFiles.Rnw:171-174
###################################################
z <- read.loci("tata", loci.sep = "\t", col.loci = 1:2, col.pop = 3, row.names = 1)
z
print(z, details = TRUE)


###################################################
### code chunk number 7: ReadingFiles.Rnw:180-181
###################################################
getAlleles(z)


###################################################
### code chunk number 8: ReadingFiles.Rnw:186-187
###################################################
attr(z, "locicol")


###################################################
### code chunk number 9: ReadingFiles.Rnw:192-193
###################################################
str(z)


###################################################
### code chunk number 10: ReadingFiles.Rnw:219-223
###################################################
file.copy(system.file("files/nancycats.dat", package = "adegenet"), getwd())
file.copy(system.file("files/nancycats.gtx", package = "adegenet"), getwd())
file.copy(system.file("files/nancycats.gen", package = "adegenet"), getwd())
file.copy(system.file("files/nancycats.str", package = "adegenet"), getwd())


###################################################
### code chunk number 11: ReadingFiles.Rnw:229-232
###################################################
A <- read.fstat("nancycats.dat", quiet = TRUE)
B <- read.genetix("nancycats.gtx", quiet = TRUE)
C <- read.genepop("nancycats.gen", quiet = TRUE)


###################################################
### code chunk number 12: ReadingFiles.Rnw:240-241
###################################################
D <- read.structure("nancycats.str", onerowperind=FALSE, n.ind=237, n.loc=9, col.lab=1, col.pop=2, ask=FALSE, quiet=TRUE)


###################################################
### code chunk number 13: ReadingFiles.Rnw:246-248
###################################################
identical(A@tab, C@tab)
identical(B@tab, D@tab)


###################################################
### code chunk number 14: ReadingFiles.Rnw:255-256
###################################################
unlink(c("nancycats.dat", "nancycats.gtx", "nancycats.gen", "nancycats.str"))


