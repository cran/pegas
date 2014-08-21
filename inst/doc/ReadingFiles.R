### R code from vignette source 'ReadingFiles.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ReadingFiles.Rnw:30-31
###################################################
options(width=60)


###################################################
### code chunk number 2: ReadingFiles.Rnw:105-110
###################################################
library(pegas)
x <- read.loci("toto", header = FALSE)
x
print(x, details = TRUE)
class(x)


###################################################
### code chunk number 3: ReadingFiles.Rnw:126-129
###################################################
y <- read.loci("titi", header = FALSE, allele.sep = "-")
print(y, details = TRUE)
identical(x, y)


###################################################
### code chunk number 4: ReadingFiles.Rnw:134-135
###################################################
args(read.loci)


###################################################
### code chunk number 5: ReadingFiles.Rnw:162-163
###################################################
print(read.loci("tutu", FALSE), TRUE)


###################################################
### code chunk number 6: ReadingFiles.Rnw:179-182
###################################################
X <- read.loci("tyty")
print(X, TRUE)
summary(X)


###################################################
### code chunk number 7: ReadingFiles.Rnw:200-203
###################################################
z <- read.loci("tata", loci.sep = "\t", col.loci = 2:3, col.pop = 4, row.names = 1)
z
print(z, details = TRUE)


###################################################
### code chunk number 8: ReadingFiles.Rnw:209-210
###################################################
getAlleles(z)


###################################################
### code chunk number 9: ReadingFiles.Rnw:215-216
###################################################
attr(z, "locicol")


###################################################
### code chunk number 10: ReadingFiles.Rnw:221-222
###################################################
str(z)


###################################################
### code chunk number 11: ReadingFiles.Rnw:230-231
###################################################
args(read.vcf)


###################################################
### code chunk number 12: ReadingFiles.Rnw:273-277
###################################################
file.copy(system.file("files/nancycats.dat", package = "adegenet"), getwd())
file.copy(system.file("files/nancycats.gtx", package = "adegenet"), getwd())
file.copy(system.file("files/nancycats.gen", package = "adegenet"), getwd())
file.copy(system.file("files/nancycats.str", package = "adegenet"), getwd())


###################################################
### code chunk number 13: ReadingFiles.Rnw:283-286
###################################################
A <- read.fstat("nancycats.dat", quiet = TRUE)
B <- read.genetix("nancycats.gtx", quiet = TRUE)
C <- read.genepop("nancycats.gen", quiet = TRUE)


###################################################
### code chunk number 14: ReadingFiles.Rnw:294-295
###################################################
D <- read.structure("nancycats.str", onerowperind=FALSE, n.ind=237, n.loc=9, col.lab=1, col.pop=2, ask=FALSE, quiet=TRUE)


###################################################
### code chunk number 15: ReadingFiles.Rnw:300-302
###################################################
identical(A@tab, C@tab)
identical(B@tab, D@tab)


###################################################
### code chunk number 16: ReadingFiles.Rnw:309-310
###################################################
unlink(c("nancycats.dat", "nancycats.gtx", "nancycats.gen", "nancycats.str"))


