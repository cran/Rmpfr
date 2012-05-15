### R code from vignette source 'RmpfrArith.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width=75)
try(Mlibrary(Rmpfr))
stopifnot(require("Rmpfr"))


###################################################
### code chunk number 2: some-groups
###################################################
getGroupMembers("Arith")
getGroupMembers("Compare")
getGroupMembers("Math")


###################################################
### code chunk number 3: Matrix-ex
###################################################
head(x <- mpfr(0:7, 64)/7)
mx <- x ; dim(mx) <- c(4,2)
mx[ 1:3, ] + c(1,10,100)


###################################################
### code chunk number 4: mxv
###################################################
t(mx) %*% 10^(1:4)


###################################################
### code chunk number 5: crossprod
###################################################
crossprod(mx)


###################################################
### code chunk number 6: apply
###################################################
(s7 <- apply(7 * mx, 2, sum))


###################################################
### code chunk number 7: all.equal
###################################################
all.equal(s7, c(6,22), tol = 1e-40) # note the tolerance!


###################################################
### code chunk number 8: sessionInfo
###################################################
toLatex(sessionInfo())


###################################################
### code chunk number 9: Rmpfr-ver-show
###################################################
packageDescription("Rmpfr")


