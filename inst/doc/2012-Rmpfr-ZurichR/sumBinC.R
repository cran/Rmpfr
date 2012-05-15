### R code from vignette source 'sumBinC.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width=75)
try(Mlibrary(Rmpfr))
stopifnot(require("Rmpfr"))


###################################################
### code chunk number 2: sumBinom-R
###################################################
sumBinom <- function(n, f, n0=0, ...) {
  k <- n0:n
  sum( choose(n, k) * (-1)^k * f(k, ...))
}
## and the same for a whole *SET* of  n  values:
sumBin.all.R <- function(n, f, n0=0, ...)
   sapply(n, sumBinom, f=f, n0=n0, ...)


###################################################
### code chunk number 3: sumBinomMpfr
###################################################
sumBinomMpfr


###################################################
### code chunk number 4: sumBin.all
###################################################
sumBin.all <- function(n, f, n0=0, precBits = 256, ...)
{
  N <- length(n)
  precBits <- rep(precBits, length = N)
  ll <- lapply(seq_len(N), function(i)
           sumBinomMpfr(n[i], f, n0=n0, precBits=precBits[i], ...))
  sapply(ll, as, "double")
}


###################################################
### code chunk number 5: sqrt-ex
###################################################
nn <- 5:80
system.time(res.R   <- sumBin.all.R(nn, f = sqrt)) ## instant!
system.time(resMpfr <- sumBin.all  (nn, f = sqrt)) ## ~2 seconds


###################################################
### code chunk number 6: sqrt-ex-2
###################################################
matplot(nn, cbind(res.R, resMpfr), type = "l", lty=1,
        ylim = extendrange(resMpfr, f = 0.25), xlab = "n",
        main = "sumBinomMpfr(n, f = sqrt)  vs.  R double precision")
legend("topleft", leg=c("double prec.", "mpfr"), lty=1, col=1:2, bty = "n")


