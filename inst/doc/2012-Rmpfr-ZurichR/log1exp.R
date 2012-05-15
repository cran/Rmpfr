### R code from vignette source 'log1exp.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(SweaveHooks= list(fig=function() par(mar=c(5.1, 4.1, 1.1, 2.1))),
        width = 75,
        digits = 7, # <-- here, keep R's default!
        prompt = "> ",
        continue="  ")
try(Mlibrary(Rmpfr))
stopifnot(require("Rmpfr"))


###################################################
### code chunk number 2: log1-exp-curve-10
###################################################
curve(-log(1 - exp(-x)), 0, 10)


###################################################
### code chunk number 3: log1-exp-curve-log
###################################################
curve(-log(1 - exp(-x)),  0, 50, log="y")


###################################################
### code chunk number 4: p-log1mexp--def
###################################################
## Nicer
p.log1mexp <- function( n=400, col=2, colA = "blue3", lwd= 1.5, ...) {
    cc <- curve(-log(1 - exp(-x)), 0, 50, 0, log="y", yaxt="n", n=n, col=col, lwd=lwd, ...)
    sfsmisc::eaxis(2, labels=FALSE)
    sfsmisc::eaxis(2, axTicks(2, log=par("ylog"))[c(TRUE,FALSE)])
    arrows(38, 1e-15,38, 4e-17, col=colA, lwd=2, length= 1/12)
    text  (38, 1e-15, "early underflow to 0", col=colA, adj=c(.2,-.5))
    invisible(cc)
}


###################################################
### code chunk number 5: log1-exp-log-niceer
###################################################
p.log1mexp()


###################################################
### code chunk number 6: log1exp-ex
###################################################
x <- -40:-35
    -log(1 - exp(x))
log(-log(1 - exp(x)))# --> -Inf values
## ok, how about more accuracy
x. <- mpfr(x, 120)
log(-log(1 - exp(x.)))# aha... looks perfect now


###################################################
### code chunk number 7: plot-log1exp-ex
###################################################
x <- seq(-40, -20, by = .5)
plot(x,x, type="n", ylab="", ann=FALSE)
lines(x, log(-log(1 - exp(x))), type = "o", col = "purple", lwd=3, cex = .6)


###################################################
### code chunk number 8: plot-log1exp-mpfr
###################################################
x <- seq(-40, -20, by = .5)
plot(x,x, type="n", ylab="", ann=FALSE)
lines(x, log(-log(1 - exp(x))), type = "o", col = "purple", lwd=3, cex = .6)
x. <- mpfr(x, 120)
lines(x, log(-log(1 - exp(x.))), col=2, lwd=1.5)


