### dgamma(): ----------------------- was part of  ./special-fun-ex.R -------------------
stopifnot(require("Rmpfr"))
require(sfsmisc)# -> eaxis(); relErrV()

(doExtras <- Rmpfr:::doExtras())
options(nwarnings = 50000, width = 99)

##                                      vvvvvvvvvvvvvvvv
## to enhance  |rel.Err| plots:  from ./special-fun-ex.R, also in ~/R/Pkgs/DPQ/tests/pow-tst.R }
drawEps.h <- function(p2 = -(53:51), side = 4, lty=3, lwd=2, col=adjustcolor(2, 1/2)) {
    abline(h = 2^p2, lty=lty, lwd=lwd, col=col)
    axis(side, las=2, line=-1, at = 2^p2,
         labels = as.expression(lapply(p2, function(p) substitute(2^E, list(E=p)))),
         col.axis = col, col=NA, col.ticks=NA)
}

(do.pdf <- !dev.interactive(orNone = TRUE))

if(do.pdf) pdf("special-fun-dgamma.pdf")

xe <- c(-2e5, -1e5, -2e4, -1e4, -2000, -1000, -500, -200, -100, -50, -20, -10)
(xe <- c(xe, -8:8, -rev(xe)))
two <- mpfr(2, 256)
## For centering at E[.], will use xP(x, shp) :
xP <- function(x, d) {
    ## cannot eliminate them, as for <mpfr> they are all finite ..
    ## x <- x[is.finite(x)]
    x - d*(x > d)
}
aEQformat <- function(xy, ...) format(xy, digits = 7, ...)
allEQ_0 <- function (target, current, ...)
    all.equal(target, current, tolerance = 0, formatFUN = aEQformat, ...)
stopIfNot <-
    if("allow.logical0" %in% names(formals(stopifnot))) { # experimental (MM only)
        stopifnot
    } else function(exprs, allow.logical0) stopifnot(exprs=exprs)
abs19 <- function(r) pmax(abs(r), 1e-19) # cut  |err| to positive {for log-plots}

for(shp in 2^c(-20, -3, -1:1, 4, 10, 14, 20, 50)) {
    cat("shape = 2^", log2(shp), ":\n-------------\n")
    d.dg  <- dgamma(xP(2 ^ xe, shp) -> x, shape=shp)
    m.dg  <- dgamma(xP(two^xe, shp), shape=shp)
    m.ldg <- dgamma(xP(two^xe, shp), shape=shp, log=TRUE)
    relE <- asNumeric(relErrV(m.dg, d.dg))
    ## Plots: do *not* observe any problems yet
    plot(x, relE, log="x", type="l",
         main = paste0("rel.Errors dgamma(., shape = 2^", log2(shp),")"))
    abline(h=0, col=adjustcolor("gray10", 1/2), lty=3, lwd=2)
    plot(x, abs19(relE), log="xy", type="l", ylim = pmax(4e-17, range(abs19(relE), finite=TRUE)))
    abline(h = 2^-(52:50), col=adjustcolor("red4",1/2), lty=3)
    ##
    stopIfNot(exprs = {
        !is.unsorted(xe)
        is.finite(m.dg)
        m.dg >= 0
        shp > 1  || all(diff(m.dg) <= 0)
        shp > 100|| all((m.dg > 0) >= (d.dg > 0))
        any(fin.d <- is.finite(d.dg))
        m.dg[!fin.d] > 1e300
        { cat("all.EQ(<mpfr>, <doubl>):", allEQ_0  (m.dg[fin.d], d.dg[fin.d]), "\n")
          shp > 100  ||                   all.equal(m.dg[fin.d], d.dg[fin.d],
                                                    tol = 1e-13) # 2.063241e-14
        }
        ## compare with log scale :
        if(any(pos.d <- m.dg > 0)) {
            cat("For non-0 <mpfr>-values;  all.EQ(log(d), d*(log)):",
                allEQ_0  (log(m.dg[pos.d]), m.ldg[pos.d]),"\n")
            ##
            all.equal(log(m.dg[pos.d]), m.ldg[pos.d], tol = 1e-14)
        } else TRUE
    })#, allow.logical0 = TRUE)
}

## NB:  dgamma(x, sh)  sh >= 1 calls
## --   dpois_raw(sh-1, x)   which then
## adds stirlerr(x.) to bd0(x., lambda) ;  where x. <- sh-1; lambda <- x
## bd0(x,L) ~= 0 iff x ~= L  <==> sh-1 ~= x  <==>  x+1 ~= sh

sh2x_gamma <- function(sh, nx, k = 12, f1 = 0.5, f2 = 1.25) {
    stopifnot(is.numeric(sh), length(sh) == 1, sh >= 0, length(k) == 1, k >= 3,
              f1 < f2, length(f1) == length(f2), length(f2) == 1)
    p2 <- 2^-(k:3)
    1 + sh* unique(sort(c(1-p2, 1, 1+p2, # <- values x very close to sh -- does *not* make any diff (????)
                          seq(f1, f2, length=nx))))
}

relEgamma <- function(sh, nx = 1001, k = 12, precBits = 256,
                      x = sh2x_gamma(sh, nx=nx, k=k)) {
    dg  <- dgamma(x, sh)
    dgM <- dgamma(mpfr(x,  precBits),
                  mpfr(sh, precBits))
    structure(cbind(x, relE = asNumeric(relErrV(dgM, dg))), shape=sh)
}

shs <- 1/32 + seq(6, 16, by = 1/8)
stopifnot(all(!is.whole(shs*2))) # shs are *not* half-integers
system.time(
    LrelE <- lapply(shs, relEgamma)
) # 7.5 sec

m.relE <- sapply(LrelE, function(m) m[,"relE"])
qrelE <- t(apply(abs(m.relE), 2, quantile, probs = c(.05, .25, .50, .75, .90, 1)))
##               ^^^^^^^^^^^

## Heureka! --- this shows quite a difference between R 4.3.3  and R-devel (R 4.4.0) !!

iS <- sort.list(qrelE[,"50%"], decreasing = TRUE)
cbind(shs, qrelE)[iS,]
## For R 4.3.3 :
##      shs             5%          25%          50%          75%          90%         100%
## 14.53125   9.410630e-15 9.815160e-15 1.023178e-14 1.065722e-14 1.092232e-14 1.138372e-14
## 15.03125   8.265317e-15 8.702900e-15 9.086072e-15 9.506915e-15 9.756106e-15 1.007928e-14
## 15.90625   6.799207e-15 7.137733e-15 7.611360e-15 8.057580e-15 8.343992e-15 8.670817e-15
## 13.53125   6.716182e-15 7.103502e-15 7.566360e-15 8.004966e-15 8.276645e-15 8.630780e-15
## 15.65625   6.031124e-15 6.389848e-15 6.803347e-15 7.261310e-15 7.527491e-15 8.031559e-15
## ..........
## ..........

myRversion <- paste(R.version.string, "--", osVersion)
if((mach <- Sys.info()[["machine"]]) != "x86_64")
    myRversion <- paste0(myRversion, "_", mach)
if(!capabilities("long.double"))
    myRversion <- paste0(myRversion, "_no_LDbl")
myRversion

rngP <- function(y, M = 1e-14) { yr <- range(y); if(yr[2] < M) yr[2] <- M; yr }

boxplot(abs19(m.relE), at = shs, log="y", ylim = c(7e-17, rngP(abs(m.relE))[2]), yaxt="n")
eaxis(2); drawEps.h(); mtext(myRversion, adj=1, cex=3/4)

matplot(shs, qrelE, type="l", log="y", yaxt="n", ylim = rngP(qrelE))
title("|relErr( dgamma(x, sh) |   for  x / sh  in [.5, 1.25]")
eaxis(2); drawEps.h(); mtext(myRversion, adj=1, cex=3/4)

## take *one* of these:
plot(abs(m.relE[, shs == 14.53125]), type="l", log="y", ylim = c(1e-16, 1.5e-14))
drawEps.h()

sh <- 14.53125
stopifnot(identical(sh, 465 / 32))

x14.5 <- sh2x_gamma(sh, nx = 21) # 21 points
xM <- mpfr(x14.5, 512)
dg1 <- stats::dgamma(x14.5, sh)
dgM <- Rmpfr::dgamma(xM,    sh)
cbind(x14.5, relE = asNumeric(relErrV(dgM, dg1))) # very "constant" ~=~  - 1e-14
                                        #            in R-devel   around  1e-16 !!
## try easier x:
sh <- 14.53125 ; stopifnot(identical(sh, 465 / 32))
x0 <- 1/4 + 8:20
xM <- mpfr(x0, 512)
dg1 <- stats::dgamma(x0, sh)
dgM <- Rmpfr::dgamma(xM, sh)
relE <- asNumeric(relErrV(dgM, dg1))
signif(cbind(x0, relE, abs(relE)), 4) # R <= 4.3.*: very "constant" ~=~  - 1e-14
## R-devel:                   |  no-long-double  == *same* numbers
##    x0       relE           |       relE
##  8.25  1.276e-16 1.276e-16 |  1.276e-16 1.276e-16
##  9.25  1.294e-16 1.294e-16 |  1.294e-16 1.294e-16
## 10.25 -1.408e-16 1.408e-16 | -1.408e-16 1.408e-16
## 11.25 -2.108e-17 2.108e-17 | -2.108e-17 2.108e-17
## 12.25 -1.306e-17 1.306e-17 | -1.306e-17 1.306e-17
## 13.25  1.464e-16 1.464e-16 |  1.464e-16 1.464e-16
## 14.25 -8.908e-17 8.908e-17 | -8.908e-17 8.908e-17
## 15.25 -5.852e-18 5.852e-18 | -5.852e-18 5.852e-18
## 16.25 -3.029e-17 3.029e-17 | -3.029e-17 3.029e-17
## 17.25  1.900e-16 1.900e-16 |  1.900e-16 1.900e-16
## 18.25 -1.511e-17 1.511e-17 | -1.511e-17 1.511e-17
## 19.25 -5.779e-17 5.779e-17 | -5.779e-17 5.779e-17
## 20.25  1.848e-16 1.848e-16 |  1.848e-16 1.848e-16

if(getRversion() >= "4.4.0") # *not* true for  R <= 4.3.3 :
    stopifnot(abs(relE) < 4e-16) # seen max = 1.900e-16

cat('Time elapsed: ', proc.time(),'\n') # "stats"
if(!interactive()) warnings()
