stopifnot(require("Rmpfr"))
(doExtras <- Rmpfr:::doExtras())
options(nwarnings = 50000, width = 99)

(do.pdf <- !dev.interactive(orNone = TRUE))
if(do.pdf) {
    pdf.options(width = 8.5, height = 6) # for all pdf plots
    pdf("special-fun.pdf")
}


## to enhance  |rel.Err| plots:  {also in ~/R/Pkgs/DPQ/tests/pow-tst.R }
drawEps.h <- function(p2 = -(53:51), side = 4, lty=3, lwd=2, col=adjustcolor(2, 1/2)) {
    abline(h = 2^p2, lty=lty, lwd=lwd, col=col)
    axis(side, las=2, line=-1, at = 2^p2,
         labels = as.expression(lapply(p2, function(p) substitute(2^E, list(E=p)))),
         col.axis = col, col=NA, col.ticks=NA)
}
mtextVersion <- function(adj = 1, col = 1) {
    mtext(osVersion, line=1, col=col, adj=adj)
    mtext(sfsmisc::shortRversion(spaces=FALSE), col=col, adj=adj)
}

all.eq.finite <- function(x,y, ...) {
    ## x = 'target'   y = 'current'
    if(any(is.finite(y[!(fx <- is.finite(x))])))
	return("current has finite values where target has not")
    if(any(is.finite(x[!(fy <- is.finite(y))])))
	return("target has finite values where current has not")
    ## now they have finite values at the same locations
    all.equal(x[fx], y[fy], ...)
}



n <- 1000
head(x <- mpfr(0:n, 100) / n)

stopifnot(exprs = {
    range(x) == 0:1
    all.equal(as.numeric(j0(x)), besselJ(as.numeric(x), 0), tol = 1e-14)
    all.equal(as.numeric(j1(x)), besselJ(as.numeric(x), 1), tol = 1e-14)
    all.equal(as.numeric(y0(x)), besselY(as.numeric(x), 0), tol = 1e-14)
    all.equal(as.numeric(y1(x)), besselY(as.numeric(x), 1), tol = 1e-14)
})

### pnorm() -> erf() : ----------------------------------------------------------
u <- 7*x - 2
stopifnot(all.equal(pnorm(as.numeric(u)),
		    as.numeric(pnorm(u)), tol = 1e-14))
## systematic random input testing:
set.seed(101)
if(doExtras) {
    nSim <- 50
    n2 <- 100
} else {
    nSim <- 10
    n2 <- 64
}
for(n in 1:nSim) {
    N <- rpois(1, lambda=n2)
    N3 <- N %/% 3
    x <- c(rnorm(N-N3), 10*rt(N3, df=1.25))# <- some large values
    m <- rnorm(N, sd = 1/32)
    s <- rlnorm(N, sd = 1/8)
    cEps <- .Machine$double.eps
    for(LOG in c(TRUE,FALSE))
	for(L.T in c(TRUE,FALSE)) {
	    p. <- pnorm( x, m=m,sd=s, log.p=LOG, lower.tail=L.T)
	    stopifnot(all.equal(p., pnorm(mpfr(x, precBits= 48), m=m,sd=s,
                                          log.p=LOG, lower.tail=L.T),
				tol = 128 * cEps))
	    stopifnot(all.equal(p., pnorm(mpfr(x, precBits= 60), m=m,sd=s,
                                          log.p=LOG, lower.tail=L.T),
				tol = 2 * cEps))
	}
    cat(".")
};cat("\n")
proc.time()

## Jerry Lewis - Aug 2, 2019
## Contrast the results of pnorm with double and mpfr inputs
x <- c(1:9, 5*(2:9), 10*(5:20)) ; x <- c(-rev(x), 0, x)
pdL <- pnorm(x, log.p=TRUE)
pdU <- pnorm(x, log.p=TRUE, lower.tail=FALSE)
stopifnot(exprs = {
    !is.unsorted(x)
    35 %in% x
    x == -rev(x) # exactly
    pdL == rev(pdU) # even exactly, currently
})
mx <- mpfr(x, precBits = 128)
pmL <- pnorm(mx, log.p=TRUE)
pmU <- pnorm(mx, log.p=TRUE, lower.tail=FALSE)
stopifnot(exprs = {
    pmL < 0 # not true for 'pdL' which underflows
    pmL == rev(pmU) # even exactly, currently
    all.equal(pmL, pdL, tol=4e-16) # 'tol=0' shows 4.46e-17
})
## some explorations :
dlp <- diff(log(-pmL))/diff(x)
n <- length(x)
x.1 <- (x[-1] + x[-n])/2
plot(x.1, dlp, type="b", ylab = "d/dx  log(-pnorm(., log=TRUE))"); mtextVersion()
plot(x.1[-1], diff(dlp)/diff(x.1), type="b", ylab = "d^2/dx^2  log(-pnorm(., log=TRUE))")
stopifnot(exprs = {
    -1 < (d2 <- diff(dlp)/diff(x.1))
    d2 < 0
    diff(d2) < 0
})
x.3 <- x.1[-c(1L,n-1L)]
plot(x.3, -diff(d2)/ diff(x.1)[-1], type="o", log="y")




### Riemann's Zeta function: ----------------------------------------------------

## -- integer arguments --
stopifnot(all(mpfrIs0(zeta(-2*(1:100)))))

k.neg <- 2*(-100:0) - 1
Z.neg <- zeta(k.neg)
plot(k.neg, abs(as.numeric(Z.neg)), type = "l", log="y")

Pi <- Const("pi", 128L)

## confirm published value of Euler's gamma to 100 digits
pub.g <-
    paste("0.5772156649", "0153286060", "6512090082", "4024310421", "5933593992",
	  "3598805767", "2348848677", "2677766467", "0936947063", "2917467495",
	  sep="")

## almost =
our.g <- Const("gamma", log2(10) * 100) # 100 digits
(ff.g <- .mpfr2str(our.g))


M <- function(x) mpfr(x, 128L)
stopifnot(all.equal(zeta( 0), -1/2,      tol = 2^-100)
	  , all.equal(zeta(-1), -1/M(12),  tol = 2^-100)
	  , all.equal(zeta(-3),  1/M(120), tol = 2^-100)
	  ## positive ones :
	  , all.equal(zeta(2),  Pi^2/6,   tol = 2^-100)
	  , all.equal(zeta(4),  Pi^4/90,  tol = 2^-100)
	  , all.equal(zeta(6),  Pi^6/945, tol = 2^-100)
	  )

### Exponential Integral Ei(.)
curve(Ei, 0,5, n=5001)
if(mpfrVersion() >= "3") { ## only available since MPFR 3.0.0
  ### Airy function Ai(.)
  curve(Ai, -10, 5, n=5001); abline(h=0,v=0, col="gray", lty=3)
}

### Utilities  hypot(), atan2() : --------------------------------------------------------------

## ======= TODO! ========

## beta(), lbeta()
## ---------------
## The simplistic "slow" versions:
B  <- function(a,b) { a <- as(a, "mpfr"); b <- as(b, "mpfr"); gamma(a)*gamma(b) / gamma(a+b) }
lB <- function(a,b) { a <- as(a, "mpfr"); b <- as(b, "mpfr"); lgamma(a)+lgamma(b) - lgamma(a+b) }

## For partly *integer* arguments
Bi1 <- function(a,b) 1/(a*chooseMpfr(a+b-1, a)) # a must be integer >= 0
Bi2 <- function(a,b) 1/(b*chooseMpfr(a+b-1, b)) # b must be integer >= 0

x <- 1:10 + 0 ; (b10 <- mpfr(x, 128L))

stopifnot(all.equal(	B(1,b10),  1/x),
	  all.equal(	B(2,b10),  1/(x*(x+1))),
	  all.equal( beta(1,b10),  1/x),
	  all.equal( beta(2,b10),  1/(x*(x+1))),
	  TRUE)

if(do.pdf) { dev.off(); pdf("special-fun-beta.pdf") }


x <- -10:10 + 0; X <- mpfr(x, 128L)
stopifnot(exprs = {
    Bi1(1,X) == (B1x <- Bi2(X,1))
    Bi1(2,X) == (B2x <- Bi2(X,2))
    Bi1(3,X) == (B3x <- Bi2(X,3))
    all.equal(B1x,  1/x,               tol= 4e-16)
    all.equal(B2x,  1/(x*(x+1)),       tol= 8e-16)
    all.equal(B3x,  2/(x*(x+1)*(x+2)), tol=16e-16)
    ## these the "poles" are all odd i.e. result in { +Inf / -Inf / NaN}
    ## are all "ok" {e.g. 1/(x*(x+1)) gives (-Inf, Inf) for x = -1:0 }
    all.eq.finite(beta(1,X),  1/x)
    all.eq.finite(beta(X,2),  1/(x*(x+1)))
    all.eq.finite(beta(3,X),  2/(x*(x+1)*(x+2)), tol=16e-16)
})

## (a,b)  *both* integer, one negative:
for(i in (-20):(-1)) {
    cat(i,":\n")
    a <- mpfr(i, 99)
    i1 <- i+1
    b. <- seq_len(-i1)
    Bab <- beta(a, b.)
    stopifnot(is.nan(beta(a, (i1:0))), is.nan(lbeta(a, (i1:0))),
	      all.equal(Bab, Bi2(a, b.),             tol=1e-20),
	      all.equal(lbeta(a, b.), log(abs(Bab)), tol=1e-20), allow.logical0 = TRUE)
}

## (a,b) all positive
c10 <- b10 + 0.25
for(a in c(0.1, 1, 1.5, 2, 20)) {
    stopifnot(all.equal( B(a,b10), (bb <-  beta(a, b10))),
	      all.equal(lB(a,b10), (lb <- lbeta(a, b10))), all.equal(lb, log(bb)),
	      all.equal( B(a,c10), (bb <-  beta(a, c10))),
	      all.equal(lB(a,c10), (lb <- lbeta(a, c10))), all.equal(lb, log(bb)),
	      TRUE)
}

## However, the speedup is *not* much (50%) when applied to vectors:
stopifnot(validObject(xx <- outer(b10, runif(20))),
	  dim(xx) == c(length(b10), 20),
	  validObject(vx <- as(xx, "mpfr")), class(vx) == "mpfr", is.null(dim(vx)))
C1 <- replicate(10, system.time(bb <<- beta(vx, vx+2)))
C2 <- replicate(10, system.time(b2 <<-    B(vx, vx+2)))
summary(1000*C1[1,]) ##  80.3 {cmath-5, 2009}
summary(1000*C2[1,]) ## 125.1 { " }
stopifnot(all.equal(bb, b2))
## and for a single number, the speedup is a factor 3:
x1 <- vx[1]; x2 <- x1+2
system.time(for(i in 1:100) bb <- beta(x1, x2))# .27
system.time(for(i in 1:100) b2 <-    B(x1, x2))# .83

## a+b is integer <= 0, but a and b are not integer:
a <- b <- .5 + -10:10
ab <- data.matrix(expand.grid(a=a, b=b, KEEP.OUT.ATTRS=FALSE))
ab <- mpfr(ab[rowSums(ab) <= 0, ], precBits = 128)
stopifnot( beta(ab[,"a"], ab[,"b"]) == 0,
	  lbeta(ab[,"a"], ab[,"b"]) == -Inf)
## was  NaN  in Rmpfr <= 0.5-2

stopifnot(all.equal(6 * beta(mpfr(1:3,99), -3.), c(-2,1,-2), tol=1e-20))
## add more checks, notably for b (> 0)  above and below the "large_b" in
## ../src/utils.c :
bb <- beta(mpfr(1:23, 128), -23)
stopifnot(all.equal(bb, Bi1(1:23, -23), tol=1e-7))
                                        # Bi1() does not get high prec for small b
## can be written via rationals:  N / D :
bn <- c(330, -360, 468, -728, 1365, -3120, 8840, -31824,
        151164, -1007760, 10581480, -232792560)
bn <- c(rev(bn[-1]), bn)
bd <- 24* as.bigz(2 * 3 * 5 * 7 * 11) * 13 * 17 * 19 * 23
stopifnot(all.equal(bb, as(bn/bd,"mpfr"), tol=0))

stopifnot(all.equal(6 * beta(mpfr(1:3,	99), -3.),
			     c(-2,1,-2),	    tol=1e-20),
	  all.equal(   lbeta(mpfr(1:3, 128), -3.),
		    log(mpfr(c( 2,1, 2), 128) / 6), tol=1e-20))

## add more checks, notably for b (> 0)  above and below the "large_b" in
## ../src/utils.c :
bb <- beta(mpfr(1:23, 128), -23)
stopifnot(all.equal(bb, Bi1(1:23, -23), tol=1e-7))
					# Bi1() does not get high prec for small b
## can be written via rationals:  N / D :
bn <- c(330, -360, 468, -728, 1365, -3120, 8840, -31824,
	151164, -1007760, 10581480, -232792560)
bn <- c(rev(bn[-1]), bn)
bd <- 24* as.bigz(2 * 3 * 5 * 7 * 11) * 13 * 17 * 19 * 23
stopifnot(all.equal(bb, as(bn/bd,"mpfr"), tol=0))

## 2) add check for 'b' >  maximal unsigned int {so C code uses different branch}
two <- mpfr(2, 128)
for(b in list(mpfr(9, 128), mpfr(5, 128)^10, two^25, two^26, two^100)) {
    a <- -(b+ (1:7))
    stopifnot(a+b == -(1:7), # just ensuring that there was no cancellation
	      is.finite( B <-  beta(a,b)), ## was NaN ..
	      is.finite(lB <- lbeta(a,b)), ## ditto
	      all.equal(log(abs(B)), lB),
	      TRUE)
}

ee <- c(10:145, 5*(30:59), 10*(30:39), 25*(16:30))
b <- mpfr(2, precBits = 10 + max(ee))^ee # enough precision {now "automatic"}
stopifnot((b+4)-b == 4, # <==> enough precision above
	  b == (b. <- as(as(b,"bigz"),"mpfr")))
(pp <- getPrec(b.))# shows why b. is not *identical* to b.
system.time(Bb <- beta(-b-4, b))# 0.334 sec
if(dev.interactive())
    plot(ee, asNumeric(log(Bb)), type="o",col=2)
lb <- asNumeric(log(Bb))
## using  coef(lm(lb ~ ee))
stopifnot(all.equal(lb, 3.175933 -3.46571851*ee, tol = 1e-5))# 4.254666 e-6


bb <- beta(           1:4,   mpfr(2,99))
stopifnot(identical(bb, beta(mpfr(2,99), 1:4)),
	  all.equal((2*bb)*cumsum(1:4), rep(1, 4), tol=1e-20),
	  getPrec(bb) == 128)


##-- The d*() density functions from ../R/special-fun.R  |  ../man/distr-etc.Rd ---

if(do.pdf) { dev.off(); pdf("special-fun-density.pdf") }

dx <- 1400+ 0:10
mx <- mpfr(dx, 120)
nx <- sort(c(c(-32:32)/2, 50*(-8:8)))

xL <- 2^(989+(0:139)/4) # "close" to double.xmax
dnbD   <- dnbinom(xL, prob=1-1/4096, size=1e307, log=TRUE)# R's own
iF <- -(130:140) # index of finite dnbD[]
dnbx8  <- dnbinom(xL, prob=1-mpfr(2, 2^ 8)^-12, size=1e307, log=TRUE)
dnbx10 <- dnbinom(xL, prob=1-mpfr(2, 2^10)^-12, size=1e307, log=TRUE)
dnbx13 <- dnbinom(xL, prob=1-mpfr(2, 2^13)^-12, size=1e307, log=TRUE)

stopifnot(exprs = {
    all.equal(dpois(dx, 1000), dpois(mx, 1000), tol = 3e-13) # 64b Lnx: 7.369e-14
    all.equal(dbinom(0:16, 16, pr = 4 / 5),
              dbinom(0:16, 16, pr = 4/mpfr(5, 128)) -> db, tol = 5e-15)# 64b Lnx: 4.3e-16
    all.equal(dnorm(     -3:3,       m=10, s=1/4),
              dnorm(mpfr(-3:3, 128), m=10, s=1/4), tol = 1e-15) # 64b Lnx: 6.45e-17
    all.equal(dnorm(nx), dnorm(mpfr(nx, 99)), tol = 1e-15)
    all.equal(dnorm(     nx,      m = 4, s = 1/4),
              dnorm(mpfr(nx, 99), m = 4, s = 1/4), tol = 1e-15)
    all.equal(dnorm(     nx,      m = -10, s = 1/4, log=TRUE),
              dnorm(mpfr(nx, 99), m = -10, s = 1/4, log=TRUE), tol = 1e-15)
    ## t-distrib. :
    all.equal(dt(nx, df=3), dt(mpfr(nx, 99), df=3), tol = 1e-15)
    all.equal(dt(     nx,      df = 0.75),
              dt(mpfr(nx, 99), df = 0.75), tol = 1e-15)
    all.equal(dt(     nx,      df = 2.5, log=TRUE),
              dt(mpfr(nx, 99), df = 2.5, log=TRUE), tol = 1e-15)
    ## negative binomial  dnbinom():
    all.equal(dnbx13, dnbx10, tol = 2^-999) # see 2^-1007, but not 2^-1008
    all.equal(dnbx13, dnbx8,  tol = 2^-238) # see 2^-239,  but not 2^-240
    all.equal(dnbx10[iF], dnbD[iF], tol = 6e-16) # R's *is* accurate here (seen 2.9e-16)
})


## plot dt() "error" of R's implementation
nx <- seq(-100, 100, by=1/8)
dtd <- dt(     nx,        df= .75)
dtM <- dt(mpfr(nx,  256), df= .75)
if(doExtras) withAutoprint({
 system.time(
  dtMx <- dt(mpfr(nx, 2048), df= .75) ) # 2.5 sec
 stopifnot(all.equal(dtMx, dtM, tol = 2^-254)) # almost all of dtM's 256 bits are correct
})
relE <- asNumeric(dtd/dtM - 1)
plot(relE ~ nx,      type="l", col=2); mtextVersion()
plot(abs(relE) ~ nx, type="l", col=2, log="y", ylim=c(5e-17, 1.5e-15))

## ============== even smaller 'df' such that lgamma1p(df) is better than lgamma(1+df) ====

require(sfsmisc)# -> eaxis(); relErrV()

u <- sort(outer(10^-(20:1), c(1,2,5))) # *not* "exact" on purpose
## .. unfinished .. exploring *when* dt() would suffer from inaccurate stirlerr()  -- would it?

nu <- 2^-(70:1)
dt10  <- dt(     10,        df=nu)
dt10M <- dt(mpfr(10, 1024), df=nu)
re10 <- asNumeric(relErrV(dt10M, dt10))

plot(re10 ~ nu, type="l", lwd=2, log="x", main = quote(rel.Err( dt(10, df==nu) )),
     xaxt="n"); eaxis(1, nintLog=20)
mtextVersion()
abline(h = (-1:1)*2^-53, lty=4, col=adjustcolor("blue", 1/2))

plot(abs(re10) ~ nu, type="l", lwd=2, log="xy",
     xlab = quote(df == nu), ylab = quote(abs(relE)),
     main = quote(abs(rel.Err( dt(10, df==nu) ))), xaxt="n", yaxt="n")
eaxis(1, nintLog=20); eaxis(2); drawEps.h()

x0 <- c(0, 10^(-5:10)) # only >= 0 should be sufficient; x0 <- c(-rev(x0),0,x0)
stopifnot(!is.unsorted(nu), # just for plotting ..
          !is.unsorted(x0))
xnu <- expand.grid(x=x0, df=nu)
dt2  <- with(xnu, dt(     x,       df=df))
dtM2 <- with(xnu, dt(mpfr(x, 512), df=df))
str(relE2 <- `attributes<-`(asNumeric(relErrV(dtM2, dt2)),
                            attr(xnu, "out.attrs")))

## consistency check that with() etc was fine:
stopifnot(identical(re10, unname(relE2[which(x0 == 10), ])))

filled.contour(x=log10(1e-7+x0), y=log10(nu), z = relE2)
filled.contour(x=log10(1e-7+x0), y=log10(nu), z = abs(relE2))
## around nu = 10^-16 is the most critical place

(pch <- c(1L:9L, 0L, letters, LETTERS)[1:ncol(relE2)])

matplot(x0+1e-7, relE2, type="b", log="x", main="rel.err{  dt(x, df=df) }")
legend("topright", legend = paste0("df=",formatC(nu,wid=3)), ncol=7,
       bty="n", lwd=1, pch=pch, col=1:6, lty=1:5, cex = 0.8)
abline(h = c(-4:4)*2^-53, lty=3, col="gray")

matplot(nu, t(relE2), type="b", log="x", main="rel.err{  dt(x, df=df) }")
legend("topright", legend = paste0("x=",formatC(x0,wid=3)), ncol=7,
       bty="n", lwd=1, pch=pch, col=1:6, lty=1:5, cex = 0.8)
abline(h = c(-4:4)*2^-53, lty=3, col="gray")

matplot(nu, pmax(abs(t(relE2)), 1e-19), type="b", log="xy", axes=FALSE, ylab = quote(abs("rel Err")),
        ylim = c(7e-17, max(abs(relE2))), main="|rel.err{ dt(x, df=df)}|")
eaxis(1, nintLog=22) ; eaxis(2, line=-1/2); drawEps.h()
legend("topright", legend = paste0("x=",formatC(x0,wid=3)), ncol=7,
       bty="n", lwd=1, pch=pch, col=1:6, lty=1:5, cex = 0.8)


1
## dnbinom() -- has mode as expected, but with huge size, the scales are "off reality" ..

### ..... TODO !

### dgamma(): ----------------------------------------------------
if(do.pdf) { dev.off(); pdf("special-fun-dgamma.pdf") }

xe <- c(-2e5, -1e5, -2e4, -1e4, -2000, -1000, -500, -200, -100, -50, -20, -10)
(xe <- c(xe, -8:8, -rev(xe)))
two <- mpfr(2, 64)
## For centering at E[.], will use xP(x, shp) :
xP <- function(x, d) x - d*(x > d)
aEQformat <- function(xy, ...) format(xy, digits = 7, ...)
allEQ_0 <- function (target, current, ...)
    all.equal(target, current, tolerance = 0, formatFUN = aEQformat, ...)
stopIfNot <-
    if("allow.logical0" %in% names(formals(stopifnot))) { # experimental (MM only)
        stopifnot
    } else function(exprs, allow.logical0) stopifnot(exprs=exprs)

for(shp in c(2^c(-20, -3, -1:1, 4, 10, 50))) {
    cat("shape = 2^", log2(shp), ":\n-------------\n")
    d.dg  <- dgamma(xP(2 ^ xe, shp), shape=shp)
    m.dg  <- dgamma(xP(two^xe, shp), shape=shp)
    m.ldg <- dgamma(xP(two^xe, shp), shape=shp, log=TRUE)
    stopIfNot(exprs = {
        !is.unsorted(xe)
        is.finite(m.dg)
        m.dg >= 0
        shp > 1  || all(diff(m.dg) <= 0)
        shp > 100|| all((m.dg > 0) >= (d.dg > 0))
        any(fin.d <- is.finite(d.dg))
        m.dg[!fin.d] > 1e300
        { cat("all.EQ(<mpfr>, <doubl>):", allEQ_0(m.dg[fin.d], d.dg[fin.d]), "\n")
          shp > 100  ||                   all.equal(m.dg[fin.d], d.dg[fin.d],
                                                    tol = 1e-13) # 2.063241e-14
        }
        ## compare with log scale :
        if(any(pos.d <- m.dg > 0)) {
            cat("all.EQ(log(d), d*(log)):",
              allEQ_0  (log(m.dg[pos.d]), m.ldg[pos.d]),"\n")
              all.equal(log(m.dg[pos.d]), m.ldg[pos.d], tol = 1e-14)
        }
    }, allow.logical0 = TRUE)
}

cat('Time elapsed: ', proc.time(),'\n') # "stats"
if(!interactive()) warnings()
