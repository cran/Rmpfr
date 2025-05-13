\name{mpfr-distr-etc}
\title{Distribution Functions with MPFR Arithmetic}
\alias{dnorm}
\alias{dbinom}
\alias{dnbinom}
\alias{dgamma}
\alias{dpois}
\alias{dt}
\alias{pgamma}
\alias{pnorm}
\alias{mpfr-distr}% <- alternative for \link
\usage{%				>>>>  ../R/special-fun.R  <<<<
dpois (x, lambda, log = FALSE, useLog = )
dbinom (x, size, prob,     log = FALSE, useLog = , warnLog = TRUE)
dnbinom(x, size, prob, mu, log = FALSE, useLog = any(x > 1e6))
dnorm (x, mean = 0, sd = 1, log = FALSE)
dgamma(x, shape, rate = 1, scale = 1/rate, log = FALSE)
dt (x, df, ncp, log = FALSE)

pgamma(q, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE,
       rnd.mode = c('N','D','U','Z','A'))
pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
}
\description{
  For some \R standard (probability) density, distribution or quantile
  functions, we provide MPFR versions.
}
\arguments{
  \item{x,q, lambda, size,prob, mu, mean,sd, shape,rate,scale, df,ncp}{%
    \code{\link{numeric}} or \code{\linkS4class{mpfr}} vectors.  % for the function call.
    All of these are \dQuote{recycled} to the length of the longest one.
    For their meaning/definition, see the corresponding standard R
    (\pkg{stats} package) function.}
  \item{log, log.p, lower.tail}{logical, see
    \code{\link[stats:Normal]{pnorm}},
    \code{\link[stats:Poisson]{dpois}}, etc.}
  \item{useLog}{\code{\link{logical}} with default depending on \code{x}
    etc, indicating if log-scale computation should be used even when
    \code{log = FALSE}, for performance or against overflow / underflow.}
  \item{warnLog}{\code{\link{logical}} indicating if the \dQuote{mismatch}
    \code{log = TRUE, useLog = FALSE} should be warned about.}% somewhat experimental
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see details of \code{\link{mpfr}}.}
}
\details{
  \code{pnorm()} is based on \code{\link{erf}()} and \code{erfc()} which
  have direct MPFR counter parts and are both reparametrizations
  of \code{pnorm}, \code{erf(x) = 2*pnorm(sqrt(2)*x)} and
  \code{erfc(x) = 2* pnorm(sqrt(2)*x, lower=FALSE)}.

  \code{pgamma(q, sh)} is based on our \code{\link{igamma}(sh, q)}, see the
  \sQuote{Warning} there!
}
\note{E.g., for \code{pnorm(*, log.p = TRUE)} to be useful, i.e., not to
  underflow or overflow, you may want to extend the exponential range of
  MPFR numbers, using \code{\link{.mpfr_erange_set}()}, see the examples.
}
\value{
  A vector of the same length as the longest of \code{x,q, ...},
  of class \code{\linkS4class{mpfr}} with the high accuracy results of
  the corresponding standard \R function.
}
\seealso{\code{\link[stats:Normal]{pnorm}},
  \code{\link[stats]{dt}},
  \code{\link[stats]{dbinom}},
  \code{\link[stats]{dnbinom}},
  \code{\link[stats]{dgamma}},
  \code{\link[stats]{dpois}} in standard package \pkg{stats}.

  \code{\link{pbetaI}(x, a,b)} is a \code{\link{mpfr}} version of
  \code{\link{pbeta}} only for \emph{integer} \code{a} and \code{b}.
}
\examples{
x <- 1400+ 0:10
print(dpois(x, 1000), digits =18) ## standard R's double precision
(px <- dpois(mpfr(x, 120), 1000))## more accuracy for the same
px. <- dpois(mpfr(x, 120), 1000, useLog=TRUE)# {failed in 0.8-8}
stopifnot(all.equal(px, px., tol = 1e-31))
dpois(0:5, mpfr(10000, 80)) ## very small exponents (underflowing in dbl.prec.)

print(dbinom(0:8, 8, pr = 4 / 5), digits=18)
      dbinom(0:8, 8, pr = 4/mpfr(5, 99)) -> dB; dB

print(dnorm(     -5:5), digits=18)
      dnorm(mpfr(-5:5, prec=99))

## For pnorm() in the extreme tails, need an exponent range
## larger than the (MPFR and Rmpfr) default:
(old_eranges <- .mpfr_erange()) # typically -/+ 2^30:
log2(abs(old_eranges))   # 30  30
.mpfr_erange_set(value = (1-2^-52)*.mpfr_erange(c("min.emin","max.emax")))
log2(abs(.mpfr_erange()))# 62  62  *if* setup -- 2023-01: *not* on Winbuilder, nor
## other Windows where long is 4 bytes (32 bit) and the erange typically cannot be extended.
tens <- mpfr(10^(4:7), 128)
pnorm(tens, lower.tail=FALSE, log.p=TRUE) # "works" (iff ...)
## "the" boundary:
pnorm(mpfr(- 38581.371, 128), log.p=TRUE) # still does not underflow {but *.372 does}
## -744261105.599283824811986753129188937418  (iff ...)
.mpfr_erange()*log(2) # the boundary
##          Emin          Emax
## -3.196577e+18  3.196577e+18 (iff ...)

## reset to previous
.mpfr_erange_set( , old_eranges)
pnorm(tens, lower.tail=FALSE, log.p=TRUE) # all but first underflow to -Inf

## dnbinom(x, size, ..)  for large (x, size): .. already after fixing R-devel dnbinom()
xx <- 6e307
sz <- 1e308
dnb <- curve(dnbinom(xx, sz, prob = x, log=TRUE), 0, 1,  n = 1024 + 1,
             xlab = quote(prob), main = sprintf("dnbinom(\%s, \%s, prob=prob, log=TRUE)", xx, sz),
             col = 2, lwd=2)
x <- dnb$x
dnbM <- dnbinom(mpfr(xx, 128), mpfr(sz, 128), prob = x, log=TRUE)
lines(x, asNumeric(dnbM), col = adjustcolor(4, 1/3), lwd=5)

## dnbinom(x, size, ..)  for large (x, size): ..
for(x.n in list(c(7e305, 1e306), c(7e306, 1e307), c(7e307, 1e308))) {
    xx <- x.n[[1]] ; sz <- x.n[[2]]               # ============ here, we saw big jumps
    dnb <- curve(dnbinom(xx, sz, prob = x, log=TRUE), 0, 1,  n = 1024 + 1,  col = 2, lwd = 2,
                xlab = quote(prob), main = sprintf("dnbinom(\%s, \%s, prob=prob, log=TRUE)", xx, sz))
    x <- dnb$x; mtext(sfsmisc::shortRversion(), adj=1, cex = 3/4)
    dnbM <- dnbinom(mpfr(xx, 128), mpfr(sz, 128), prob = x, log=TRUE)
    lines(x, asNumeric(dnbM), col = adjustcolor(4, 1/3), lwd=5)
    if(dev.interactive()) Sys.sleep(1.5)
}

## pgamma() {and when igamma() is available}:
x <- c(10^(-20:-1), .5, 1:20, 10^(2:20))
xM <- mpfr(x, precBits = 128)
## CAREFUL --- some of these take *infinite* time ...
## subset , as "... infinite time ..."
iOk <- 1e-3 <= abs(x) & abs(x) <= 100
## x = 1e-3 is where our pgamma() {from igamma()} becomes very inaccurate
xm <- xM <- xM[iOk]; x <- x[iOk]
 ## sh.v <- c(1e-100, 1e-20, 1e-10, .5, 1,2,5, 10^c(1:10, 100, 300))
 ## sh.v <- c(1e-100, 1e-11,  1e-4, .5, 1,2,5, 10^c(1:5, 10, 100)) # less extreme ..
sh.v <- c(1e-100, 1e-11, 1e-4, .5, 1,2,5, 10^c(1:5,7)) # much less extreme than above ..
FT <- c("F", "T") # for printing
for(scale in c(1/2, 2))
 for(sh in sh.v) {
   cat(sprintf("scale = \%4.3g, shape= \%9g: ", scale, sh))
   stim <- system.time(
      for(ltail in c(FALSE, TRUE))
          for(lg in c(FALSE,TRUE)) {
              ae <- all.equal(pgamma(xM, sh, scale=scale, lower.tail=ltail, log.p=lg),
                              pgamma(x , sh, scale=scale, lower.tail=ltail, log.p=lg))
              if(!isTRUE(ae))
                  cat(sprintf(" ltail=\%s, lg=\%s: NOT eq.: \%s", FT[1+ltail], FT[1+lg], ae))
          }
   )
   cat(" user.time: ", stim[["user.self"]], "\n")
 } # for (sh ..)
## scale = 0.5, shape= 1e-100:  user.time:  0.292 
## scale = 0.5, shape=  1e-11:  user.time:  0.081 
## scale = 0.5, shape= 0.0001:  user.time:  0.051 
## scale = 0.5, shape=    0.5:  user.time:  0.05 
## scale = 0.5, shape=      1:  user.time:  0.031 


## scale = 0.5, shape=      2:  user.time:  0.032 
## scale = 0.5, shape=      5:  user.time:  0.031 
## scale = 0.5, shape=     10:  user.time:  0.032 
## scale = 0.5, shape=    100:  ltail=T, lg=T: NOT eq.: Mean abs diff: Inf user.time:  0.029 
## scale = 0.5, shape=   1000:  ltail=T, lg=T: NOT eq.: Mean abs diff: Inf user.time:  0.02 
## scale = 0.5, shape=  10000:  ltail=T, lg=T: NOT eq.: Mean abs diff: Inf user.time:  0.019 
## scale = 0.5, shape= 100000:  ltail=T, lg=T: NOT eq.: Mean abs diff: Inf user.time:  0.022 
## scale = 0.5, shape=  1e+10:  ltail=F, lg=F: NOT eq.: Numeric: lengths (0, 24) differ 
##                              ltail=F, lg=T: NOT eq.: Mean absolute difference: Inf
##                              ltail=T, lg=F: NOT eq.: 'is.NA' ...: 0 in current 24 in target
##                              ltail=T, lg=T: NOT eq.: 'is.NA' ...: 0 in current 24 in target
##                              user.time:  0.021 
## scale = 0.5, shape= 1e+100: gamma_inc.c:290: MPFR assertion failed: 
##					!(__builtin_expect(!!((flags) & (2)), 0))
## On Windows (erange etc): alread shape = 1e10  leads to the above MPFR assertion fail !!
}
\keyword{distribution}
