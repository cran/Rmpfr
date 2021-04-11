\name{mpfr-distr-etc}
\title{Distribution Functions etc (MPFR)}
\alias{pnorm}
\alias{dnorm}
\alias{dbinom}
\alias{dnbinom}
\alias{dgamma}
\alias{dpois}
\alias{mpfr-distr}% <- alternative for \link
\usage{%				>>>>  ../R/special-fun.R  <<<<
dpois (x, lambda, log = FALSE, useLog = )
dbinom (x, size, prob,     log = FALSE, useLog = )
dnbinom(x, size, prob, mu, log = FALSE, useLog = any(x > 1e6))
dnorm (x, mean = 0, sd = 1, log = FALSE)
dgamma(x, shape, rate = 1, scale = 1/rate, log = FALSE)

pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
}
\description{
  For some \R standard (probability) density, distribution or quantile
  functions, we provide MPFR versions.
}
\details{
  \code{pnorm()} is based on \code{\link{erf}()} and \code{erfc()} which
  have direct MPFR counter parts and are both reparametrizations
  of \code{pnorm}, \code{erf(x) = 2*pnorm(sqrt(2)*x)} and
  \code{erfc(x) = 2* pnorm(sqrt(2)*x, lower=FALSE)}.
}
\note{E.g., for \code{pnorm(*, log.p = TRUE)} to be useful, i.e., not to
  underflow or overflow, you may want to extend the exponential range of
  MPFR numbers, using \code{\link{.mpfr_erange_set}()}, see the examples.
}
\arguments{
  \item{x,q, lambda, size,prob, mu, mean,sd, shape,rate,scale}{%
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
}
\value{
  A vector of the same length as the longest of \code{x,q, ...},
  of class \code{\linkS4class{mpfr}} with the high accuracy results of
  the corresponding standard \R function.
}
\seealso{\code{\link[stats:Normal]{pnorm}},
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
dpois(mpfr(x, 120), 1000)## more accuracy for the same
dpois(0:5, mpfr(10000, 80)) ## very small exponents (underflowing in dbl.prec.)

print(dbinom(0:8, 8, pr = 4 / 5), digits=18)
      dbinom(0:8, 8, pr = 4/mpfr(5, 99)) -> dB; dB

print(dnorm(     -5:5), digits=18)
      dnorm(mpfr(-5:5, prec=99))

## For pnorm() in the extreme tails, need an exponent range
## larger than the (MPFR and Rmpfr) default:
(old_eranges <- .mpfr_erange()) # typically -/+ 2^30
.mpfr_erange_set(value = (1-2^-52)*.mpfr_erange(c("min.emin","max.emax")))
tens <- mpfr(10^(4:7), 128)
pnorm(tens, lower.tail=FALSE, log.p=TRUE) # "works"
## reset to previous
.mpfr_erange_set( , old_eranges)
pnorm(tens, lower.tail=FALSE, log.p=TRUE) # all but first underflow to -Inf
}
\keyword{distribution}
