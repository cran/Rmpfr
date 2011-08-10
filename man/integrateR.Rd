\name{integrateR}
\title{One-Dimensional Numerical Integration - in pure R}
\alias{integrateR}
\alias{print.integrate}
\alias{show,integrate-method}
\description{
  Numerical integration of one-dimensional functions in pure \R,
  with care so it also works for \code{"mpfr"}-numbers.

  Currently, only classical Romberg integration of order \code{ord} is
  available.
}
\usage{
integrateR(f, lower, upper, \dots, ord = NULL,
           rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
           verbose = FALSE)
}
\arguments{
  \item{f}{an \R function taking a numeric or \code{"mpfr"} first
    argument and returning a numeric (or \code{"mpfr"}) vector of the
    same length.  Returning a non-finite element will generate an error.
  }
  \item{lower, upper}{the limits of integration.  Currently \emph{must}
    be finite.}
  \item{\dots}{additional arguments to be passed to \code{f}.}
  \item{ord}{integer, the order of Romberg integration to be used.  If
  this is \code{NULL}, as per default, the order is increased until
  convergence, see \code{rel.tol} and \code{abs.tol}.}
  \item{rel.tol}{relative accuracy requested.  The default is 1.2e-4.}
  \item{abs.tol}{absolute accuracy requested.}
  \item{verbose}{logical or integer, indicating if and how much
    information should be printed during computation.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  \code{rel.tol} cannot be less than \code{max(50*.Machine$double.eps,
    0.5e-28)} if \code{abs.tol <= 0}.
}
\note{
  \code{f} must accept a vector of inputs and produce a vector of function
  evaluations at those points.  The \code{\link{Vectorize}} function
  may be helpful to convert \code{f} to this form.

  Note that the default tolerances (\code{rel.tol}, \code{abs.tol}) are
  not very accurate, but the same as for \code{\link{integrate}}, which
  however often returns considerably more accurate results than
  requested. This is typically \emph{not} the case for
  \code{integrateR()}.
}
\value{
  A list of class \code{"integrate"} with components
  \item{value}{the final estimate of the integral.}
  \item{abs.error}{estimate of the modulus of the absolute error.}
  \item{subdivisions}{for Romberg, the number of function evaluations.}
  \item{message}{\code{"OK"} or a character string giving the error message.}
  \item{call}{the matched call.}
}
\references{
  Bauer, F.L. (1961)
  Algorithm 60 -- Romberg Integration,
  \emph{Communications of the ACM} \bold{4}(6), p.255.
}
\author{Martin Maechler}
\seealso{
  \R's standard, \code{\link{integrate}}, is much more adaptive,
  also allowing infinite integration boundaries, and typically
  considerably faster for a given accuracy.

  We use the same (actually a copy) \code{print} S3 method,
  \code{\link{print.integrate}()}, as provided by \R.% 'stats' package
}
\examples{

## See more  from  ?integrate
## this is in the region where  integrate() can get problems:
integrateR(dnorm,0,2000)
integrateR(dnorm,0,2000, rel.tol=1e-15)
integrateR(dnorm,0,2000, rel.tol=1e-15, verbose=TRUE)

## Demonstrating that 'subdivisions' is correct:
Exp <- function(x) { .N <<- .N+ length(x); exp(x) }
.N <- 0; str(integrateR(Exp, 0,1, rel.tol=1e-10), digits=15); .N

### Using high-precision functions -----

## Polynomials are very nice:
integrateR(function(x) (x-2)^4 - 3*(x-3)^2, 0, 5, verbose=TRUE)
# n= 1, 2^n=        2 | I =            46.04, abs.err =      98.9583
# n= 2, 2^n=        4 | I =               20, abs.err =      26.0417
# n= 3, 2^n=        8 | I =               20, abs.err =  7.10543e-15
## 20 with absolute error < 7.1e-15
I <- integrateR(function(x) (x-2)^4 - 3*(x-3)^2, 0, mpfr(5,128),
                rel.tol = 1e-20, verbose=TRUE)
I ; I$value  ## all fine

## with floats:
integrateR(exp,      0     , 1, rel.tol=1e-15, verbose=TRUE)
## with "mpfr":
(I <- integrateR(exp, mpfr(0,200), 1, rel.tol=1e-25, verbose=TRUE))
(I.true <- exp(mpfr(1, 200)) - 1)
## true absolute error:
stopifnot(print(as.numeric(I.true - I$value)) < 4e-25)
}
\keyword{math}
\keyword{utilities}

