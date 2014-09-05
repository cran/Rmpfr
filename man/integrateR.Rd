\name{integrateR}
\title{One-Dimensional Numerical Integration - in pure R}
\alias{integrateR}
\alias{print.integrateR}
\alias{show,integrateR-method}
\description{
  Numerical integration of one-dimensional functions in pure \R,
  with care so it also works for \code{"mpfr"}-numbers.

  Currently, only classical Romberg integration of order \code{ord} is
  available.
}
\usage{
integrateR(f, lower, upper, \dots, ord = NULL,
           rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
           max.ord = 19, verbose = FALSE)
}
\arguments{
  \item{f}{an \R function taking a numeric or \code{"mpfr"} first
    argument and returning a numeric (or \code{"mpfr"}) vector of the
    same length.  Returning a non-finite element will generate an error.
  }
  \item{lower, upper}{the limits of integration.  Currently \emph{must}
    be finite.  Do use \code{"mpfr"}-numbers to get higher than double
    precision, see the examples.}
  \item{\dots}{additional arguments to be passed to \code{f}.}
  \item{ord}{integer, the order of Romberg integration to be used.  If
    this is \code{NULL}, as per default, and either \code{rel.tol} or
    \code{abs.tol} are specified, the order is increased until
    convergence.}
  \item{rel.tol}{relative accuracy requested.  The default is 1.2e-4,
    about 4 digits only, see the Note.}
  \item{abs.tol}{absolute accuracy requested.}
  \item{max.ord}{only used, when neither \code{ord} or one of
    \code{rel.tol}, \code{abs.tol} are specified:  Stop Romberg
    iterations after the order reaches \code{max.ord}; may prevent
    infinite loops or memory explosion.}
  \item{verbose}{logical or integer, indicating if and how much
    information should be printed during computation.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  For convergence, \emph{both} relative and absolute changes must be
  smaller than \code{rel.tol} and \code{abs.tol}, respectively.

  \code{rel.tol} cannot be less than \code{max(50*.Machine$double.eps,
    0.5e-28)} if \code{abs.tol <= 0}.
}
\note{
  \code{f} must accept a vector of inputs and produce a vector of function
  evaluations at those points.  The \code{\link{Vectorize}} function
  may be helpful to convert \code{f} to this form.

  If you want to use higher accuracy, you \emph{must} set \code{lower} or
  \code{upper} to \code{"\link{mpfr}"} numbers (and typically lower the
  relative tolerance, \code{rel.tol}), see also the examples.

  Note that the default tolerances (\code{rel.tol}, \code{abs.tol}) are
  not very accurate, but the same as for \code{\link{integrate}}, which
  however often returns considerably more accurate results than
  requested.  This is typically \emph{not} the case for
  \code{integrateR()}.
}
\value{
  A list of class \code{"integrateR"} (as from standard \R's
  \code{\link{integrate}()}) with a \code{\link{print}} method and components
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
}
\note{
  We use practically the same \code{print} S3 method as
  \code{\link{print.integrate}}, provided by \R,% 'stats' package
  with a difference when the \code{message} component is not \code{"Ok"}.
}
\examples{

## See more  from  ?integrate
## this is in the region where  integrate() can get problems:
integrateR(dnorm,0,2000)
integrateR(dnorm,0,2000, rel.tol=1e-15)
(Id <- integrateR(dnorm,0,2000, rel.tol=1e-15, verbose=TRUE))
Id$value == 0.5 # exactly

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
## Now, using higher accuracy:
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

## Want absolute tolerance check only (=> set 'rel.tol' very high, e.g. 1):
(Ia <- integrateR(exp, mpfr(0,200), 1, abs.tol = 1e-6, rel.tol=1, verbose=TRUE))

## Set 'ord' (but no  '*.tol') --> Using 'ord'; no convergence checking
(I <- integrateR(exp, mpfr(0,200), 1,  ord = 13, verbose=TRUE))

}
\keyword{math}
\keyword{utilities}
