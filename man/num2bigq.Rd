\name{num2bigq}
\alias{num2bigq}
\title{BigQ / BigRational Approximation of Numbers}
\description{
  \code{num2bigq(x)} searches for \dQuote{small} denominator
  \code{bigq} aka \sQuote{bigRational} approximations to numeric or
  \code{"mpfr"} \code{x}.

  It uses the same continued fraction approximation as package
  \CRANpkg{MASS}' \code{\link[MASS]{fractions}()}, but using big integer,
  rational and mpfr-arithmetic from packages \CRANpkg{Rmpfr} and \CRANpkg{gmp}.
}
\usage{
num2bigq(x, cycles = 50L, max.denominator = 2^25, verbose = FALSE)
}
\arguments{
  \item{x}{numeric or mpfr-number like}
  \item{cycles}{a positive integer, the maximal number of approximation
    cycles, or equivalently, continued fraction terms to be used.}
  \item{max.denominator}{an \emph{approximate} bound on the maximal
    denominator used in the approximation.
    If small, the algorithm may use less than \code{cycles} cycles.
  }
  \item{verbose}{a logical indicating if some intermediate results should
    be printed during the iterative approximation.}
}
%% \details{
%% }
\value{
  a big rational, i.e.,  \code{\link[gmp]{bigq}} (from \CRANpkg{gmp})
  vector of the same length as \code{x}.
}
%% \references{

%% }
\author{Bill Venables and Brian Ripley, for the algorithm in \code{\link[MASS]{fractions}()};
  Martin Maechler, for the port to use \pkg{Rmpfr} and \code{gmp} arithmetic.}
\seealso{
  \code{\link{.mpfr2bigq}()} seems similar but typically uses much larger
  denominators in order to get full accuracy.
}
\examples{
num2bigq(0.33333)

num2bigq(pi, max.denominator = 200) #  355/113
num2bigq(pi) # much larger
num2bigq(pi, cycles=10) # much larger
}
\keyword{arith}
\keyword{math}
