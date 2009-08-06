\name{bessel_n}
\title{Bessel functions of Integer Order in multiple precisions}
\alias{j0}
\alias{j1}
\alias{jn}
\alias{y0}
\alias{y1}
\alias{yn}
\description{
  Bessel functions of integer orders, provided via arbitrary precision
  algorithms from the MPFR library.
}
\usage{
j0(x)
j1(x)
jn(n, x)
y0(x)
y1(x)
yn(n, x)
}
\arguments{
  \item{x}{a \code{\link{numeric}} or \code{\linkS4class{mpfr}} vector.}
  \item{n}{non-negative integer (vector).}
}
\value{
  Computes multiple precision versions of the Bessel functions of
  \emph{integer} order, \eqn{J_n(x)}{J[n](x)} and \eqn{Y_n(x)}{Y[n](x)}.
}
\seealso{\code{\link{besselJ}}, and \code{\link{besselY}} compute the
  same bessel functions but for arbitrary \emph{real} order and only
  precision of a bit more than ten digits.
}
\examples{
x <- (0:100)/8 # (have exact binary representation)
stopifnot(  all.equal(besselY(x, 0), bY0 <- y0(x))
          , all.equal(besselJ(x, 1), bJ1 <- j1(x))
          , all.equal(yn(0,x), bY0)
          , all.equal(jn(1,x), bJ1)
         )
}
\keyword{math}
