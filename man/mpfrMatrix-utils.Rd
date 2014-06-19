\name{mpfrMatrix-utils}
%\alias{det}% "the function" (our copy of base::det)
\alias{determinant.mpfrMatrix}
%
\title{Functions for mpfrMatrix Objects}
\description{
  \code{determinant(x, ..)} computes the determinant of the mpfr square
  matrix \code{x}.  May work via coercion to \code{"numeric"}, i.e., compute
  \code{determinant(\link{asNumeric}(x), logarithm)}, if
  \code{asNumeric} is true, by default, if the dimension is larger
  than three.  Otherwise, use precision \code{precBits} for the
  \dQuote{accumulator} of the result, and use the
  recursive mathematical definition of the determinant (with
  computational complexity \eqn{n!}, where \eqn{n} is the matrix
  dimension, i.e., \bold{very} inefficient for all but small matrices!)
}
\usage{
\S3method{determinant}{mpfrMatrix}(x, logarithm = TRUE,
            asNumeric = (d[1] > 3), precBits = max(.getPrec(x)), \dots)
}
\arguments{
  \item{x}{an \code{\linkS4class{mpfrMatrix}} object of \emph{square} dimension.}
  \item{logarithm}{logical indicating if the \code{\link{log}} of the
    absolute determinant should be returned.
  }
  \item{asNumeric}{logical .. .. if rather
    \code{determinant(\link{asNumeric}(x), ...)} should be computed.
  }
  \item{precBits}{the number of binary digits for the result (and the
    intermediate accumulations).}
  \item{\dots}{unused (potentially further arguments passed to methods).}
}
\value{
  as \code{\link{determinant}()}, an object of S3 class \code{"det"}, a
  \code{\link{list}} with components
  \item{modulus}{the (logarithm of) the absolute value
    (\code{\link{abs}}) of the determinant of \code{x}.}
  \item{sign}{the sign of the determinant.}
}
%% \details{
%% }
%% \references{
%% }
\seealso{
  \code{\link{determinant}} in base \R, which relies on a fast LU decomposition.
  \code{\linkS4class{mpfrMatrix}}
}
\author{Martin Maechler}
\examples{
m6 <- mpfrArray(1:6, prec=128, dim = c(2L, 3L))
m6
S2 <- m6[,-3] # 2 x 2
S3 <- rbind(m6, c(1:2,10))
det(S2)
str(determinant(S2))
det(S3)
stopifnot(all.equal(det(S2), det(asNumeric(S2)), tolerance=1e-15),
	  all.equal(det(S3), det(asNumeric(S3)), tolerance=1e-15))
}
\keyword{array}
