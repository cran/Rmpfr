\name{mpfr.utils}
\title{MPFR Number Utilities}
\alias{c.mpfr}
\alias{mpfr.is.0}
\alias{mpfr.is.integer}
\alias{mpfrVersion}
\alias{atan2,mpfr,ANY-method}
\alias{atan2,ANY,mpfr-method}
\alias{atan2,ANY,mpfrArray-method}
\alias{atan2,mpfrArray,ANY-method}
\alias{atan2,mpfrArray,mpfrArray-method}
\alias{atan2,mpfr,mpfr-method}
\alias{hypot}
%
\description{
  \code{mpfrVersion()} returns the version of the MPFR library which
  \pkg{Rmpfr} is currently linked to.

  \code{\link{c}(x,y,...)} can be used to combine MPFR numbers in the
  same way as regular numbers \bold{IFF} the first argument \code{x} is
  of class \code{\linkS4class{mpfr}}.

  \code{mpfr.is.0(.)} uses the MPFR library in the documented way to
  check if (a vector of) MPFR numbers are zero.

  \code{mpfr.is.integer(x)} uses the MPFR library in the documented way to
  check if (a vector of) MPFR numbers is integer \emph{valued}.  This is
  equivalent to \code{x == round(x)}, but \emph{not} at all to
  \code{is.integer(as(x, "numeric"))}.\cr
  You should typically rather use \code{\link{is.whole}(x)} instead.

  \code{hypot(x,y)} computes the hypothenuse length \eqn{z} in a rectangular
  triangle with \dQuote{leg} side lengths \eqn{x} and \eqn{y}, i.e.,
  \deqn{z = hypot(x,y) = \sqrt{x^2 + y^2},}{z = hypot(x,y) = sqrt(x^2 + y^2),}
  in a numerically stable way.
}
\usage{
mpfrVersion()
mpfr.is.0(x)
mpfr.is.integer(x)
\method{c}{mpfr}(...)

hypot(x,y)
}
\arguments{
  \item{x,y}{an object of class \code{\linkS4class{mpfr}}.}
  \item{\dots}{further \code{\linkS4class{mpfr}} class objects or simple
  numbers (\code{\link{numeric}} vectors) which are coerced to
  \code{mpfr} with default precision of 128 bits.}
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\section{Methods}{
  \describe{
    \item{atan2}{\code{signature(y = "mpfr", x = "ANY")}, and}
    \item{atan2}{\code{signature(x = "ANY", y = "mpfr")}: compute the
     arc-tangent of two arguments: \code{atan2(y, x)} returns the angle
     between the x-axis and the vector from the origin to \eqn{(x, y)},
     i.e., for positive arguments \code{atan2(y, x) == \link{atan}(y/x)}.}
  }
}
\value{
  \code{mpfr.is.0} returns a logical vector of length \code{length(x)}
  with values \code{TRUE} iff the corresponding \code{x[i]} is an MPFR
  representation of zero (\code{0}).
  \cr Similarly, \code{mpfr.is.integer} returns a logical vector of
  length \code{length(x)}.

  \code{mpfrVersion} returns an object of S3 class
  \code{"\link{numeric_version}"}, so it can be used in comparisons.

  The other functions return MPFR number (vectors), i.e., extending
  class \code{\linkS4class{mpfr}}.
}
\seealso{
  \code{\link{erf}} for special mathematical functions on MPFR;
  The class description \code{\linkS4class{mpfr}} mentioning the
  generic arithmetic and mathematical functions for which \code{"mpfr"}
  methods are available.
}
\examples{
mpfrVersion()

(x <- c(Const("pi", 64), mpfr(-2:2, 64)))
mpfr.is.0(x)        # one of them is
str(x[mpfr.is.0(x)])# but it may not have been obvious..

xy <- expand.grid(x = -2:2, y = -2:2) ; x <- xy[,"x"] ; y <- xy[,"y"]
a2. <- atan2(y,x)
stopifnot(all.equal(a2., atan2(as(y,"mpfr"), x)),
          mpfr.is.integer(mpfr(2, 500) ^ (1:200)))
}
\keyword{arith}
