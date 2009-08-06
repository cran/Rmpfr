\name{mpfr-utils}
\Rdversion{1.1}
\alias{getPrec}
\alias{mpfr_default_prec}
\alias{print.mpfrArray}
\alias{toNum}
\title{
  Rmpfr -- Utilities for Precision Setting, Printing, etc
}
\description{
  This page documents utilities from package \pkg{Rmpfr} which are
  typically not called by the user.  In some case the may come handy.
}
\usage{
getPrec(x)
mpfr_default_prec(prec)
\method{print}{mpfrArray}(x, digits = NULL, ...)
toNum(from)
}
\arguments{
  \item{x, from}{an \R object of class \code{"\linkS4class{mpfr}"}, or
    \code{"\linkS4class{mpfrArray}"}, respectively.}
  \item{prec}{a positive integer, or missing.}
  \item{digits, \dots}{further arguments to print methods.}
}
\details{
 ..
}
\value{
  \code{getPrec(x)} returns a \code{\link{integer}} vector of the same
  length as \code{x},

  \code{mpfr_default_prec()} returns the current MPFR default precision,
  an \code{\link{integer}}.  This is currently % ??
  not made use of, in all of package \pkg{Rmpfr}, where functions have
  their own default precision where needed.
  \cr
  \code{mpfr_default_prec(prec)} \emph{sets} the current MPFR default
  precision and returns the previous one; see above.

  \code{toNum(m)} returns a numeric \code{\link{array}} or
  \code{\link{matrix}}, when \code{m} is of class
  \code{"\linkS4class{mpfrArray}"} or \code{"\linkS4class{mpfrMatrix}"},
  respectively.  It should be equivalent to \code{as(m, "array")} or
  \code{... "matrix"}.
}
\seealso{
  Start using \code{\link{mpfr}(..)}, and compute with these numbers.
}
\examples{
getPrec(as(c(1,pi), "mpfr")) # 128 for both

(opr <- mpfr_default_prec()) ## typically  53, the MPFR system default
stopifnot(opr == (oprec <- mpfr_default_prec(70)),
          70  == mpfr_default_prec())
## and reset it:
mpfr_default_prec(opr)
}
\keyword{utilities}
