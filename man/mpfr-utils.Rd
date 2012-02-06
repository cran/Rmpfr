\name{mpfr-utils}
\alias{getD}
\alias{getPrec}
\alias{mpfr_default_prec}
\alias{mpfr2array}
\alias{print.mpfr}
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
getPrec(x, base = 10, doNumeric = TRUE, is.mpfr = NA)
getD(x)
mpfr_default_prec(prec)
\S3method{print}{mpfrArray}(x, digits = NULL, drop0trailing = FALSE,
      right = TRUE, \dots)
\S3method{print}{mpfr}(x, digits = NULL, drop0trailing = TRUE,
      right = TRUE, \dots)
toNum(from)
mpfr2array(x, dim, dimnames = NULL, check = FALSE)
}
\arguments{
  \item{x, from}{an \R object of class \code{"\linkS4class{mpfr}"}, or
    \code{"\linkS4class{mpfrArray}"}, respectively.}
  \item{base}{(only when \code{x} is \code{\link{character}}) the base
    with respect to which \code{x[i]} represent numbers; \code{base}
    \eqn{b} must fulfill \eqn{2 \le b \le 36}{2 <= b <= 36}.}
  \item{doNumeric}{logical indicating \code{\link{integer}} or
    \code{\link{double}} typed \code{x} should be accepted and a default
    precision be returned.  Should typically be kept at default \code{TRUE}.}
  \item{is.mpfr}{logical indicating if \code{\link{class}(x)} is already
    known to be \code{"mpfr"}; typically should be kept at default, \code{NA}.}
  \item{prec}{a positive integer, or missing.}
  \item{drop0trailing}{logical indicating if trailing \code{"0"}s should
    be omitted.}
  \item{right}{logical indicating \code{print()}ing should right
    justify the strings; see \code{\link{print.default}()} to which it
    is passed.}
  \item{digits, \dots}{further arguments to print methods.}
  \item{dim, dimnames}{for \code{"\linkS4class{mpfrArray}"}
    construction.}
  \item{check}{logical indicating if the mpfrArray construction should
    happen with internal safety check.  Previously, the implicit default
    used to be true.}
}
\details{
  The \code{print} method is currently built on the \code{\link{format}}
  method for class \code{\linkS4class{mpfr}}.  This, currently does
  \emph{not} format columns jointly which leads to suboptimally looking
  output.  There are plans to change this.
}
\value{
  \code{getPrec(x)} returns a \code{\link{integer}} vector of the same
  length as \code{x}.
  If you need to \emph{change} the precision of \code{x}, i.e., need
  something like \dQuote{setPrec}, use \code{\link{roundMpfr}()}.

  \code{getD(x)} is intended to be a fast version of \code{x@.Data},
  and should not be used outside of lower level functions.

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

  \code{mpfr2array()} a slightly more flexible alternative to
  \code{dim(.) <- dd}.
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

## Printing of "MPFR" matrices is less nice than R's usual matrix printing:
m <- outer(c(1, 3.14, -1024.5678), c(1, 1e-3, 10,100))
m[3,3] <- round(m[3,3])
m
mpfr(m, 50)

mpfr2array(Bernoulli(1:6, 60), c(2,3),
           dimnames = list(LETTERS[1:2], letters[1:3]))
}
\keyword{utilities}
