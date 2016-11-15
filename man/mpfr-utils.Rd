\name{mpfr-utils}
\title{Rmpfr -- Utilities for Precision Setting, Printing, etc}
\alias{getPrec}
\alias{getD}
\alias{mpfr_default_prec}
\alias{mpfr2array}
\alias{mpfrImport}
\alias{mpfrXport}
\alias{.mpfr2list}% for now only mentioned in example
%\alias{.mpfr1tolist}% not exported
\alias{print.mpfr}
\alias{print.mpfrArray}
\alias{toNum}
\description{
  This page documents utilities from package \pkg{Rmpfr} which are
  typically not called by the user, but may come handy in some
  situations.
}
\usage{
getPrec(x, base = 10, doNumeric = TRUE, is.mpfr = NA, bigq. = 128L)
getD(x)
mpfr_default_prec(prec)
\S3method{print}{mpfrArray}(x, digits = NULL, drop0trailing = FALSE,
      right = TRUE, \dots)
\S3method{print}{mpfr}(x, digits = NULL, drop0trailing = TRUE,
      right = TRUE, \dots)
toNum(from, rnd.mode = c('N','D','U','Z','A'))
mpfr2array(x, dim, dimnames = NULL, check = FALSE)

mpfrXport(x)
mpfrImport(mxp)
}
\arguments{
  \item{x, from}{typically, an \R object of class \code{"\linkS4class{mpfr}"}, or
    \code{"\linkS4class{mpfrArray}"}, respectively.  For
    \code{getPrec()}, any number-like \R object, or \code{\link{NULL}}.}
  \item{base}{(only when \code{x} is \code{\link{character}}) the base
    with respect to which \code{x[i]} represent numbers; \code{base}
    \eqn{b} must fulfill \eqn{2 \le b \le 62}{2 <= b <= 62}.}
  \item{doNumeric}{logical indicating \code{\link{integer}} or
    \code{\link{double}} typed \code{x} should be accepted and a default
    precision be returned.  Should typically be kept at default \code{TRUE}.}
  \item{is.mpfr}{logical indicating if \code{\link{class}(x)} is already
    known to be \code{"mpfr"}; typically should be kept at default, \code{NA}.}
  \item{bigq.}{for \code{getPrec()}, the precision to use for a big
    rational (class \code{"bigq"}); if not specified gives warning when used.}
  \item{prec}{a positive integer, or missing.}
  \item{drop0trailing}{logical indicating if trailing \code{"0"}s should
    be omitted.}
  \item{right}{logical indicating \code{print()}ing should right
    justify the strings; see \code{\link{print.default}()} to which it
    is passed.}
  \item{digits, \dots}{further arguments to print methods.}
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see details of
    \code{\link{mpfr}}.}
  \item{dim, dimnames}{for \code{"\linkS4class{mpfrArray}"}
    construction.}
  \item{check}{logical indicating if the mpfrArray construction should
    happen with internal safety check.  Previously, the implicit default
    used to be true.}
  \item{mxp}{an \code{"mpfrXport"} object, as resulting from
    \code{mpfrXport()}.}
}
\details{
  The \code{print} method is currently built on the \code{\link{format}}
  method for class \code{\linkS4class{mpfr}}.  This, currently does
  \emph{not} format columns jointly which leads to suboptimally looking
  output.  There are plans to change this.
}
\value{
  \code{getPrec(x)} returns a \code{\link{integer}} vector of the same
  length as \code{x} when that is positive, whereas \code{getPrec(NULL)}
  returns \code{mpfr_default_prec()}, see below.
  % the current MPFR default precision, see   \code{mpfr_default_prec()}
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
  \code{... "matrix"}.   Note that the slightly more general
  \code{\link{asNumeric}()} is preferred now.

  \code{mpfr2array()} a slightly more flexible alternative to
  \code{dim(.) <- dd}.
}
\note{
  \code{mpfrXport()} and \code{mpfrImport()} are \bold{experimental} and
  used to explore reported platform incompatibilities of
  \code{\link{save}()}d and \code{\link{load}()}ed \code{"mpfr"}
  objects between Windows and non-Windows platforms.

  In other words, the format of the result of \code{mpfrXport()} and
  hence the \code{mxp} argument to \code{mpfrImport()} are considered
  internal, not part of the API and subject to change.
}
\seealso{
  Start using \code{\link{mpfr}(..)}, and compute with these numbers.

  \code{\link{mpfrArray}(x)} is for numeric (\dQuote{non-mpfr})
  \code{x}, whereas \code{mpfr2array(x)} is for \code{"mpfr"} classed
  \code{x}, only.
}
\examples{
getPrec(as(c(1,pi), "mpfr")) # 128 for both

(opr <- mpfr_default_prec()) ## typically  53, the MPFR system default
stopifnot(opr == (oprec <- mpfr_default_prec(70)),
          70  == mpfr_default_prec())
## and reset it:
mpfr_default_prec(opr)

## Explore behavior of rounding modes 'rnd.mode':
x <- mpfr(10,99)^512 # too large for regular (double prec. / numeric):
sapply(c("N", "D", "U", "Z", "A"), function(RM)
       sapply(list(-x,x), function(.) toNum(., RM)))
##    N             D              U              Z    A
## -Inf          -Inf -1.797693e+308 -1.797693e+308 -Inf
##  Inf 1.797693e+308            Inf  1.797693e+308  Inf

## Printing of "MPFR" matrices is less nice than R's usual matrix printing:
m <- outer(c(1, 3.14, -1024.5678), c(1, 1e-3, 10,100))
m[3,3] <- round(m[3,3])
m
mpfr(m, 50)

B6 <- mpfr2array(Bernoulli(1:6, 60), c(2,3),
                 dimnames = list(LETTERS[1:2], letters[1:3]))
B6

## Looking at internal representation [for power users only!]:

i8 <- mpfr(-2:5, 32)
x4 <- mpfr(c(NA, NaN, -Inf, Inf), 32)
## The output of the following depends on the GMP "numb" size
## (32 bit vs. 64 bit), and may be even more platform specifics:
str( .mpfr2list(i8) )
str( .mpfr2list(x4) )

str(xp4 <- mpfrXport(x4))
stopifnot(identical(x4, mpfrImport(mpfrXport(x4))),
          identical(i8, mpfrImport(mpfrXport(i8))))
if(FALSE) ## FIXME: not yet working:
  stopifnot(identical(B6, mpfrImport(mpfrXport(B6))))
}
\keyword{utilities}
