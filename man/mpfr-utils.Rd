\name{mpfr-utils}
\title{Rmpfr -- Utilities for Precision Setting, Printing, etc}
\alias{getPrec}
\alias{.getPrec}
\alias{getD}
\alias{mpfr_default_prec}
\alias{mpfr2array}
\alias{mpfrImport}
\alias{mpfrXport}
%\alias{.mpfr1tolist}% not exported
\alias{print.mpfr}
\alias{print.mpfrArray}
\alias{toNum}
\alias{.mpfr2list}
%-
\alias{.mpfr_formatinfo}
\alias{.mpfr2exp}
%-
\alias{.mpfr_erange}
\alias{.mpfr_erange_set}
\alias{.mpfr_erange_kinds}
\alias{.mpfr_erange_is_int}
\alias{.mpfr_maxPrec}
\alias{.mpfr_minPrec}
\description{
  This page documents utilities from package \pkg{Rmpfr} which are
  typically not called by the user, but may come handy in some
  situations.

  Notably, the (base-2) maximal (and minimal) precision and the
  \dQuote{erange}, the range of possible (base-2) exponents of
  \code{\link{mpfr}}-numbers can be queried and partly extended.
}
\usage{
getPrec(x, base = 10, doNumeric = TRUE, is.mpfr = NA, bigq. = 128L)
.getPrec(x)
getD(x)
mpfr_default_prec(prec)
\S3method{print}{mpfrArray}(x, digits = NULL, drop0trailing = FALSE,
      right = TRUE,
      max.digits = getOption("Rmpfr.print.max.digits", 999L),
      exponent.plus = getOption("Rmpfr.print.exponent.plus", TRUE),
      \dots)
\S3method{print}{mpfr}(x, digits = NULL, drop0trailing = TRUE,
      right = TRUE,
      max.digits = getOption("Rmpfr.print.max.digits", 999L),
      exponent.plus = getOption("Rmpfr.print.exponent.plus", TRUE),
      \dots)
toNum(from, rnd.mode = c('N','D','U','Z','A'))
mpfr2array(x, dim, dimnames = NULL, check = FALSE)

.mpfr2list(x, names = FALSE)

mpfrXport(x, names = FALSE)
mpfrImport(mxp)

.mpfr_formatinfo(x)
.mpfr2exp(x)

.mpfr_erange(kind = c("Emin", "Emax"), names = TRUE)
.mpfr_erange_set(kind = c("Emin", "Emax"), value)
.mpfr_erange_kinds
.mpfr_erange_is_int()
.mpfr_maxPrec()
.mpfr_minPrec()
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
  \item{max.digits}{a number (possibly \code{Inf}) to limit the number
    of (mantissa) digits to be printed, simply passed to
    \code{\link{formatMpfr}()}.  The default is finite to protect from
    printing very long strings which is often undesirable, notably in
    \code{\link{interactive}} use.}
  \item{exponent.plus}{logical, simply passed to \code{\link{formatMpfr}()}.
    Was \code{FALSE} hardwired in Rmpfr versions before 0.8-0, and hence is
    allowed to be tweaked by an \code{\link{options}()} setting.}
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see details of
    \code{\link{mpfr}}.}
  \item{dim, dimnames}{for \code{"\linkS4class{mpfrArray}"}
    construction.}
  \item{check}{logical indicating if the mpfrArray construction should
    happen with internal safety check.  Previously, the implicit default
    used to be true.}
  \item{names}{(for \code{.mpfr2list()}) \code{\link{logical}} or
    \code{\link{character}} vector, indicating if the list returned
    should have \code{\link{names}}.  If character, it specifies the
    names; if true, the names are set to \code{format(x)}.}
  \item{mxp}{an \code{"mpfrXport"} object, as resulting from
    \code{mpfrXport()}.}
  \item{kind}{a \code{\link{character}} string or vector, specifying the
    kind of \dQuote{erange} value; must be an element of
    \code{.mpfr_erange_kinds}, i.e., one of \code{"Emin"}, \code{"Emax"},
    \code{"min.emin"}, \code{"max.emin"}, \code{"min.emax"}, \code{"max.emax"}.}
  \item{value}{\code{\link{numeric}}, for \code{.mpfr_erange_set()} one
    number per \code{kind}.  Must be in range specified by the
    \code{*."emin"} and \code{*."emax"} erange values.}
}
\details{
  The \code{print} method is currently built on the \code{\link{format}}
  method for class \code{\linkS4class{mpfr}}.  This, currently does
  \emph{not} format columns jointly which leads to suboptimally looking
  output.  There are plans to change this.

  Note that \code{\link{formatMpfr}()} which is called by \code{print()}
  (or \code{show()} or \R's implicit printing) uses \code{max.digits =
    Inf}, differing from our \code{print()}'s default on purpose.
  If you do want to see the full accuracy even in cases it is large, use
  \code{\link{options}(Rmpfr.print.max.digits = Inf)} or
  \code{(.. = 1e7)}, say.

  The \code{.mpfr_erange*} functions (and variable) allow to query and set
  the allowed range of values for the base-2 \emph{exponents} of
  \code{"mpfr"} numbers.
  See the examples below and GNU MPFR library documentation on the C functions
  \code{mpfr_get_emin()}, \code{mpfr_set_emin(.)},
  \code{mpfr_get_emin_min()}, and \code{mpfr_get_emin_max()},
  (and those four with \sQuote{_emin} replaced by \sQuote{_emax} above).
}
\value{
  \code{getPrec(x)} returns a \code{\link{integer}} vector of length one or
  the same length as \code{x} when that is positive, whereas
  \code{getPrec(NULL)} returns \code{mpfr_default_prec()}, see below.
  % the current MPFR default precision, see   \code{mpfr_default_prec()}
  If you need to \emph{change} the precision of \code{x}, i.e., need
  something like \dQuote{setPrec}, use \code{\link{roundMpfr}()}.

  \code{.getPrec(x)} is a simplified version of \code{getPrec()} which only
  works for \code{"mpfr"} objects \code{x}.

  \code{getD(x)} is intended to be a fast version of \code{x@.Data},
  and should not be used outside of lower level functions.

  \code{mpfr_default_prec()} returns the current MPFR default precision,
  an \code{\link{integer}}.  This is currently % ??
  not made use of much in package \pkg{Rmpfr}, where functions have
  their own default precision where needed, and otherwise we'd rather not
  be dependent of such a \emph{global} setting.
  \cr
  \code{mpfr_default_prec(prec)} \emph{sets} the current MPFR default
  precision and returns the previous one; see above.

  \code{.mpfr_maxPrec()} and (less interestingly) \code{.mpfr_minPrec()}
  give the maximal and minimal base-2 precision allowed in the current
  version of the MPFR library linked to by \R package \pkg{Rmpfr}.
  The maximal precision is typically \eqn{2^{63}}, i.e.,
  \preformatted{ all.equal(.mpfr_maxPrec(), 2^63) }
  is typically true.

  \code{toNum(m)} returns a numeric \code{\link{array}} or
  \code{\link{matrix}}, when \code{m} is of class
  \code{"\linkS4class{mpfrArray}"} or \code{"\linkS4class{mpfrMatrix}"},
  respectively.  It should be equivalent to \code{as(m, "array")} or
  \code{... "matrix"}.   Note that the slightly more general
  \code{\link{asNumeric}()} is preferred now.

  \code{mpfr2array()} a slightly more flexible alternative to
  \code{dim(.) <- dd}.

  \code{.mpfr2exp(x)} returns the base-2 (integer valued) exponents of
  \code{x}, i.e., it is the \R interface to MPFR C's \code{mpfr_get_exp()}.
  The result is \code{\link{integer}}  iff \code{.mpfr_erange_is_int()}
  is true, otherwise \code{\link{double}}.  Note that the MPFR (4.0.1)
  manual says about \code{mpfr_get_exp()}: \emph{The behavior for NaN,
    infinity or zero is undefined}.

  \code{.mpfr_erange_is_int()} returns \code{TRUE} iff the
  \code{.mpfr_erange(c("Emin","Emax"))} range lies inside the range of \R's
  \code{\link{integer}} limits, i.e., has absolute values not larger than
  \code{\link{.Machine}$integer.max} (\eqn{ = 2^{31} - 1}).

  \code{.mpfr_formatinfo(x)} returns conceptually a subset of
  \code{\link{.mpfr2str}()}'s result, a list with three components
  \item{exp}{the base-2 exponents of \code{x}, identical to \code{.mpfr2exp(x)}.}
  \item{finite}{\code{\link{logical}} identical to \code{\link{is.finite}(x)}.}
  \item{is.0}{\code{\link{logical}} indicating if the corresponding
    \code{x[i]} is zero; identical to \code{\link{mpfrIs0}(x)}.}
  (Note that \code{\link{.mpfr2str}(x, .., base)$exp} is wrt \code{base}
  \emph{and} is not undefined but
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

## Ranges of (base 2) exponents of MPFR numbers:
.mpfr_erange() # the currently active range of possible base 2 exponents:

## A factory fresh setting fulfills
.mpfr_erange(c("Emin","Emax")) == c(-1,1) * (2^30 - 1)

## There are more 'kind's, the latter 4 showing how you could change the first two :
.mpfr_erange_kinds
.mpfr_erange(.mpfr_erange_kinds)
eLimits <- .mpfr_erange(c("min.emin", "max.emin", "min.emax", "max.emax"))
## Typically true in "current" MPFR versions:
eLimits == c(-1,1, -1,1) * (2^62 - 1)


## Looking at internal representation [for power users only!]:

i8 <- mpfr(-2:5, 32)
x4 <- mpfr(c(NA, NaN, -Inf, Inf), 32)
## The output of the following depends on the GMP "numb" size
## (32 bit vs. 64 bit), and may be even more platform specifics:
str( .mpfr2list(i8) )
str( .mpfr2list(x4, names = TRUE) )

str(xp4 <- mpfrXport(x4, names = TRUE))
stopifnot(identical(x4, mpfrImport(mpfrXport(x4))),
          identical(i8, mpfrImport(mpfrXport(i8))))
if(FALSE) ## FIXME: not yet working:
  stopifnot(identical(B6, mpfrImport(mpfrXport(B6))))
}
\keyword{utilities}
