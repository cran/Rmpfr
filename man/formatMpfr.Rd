\name{formatMpfr}
\Rdversion{1.1}
\alias{formatMpfr}
\alias{formatN.mpfr}
\title{Formatting MPFR (multiprecision) Numbers}
\description{
  Flexible formatting of \dQuote{multiprecision numbers}, i.e., objects
  of class \code{\linkS4class{mpfr}}.  \code{formatMpfr()} is also the
  \code{mpfr} method of the generic \code{\link{format}} function.

  The \code{formatN()} methods for \code{\linkS4class{mpfr}} numbers
  renders them differently than their double precision equivalents, by
  appending \code{"_M"}.
}
\usage{
formatMpfr(x, digits = NULL, trim = FALSE, scientific = NA,
	   showNeg0 = TRUE,
	   big.mark = "", big.interval = 3L,
	   small.mark = "", small.interval = 5L, decimal.mark = ".",
	   zero.print = NULL, drop0trailing = FALSE, ...)
\S3method{formatN}{mpfr}(x, drop0trailing = TRUE, \dots)
}
\arguments{
  \item{x}{an MPFR number (vector or array).}
  \item{digits}{how many significant digits are to be used for
    numeric and complex \code{x}.  The default, \code{NULL}, uses
    enough digits to represent the full precission.}
  \item{trim}{logical; if \code{FALSE}, numbers are right-justified to a
    common width: if \code{TRUE} the leading blanks for justification
    are suppressed.}
  \item{scientific}{either a logical specifying whether
    MPFR numbers should be encoded in scientific
    format, or an integer penalty (see \code{\link{options}("scipen")}).
    Missing values correspond to the current default penalty.}
  \item{showNeg0}{logical indicating if \dQuote{\bold{neg}ative} zeros
    should be shown with a \code{"-"}.  The default, \code{TRUE} is
    intentially different from \code{\link{format}(<numeric>)}.}
  \item{big.mark, big.interval, small.mark,
    small.interval, decimal.mark, zero.print, drop0trailing}{%
    used for prettying decimal sequences, these are passed to
    \code{\link{prettyNum}} and that help page explains the details.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  a character vector or array, say \code{cx}, of the same length as
  \code{x}.  Since Rmpfr version 0.5-3 (2013-09), if \code{x} is an
  \code{\linkS4class{mpfrArray}}, then \code{cx} is a character
  \code{\link{array}} with the same \code{\link{dim}} and
  \code{\link{dimnames}} as \code{x}.
}
\references{
  The MPFR manual's description of \samp{mpfr_get_str()} which is the
  C-internal workhorse for the (internal) \R function \code{.mpfr2str()}
  on which \code{formatMpfr} builds.
}
\author{Martin Maechler}
\seealso{
  \code{\link{mpfr}} for creation and
  the \code{\linkS4class{mpfr}} class description with its many methods.
  The \code{\link{format}} generic, and the \code{\link{prettyNum}}
  utility on which \code{formatMpfr} is based as well.
  The S3 generic function \code{\link[gmp]{formatN}} from package
  \pkg{gmp}.
}
\examples{
 ## Printing of MPFR numbers  uses formatMpfr() internally.
 ## Note how each components uses the "necessary" number of digits:
 ( x3 <- c(Const("pi", 168), mpfr(pi, 140), 3.14) )
 format(x3[3], 15)
 format(x3[3], 15, drop0 = TRUE)# "3.14" .. dropping the trailing zeros
 x3[4] <- 2^30
 x3[4] # automatically drops trailing zeros
 format(x3[1], dig = 41, small.mark = "'") # (41 - 1 = ) 40 digits after "."

 rbind(formatN(           x3,  digits = 15),
       formatN(as.numeric(x3), digits = 15))

 (Zero <- mpfr(c(0,1/-Inf), 20)) # 0 and "-0"
 xx <- c(Zero, 1:2, Const("pi", 120), -100*pi, -.00987)
 format(xx, digits = 2)
 format(xx, digits = 1, showNeg0 = FALSE)# "-0" no longer shown
}
\keyword{character}
\keyword{print}
