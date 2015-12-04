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
	   base = 10, showNeg0 = TRUE,
	   big.mark = "", big.interval = 3L,
	   small.mark = "", small.interval = 5L,
           decimal.mark = ".",
           exponent.char = if(base <= 14) "e" else if(base <= 36) "E" else "|e",
           zero.print = NULL, drop0trailing = FALSE, ...)
\S3method{formatN}{mpfr}(x, drop0trailing = TRUE, \dots)
}
\arguments{
  \item{x}{an MPFR number (vector or array).}
  \item{digits}{how many significant digits (in the \code{base} chosen!)
    are to be used in the result.  The default, \code{NULL}, uses enough
    digits to represent the full precision, often one or two digits more
    than \dQuote{you} would expect.}
  \item{trim}{logical; if \code{FALSE}, numbers are right-justified to a
    common width: if \code{TRUE} the leading blanks for justification
    are suppressed.}
  \item{scientific}{either a logical specifying whether
    MPFR numbers should be encoded in scientific
    format, or an integer penalty (see \code{\link{options}("scipen")}).
    Missing values correspond to the current default penalty.}
  \item{base}{an integer in \eqn{2,3,..,62}; the base (\dQuote{basis})
    in which the numbers should be represented.  Apart from the default
    base 10, binary (\code{base = 2}) or hexadecimal (\code{base = 16}) are particularly
    interesting.}
  \item{showNeg0}{logical indicating if \dQuote{\bold{neg}ative} zeros
    should be shown with a \code{"-"}.  The default, \code{TRUE} is
    intentially different from \code{\link{format}(<numeric>)}.}
  \item{exponent.char}{the \dQuote{exponent} character to be used in
    scientific notation.  The default takes into account that for
    \code{base} \eqn{B \ge 15}{B >= 15}, \code{"e"} is part of the
    (mantissa) digits and the same is true for \code{"E"} when
    \eqn{B \ge 37}{B >= 37}.}
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

  Note that in scientific notation, the integer exponent is always in
  \emph{decimal}, i.e., base 10 (even when \code{base} is not 10), but
  of course meaning \code{base} powers, e.g., in base 32,
  \code{"u.giE3"}is the same as \code{"ugi0"} which is \eqn{32^3} times
  \code{"u.gi"}. This is in contrast, e.g., with
  \code{\link{sprintf}("\%a", x)} where the powers after \code{"p"} are
  powers of \eqn{2}.
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

## Output in other bases :
formatMpfr(mpfr(10^6, 40), base=32, drop0trailing=TRUE)
## "ugi0"
mpfr("ugi0", base=32) #-> 1'000'000
\dontshow{
stopifnot(
   identical("ugi0", formatMpfr(mpfr(10^6, 40), base=32, drop0trailing=TRUE)),
   mpfr("ugi0", base=32) == 10^6)
}

i32 <- mpfr(1:32, precBits = 64)
format(i32,   base=  2, drop0trailing=TRUE)
format(i32,   base= 16, drop0trailing=TRUE)
format(1/i32, base=  2, drop0trailing=TRUE)# using scientific notation for [17..32]
format(1/i32, base= 32)
format(1/i32, base= 62, drop0trailing=TRUE)
format(mpfr(2, 64)^-(1:16), base=16, drop0trailing=TRUE)
}
\keyword{character}
\keyword{print}
