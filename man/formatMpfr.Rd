\name{formatMpfr}
\Rdversion{1.1}
\alias{formatMpfr}
\title{Formatting MPFR (multiprecision) Numbers}
\description{
  Flexible formatting of \dQuote{multiprecision numbers}, i.e., objects
  of class \code{\linkS4class{mpfr}}.  This is also the \code{mpfr}
  method of the generic \code{\link{format}} function.
}
\usage{
formatMpfr(x, digits = NULL, trim = FALSE, scientific = NA,
           big.mark = "", big.interval = 3L,
           small.mark = "", small.interval = 5L, decimal.mark = ".",
           zero.print = NULL, drop0trailing = FALSE, ...)
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
  \item{big.mark, big.interval, small.mark,
    small.interval, decimal.mark, zero.print, drop0trailing}{%
    used for prettying decimal sequences, these are passed to
    \code{\link{prettyNum}} and that help page explains the details.}
  \item{\dots}{further arguments passed to or from other methods.}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  a character vector of the same length as \code{x}.
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
}
\keyword{character}
\keyword{print}

