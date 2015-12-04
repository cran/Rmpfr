\name{formatHex}
\title{Flexibly Format Numbers in Binary, Hex and Decimal Format}
\alias{formatHex}
\alias{formatBin}
\alias{formatDec}
%% \alias{formatHexInternal}
\alias{print.Bcharacter}% not mentioned
\description{
  Show numbers in binary, hex and decimal format.  The resulting
  character-like objects can be back-transformed to \code{"mpfr"}
  numbers via \code{\link{mpfr}()}.
}
\usage{
formatHex(x, precBits = min(getPrec(x)), style = "+")

formatBin(x, precBits = min(getPrec(x)), scientific = TRUE,
          left.pad = "_", right.pad = left.pad, style = "+")
%% FIXME:  Currently needed because decimal point alignment.
%% -----   TODO: Allow formatMpfr(.) to do decimal alignment -- even by default
formatDec(x, precBits = min(getPrec(x)), digits = decdigits,
          nsmall = NULL, scientific = FALSE, style = "+", ...)
%% ## not exported
%% formatHexInternal(x, precBits = min(Rmpfr::getPrec(x)), \dots)
}
\arguments{
  \item{x}{a \code{numeric} or \code{\link{mpfr}} \R object.}
  \item{precBits}{integer, the number of bits of precision, typically
    derived from \code{x}, see \code{\link{getPrec}}.  Numeric, i.e.,
    double precision numbers have 53 bits.  For more detail, see
    \code{\link[Rmpfr]{mpfr}}.}
  \item{style}{a single character, to be used in \code{\link{sprintf}}'s
    format (\code{fmt}), immediately after the "%".  The default always
    sets a sign in the output, i.e., \code{"+"} or \code{"-"}, where as
    \code{style = " "} may seem more standard.}
  \item{scientific}{\code{\link{logical}} indicating that
    \code{formatBin} should display the binary representation in
    scientific notation (\code{mpfr(3, 5)} is displayed as
    \code{+0b1.1000p+1}).  When \code{FALSE}, \code{formatBin} will
    display the binary representation in regular format shifted to align
    binary points (\code{mpfr(3, 5)} is displayed \code{+0b11.000}).}
  \item{\dots}{additional optional arguments.
    %% \code{formatHexInternal},
    \code{formatHex}, \code{formatBin}: \code{precBits} is the only
    \code{\dots} argument acted on.  Other \code{\dots} arguments are ignored.

    \code{formatDec}: \code{precBits} is acted on.  Any argument
    accepted by \code{\link[base]{format}} (except \code{nsmall}) is
    acted on.  Other \code{\dots} arguments are ignored.
  }
  \item{left.pad, right.pad}{characters (one-character strings) that
    will be used for left- and right-padding of the formatted string
    when \code{scientific=FALSE}.  \emph{Do not change these unless for
    display-only purpose !!}}
  \item{nsmall}{only used when \code{scientific} is false, then passed
    to \code{\link{format}()}.  If \code{NULL}, the default is computed
    from the range of the non-zero values of \code{x}.}
  \item{digits}{integer; the number of decimal digits displayed is the
    larger of this argument and the internally generated value that is a
    function of \code{precBits}.  This is related to but different than
    \code{digits} in \code{\link{format}}.}
}
\details{
  For the hexadecimal representation, when the precision is not larger
  than double precision, \code{\link[base]{sprintf}()} is used directly,
  otherwise \code{\link{formatMpfr}()} is used and post processed.
  For the binary representation, the hexadecimal value is calculated and
  then edited by
  substitution of the binary representation of the hex characters coded
  in the \code{HextoBin} vector.  For binary with \code{scientific=FALSE}, the
  result of the \code{scientific=TRUE} version is edited to align binary
  points.  For the decimal representation, the hexadecimal value is
  calculated with the specified precision and then sent to the
  \code{format} function for \code{scientific=FALSE} or to the sprintf
  function for \code{scientific=TRUE}.
%  \code{mpfrBchar} reads the output of \code{formatBin}.
%  \code{mpfrHchar} reads the output of \code{formatHex}.
}
\value{
  a character vector (or matrix) like \code{x}, say \code{r}, containing
  the formatted represention of \code{x}, with a \code{\link{class}}
  (unless \code{left.pad} or \code{right.pad} were not \code{"_"}).  In
  that case, \code{formatHex()} returns class \code{"Hcharacter"} and
  \code{formatBin()} gives class \code{"Bcharacter"}; for these,
  \code{\link{mpfr}(.)} has methods and will basically return \code{x},
  i.e., work as \emph{inverse} function.
}
\references{
  R FAQ 7.31:  Why doesn't R think these numbers are equal?
  \code{system.file("../../doc/FAQ")}
}
\author{Richard M. Heiberger \email{rmh@temple.edu}, with minor tweaking
  by Martin M.}
\seealso{
  \code{\link[Rmpfr]{mpfr}}, \code{\link[base]{sprintf}}
}
\examples{
FourBits <- mpfr(matrix(0:31, 8, 4, dimnames = list(0:7, c(0,8,16,24))),
                 precBits=4) ## 4 significant bits
FourBits

formatHex(FourBits)
formatBin(FourBits, style = " ")
formatBin(FourBits, scientific=FALSE)
formatDec(FourBits)

FBB <- formatBin(FourBits)
(nFBB <- mpfr(FBB))
stopifnot(all.equal(nFBB, FourBits, tol=0))

FBH <- formatHex(FourBits)
(nFBH <- mpfr(FBH))
stopifnot(all.equal(nFBH, FourBits, tol=0))

TenFrac <- matrix((1:10)/10, dimnames=list(1:10, expression(1/x)))
TenFrac9 <- mpfr(TenFrac, precBits=9) ## 9 significant bits
TenFrac9
formatHex(TenFrac9)
formatBin(TenFrac9)
formatBin(TenFrac9, scientific=FALSE)
formatDec(TenFrac9)
formatDec(TenFrac9, precBits=10)
}
\keyword{ arith }
