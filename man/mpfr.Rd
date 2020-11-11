\name{mpfr}
\title{Create "mpfr" Numbers (Objects)}
\alias{mpfr}
\alias{is.mpfr}
\alias{mpfr.default}
\alias{mpfr.mpfr}% and more
\alias{Const}
\description{
  Create multiple (i.e. typically \emph{high}) precision numbers, to be
  used in arithmetic and mathematical computations with \R.
}
\usage{
mpfr(x, precBits, \dots)
\S3method{mpfr}{default}(x, precBits, base = 10,
     rnd.mode = c("N","D","U","Z","A"), scientific = NA, \dots)

Const(name = c("pi", "gamma", "catalan", "log2"), prec = 120L,
      rnd.mode = c("N","D","U","Z","A"))

is.mpfr(x)
}
\arguments{
  \item{x}{a \code{\link{numeric}}, \code{\linkS4class{mpfr}},
      \code{\link[gmp]{bigz}}, \code{\link[gmp]{bigq}}, or
      \code{\link{character}} vector or \code{\link{array}}.}
  \item{precBits, prec}{a number, the maximal precision to be used, in
    \bold{\emph{bits}}; i.e. \code{53} corresponds to double precision.
    Must be at least 2.  If \code{\link{missing}},
    \code{\link{getPrec}(x)} determines a default precision.}
  \item{base}{(only when \code{x} is \code{\link{character}}) the base
    with respect to which \code{x[i]} represent numbers; \code{base}
    \eqn{b} must fulfill \eqn{2 \le b \le 62}{2 <= b <= 62}.}
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see details.}
  \item{scientific}{(used only when \code{x} is the result of
    \code{\link{formatBin}()}, i.e., of class \code{"Bcharacter"}:)
    logical indicating that the binary representation of \code{x} is in
    scientific notation.  When \code{TRUE}, \code{mpfr()} will substitute
    \code{0} for \code{_}; when \code{NA}, \code{mpfr()} will guess, and
    use \code{TRUE} when finding a \code{"p"} in \code{x}; see also
    \code{\link{formatBin}}.}
  \item{name}{a string specifying the mpfrlib - internal constant
    computation.  \code{"gamma"} is Euler's gamma \eqn{(\gamma)}{}, and
    \code{"catalan"} Catalan's constant.}
  \item{\dots}{potentially further arguments passed to and from methods.}
}
\value{
  an object of (S4) class \code{\linkS4class{mpfr}}, or for
  \code{mpfr(x)} when \code{x} is an array,
  \code{\linkS4class{mpfrMatrix}}, or \code{\linkS4class{mpfrArray}}
  which the user should just as a normal numeric vector or array.

  \code{is.mpfr()} returns \code{TRUE} or \code{FALSE}.
}
\details{
  The \code{"\linkS4class{mpfr}"} method of \code{mpfr()} is a simple
      wrapper around \code{\link{roundMpfr}()}.

  MPFR supports the following rounding modes,
  %% the following is +- cut & paste from the MPFR manual
  \describe{
    \item{GMP_RND\bold{N}:}{round to \bold{n}earest (roundTiesToEven in
      IEEE 754-2008).}
    \item{GMP_RND\bold{Z}:}{round toward \bold{z}ero (roundTowardZero in
      IEEE 754-2008).}
    \item{GMP_RND\bold{U}:}{round toward plus infinity (\dQuote{Up},
      roundTowardPositive in IEEE 754-2008).}
    \item{GMP_RND\bold{D}:}{round toward minus infinity (\dQuote{Down},
      roundTowardNegative in IEEE 754-2008).}
    \item{GMP_RND\bold{A}:}{round \bold{a}way from zero (new since MPFR 3.0.0).}
  }
  The \sQuote{round to nearest} (\code{"N"}) mode, the default here,
  works as in the IEEE 754 standard: in case the number to be rounded
  lies exactly in the middle of two representable numbers, it is rounded
  to the one with the least significant bit set to zero.  For example,
  the number 5/2, which is represented by (10.1) in binary, is rounded
  to (10.0)=2 with a precision of two bits, and not to (11.0)=3.  This
  rule avoids the "drift" phenomenon mentioned by Knuth in volume 2 of
  The Art of Computer Programming (Section 4.2.2).

  When \code{x} is \code{\link{character}}, \code{mpfr()}
  will detect the precision of the input object. %% FIXME examples
}
\references{
  The MPFR team. (201x).
  \emph{GNU MPFR -- The Multiple Precision Floating-Point Reliable
    Library}; see \url{https://www.mpfr.org/mpfr-current/#doc} or directly
  \url{https://www.mpfr.org/mpfr-current/mpfr.pdf}.
}
\author{Martin Maechler}
\seealso{The class documentation \code{\linkS4class{mpfr}} contains more
  details.  Use \code{\link{asNumeric}} to transform back to double
  precision ("\code{\link{numeric}}").
}
\examples{
mpfr(pi, 120) ## the double-precision pi "translated" to 120-bit precision

pi. <- Const("pi", prec = 260) # pi "computed" to correct 260-bit precision
pi. # nicely prints 80 digits [260 * log10(2) ~= 78.3 ~ 80]

Const("gamma",   128L) # 0.5772...
Const("catalan", 128L) # 0.9159...

x <- mpfr(0:7, 100)/7 # a more precise version of  k/7, k=0,..,7
x
1 / x

## character input :
mpfr("2.718281828459045235360287471352662497757") - exp(mpfr(1, 150))
## ~= -4 * 10^-40
## Also works for  NA, NaN, ... :
cx <- c("1234567890123456", 345, "NA", "NaN", "Inf", "-Inf")
mpfr(cx)

## with some 'base' choices :
print(mpfr("111.1111", base=2)) * 2^4

mpfr("af21.01020300a0b0c", base=16)
##  68 bit prec.  44833.00393694653820642

mpfr("ugi0", base = 32) == 10^6   ## TRUE

## --- Large integers from package 'gmp':
Z <- as.bigz(7)^(1:200)
head(Z, 40)
## mfpr(Z) by default chooses the correct *maximal* default precision:
mZ. <- mpfr(Z)
## more efficiently chooses precision individually
m.Z <- mpfr(Z, precBits = frexpZ(Z)$exp)
## the precBits chosen are large enough to keep full precision:
stopifnot(identical(cZ <- as.character(Z),
                    as(mZ.,"character")),
          identical(cZ, as(m.Z,"character")))

## compare mpfr-arithmetic with exact rational one:
stopifnot(all.equal(mpfr(as.bigq(355,113), 99),
                    mpfr(355, 99) / 113,	tol = 2^-98))

## look at different "rounding modes":
sapply(c("N", "D","U","Z","A"), function(RND)
       mpfr(c(-1,1)/5, 20, rnd.mode = RND), simplify=FALSE)

symnum(sapply(c("N", "D","U","Z","A"),
              function(RND) mpfr(0.2, prec = 5:15, rnd.mode = RND) < 0.2 ))
}
\keyword{classes}
