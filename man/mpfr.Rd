\name{mpfr}
\alias{mpfr}
\alias{Const}
\title{Create "mpfr" Numbers (Objects)}
\description{
  Create multiple (i.e. typically \emph{high}) precision numbers, to be
  used in arithmetic and mathematical computations with \R.
}
\usage{
mpfr(x, precBits, base = 10, rnd.mode = c("N","D","U","Z"))
Const(name = c("pi", "gamma", "catalan"), prec = 120L)
}
\arguments{
  \item{x}{a \code{\link{numeric}} or \code{\link{character}} vector.}
  \item{precBits, prec}{a number, the maximal precision to be used, in
    \bold{\emph{bits}}; i.e. \code{53} corresponds to double precision.
    Must be at least 2.}
  \item{base}{(only when \code{x} is \code{\link{character}}) the base
    with respect to which \code{x[i]} represent numbers; \code{base}
    \eqn{b} must fulfill \eqn{2 \le b \le 36}{2 <= b <= 36}.}
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see details.}
  \item{name}{a string specifying the mpfrlib - internal constant
    computation.  \code{"gamma"} is Euler's gamma \eqn{(\gamma)}{}, and
    \code{"catalan"} Catalan's constant.}
}
\value{
  an object of (S4) class \code{\linkS4class{mpfr}} which the user
  should just as a normal numeric vector.
}
\details{
  MPFR supports four rounding modes,
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
  }
  The \sQuote{round to nearest} (\code{"N"}) mode, the default here,
  works as in the IEEE 754 standard: in case the number to be rounded
  lies exactly in the middle of two representable numbers, it is rounded
  to the one with the least significant bit set to zero.  For example,
  the number 5/2, which is represented by (10.1) in binary, is rounded
  to (10.0)=2 with a precision of two bits, and not to (11.0)=3.  This
  rule avoids the "drift" phenomenon mentioned by Knuth in volume 2 of
  The Art of Computer Programming (Section 4.2.2).

}
\references{
  The MPFR team. (2009).
  \emph{GNU MPFR -- The Multiple Precision Floating-Point Reliable
    Library}; Edition 2.4.2, November 2009;
  see \url{http://www.mpfr.org/mpfr-current/#doc} or directly
  \url{http://www.mpfr.org/mpfr-current/mpfr.pdf}.
}
\author{Martin Maechler}
\seealso{The class documentation \code{\linkS4class{mpfr}} contains more
  details.
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

## with some 'base' choices :
print(mpfr("111.1111", base=2)) * 2^4

mpfr("af21.01020300a0b0c", base=16)
##  68 bit prec.  44833.00393694653820642

## look at different "rounding modes":
sapply(c("N", "D","U","Z"), function(RND)
       mpfr(c(-1,1)/5, 20, rnd.mode = RND), simplify=FALSE)

symnum(sapply(c("N", "D","U","Z"),
              function(RND) mpfr(0.2, prec = 5:15, rnd.mode = RND) < 0.2 ))
}
\keyword{classes}
