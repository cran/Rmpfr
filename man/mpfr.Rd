\name{mpfr}
\alias{mpfr}
\alias{Const}
\title{Create "mpfr" Numbers (Objects)}
\description{
  Create multiple (i.e. typically \emph{high}) precision numbers, to be
  used in arithmetic and mathematical computations with \R.
}
\usage{
mpfr(x, precBits, base = 10)
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
  \item{name}{a string specifying the mpfrlib - internal constant
    computation.  \code{"gamma"} is Euler's gamma \eqn{(\gamma)}{}, and
    \code{"catalan"} Catalan's constant.}
}
\value{
  an object of (S4) class \code{\linkS4class{mpfr}} which the user
  should just as a normal numeric vector.
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
}
\keyword{classes}
