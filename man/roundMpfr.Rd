\name{roundMpfr}
\alias{roundMpfr}
\title{Rounding to binary bits, "mpfr-internally"}
\description{Rounding to binary bits, not decimal digits. Closer to the number
  representation, this also allows to increase or decrease a number's
  precBits.
}
\usage{
roundMpfr(x, precBits)
}
\arguments{
  \item{x}{an mpfr number (vector)}
  \item{precBits}{integer specifying the desired precision in bits.}
}
% \details{
% }
\value{an mpfr number as \code{x} but with the new 'precBits' precision}
\seealso{
  The \code{\linkS4class{mpfr}} class group method \code{Math2}
  implements a method for \code{\link{round}(x, digits)} which rounds to
  \emph{decimal} digits.
}
\examples{
 (p1 <- Const("pi", 100)) # 100 bit prec
 roundMpfr(p1, 120) # 20 bits more, but "random noise"
 Const("pi", 120)   # same "precision", but really precise
\dontshow{
 stopifnot(120 == getPrec(roundMpfr(p1, 120)),
            75 == getPrec(roundMpfr(p1,  75)))
}
}
\keyword{arith}
