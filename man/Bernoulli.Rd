\name{Bernoulli}
\alias{Bernoulli}
\title{Bernoulli Numbers in Arbitrary Precision}
\description{
  Computes the Bernoulli numbers in the desired (binary) precision.
  The computation happens via the \code{\link{zeta}} function and the
  formula
  \deqn{B_k = -k \zeta(1 - k),}
  and hence the only non-zero odd Bernoulli number is \eqn{B_1 = +1/2}.
  (Another tradition defines it, equally sensibly, as \eqn{-1/2}.)
}
\usage{
Bernoulli(k, precBits = 128)
}
\arguments{
  \item{k}{non-negative integer vector}
  \item{precBits}{the precision in \emph{bits} desired.}
}
\value{
  an \code{\linkS4class{mpfr}} class vector of the same length as
  \code{k}, with i-th component the \code{k[i]}-th Bernoulli number.
}
\references{\url{https://en.wikipedia.org/wiki/Bernoulli_number}
}
\author{Martin Maechler}
\seealso{\code{\link{zeta}} is used to compute them.

  The next version of package \CRANpkg{gmp} is to contain %\code{\link[gmp]{BernoulliQ}()}
  \code{BernoulliQ()}, providing exact Bernoulli numbers as
  big rationals (class \code{"bigq"}).
}
\examples{
\dontshow{sessionInfo()
 .libPaths()
 packageDescription("gmp")}
Bernoulli(0:10)
plot(as.numeric(Bernoulli(0:15)), type = "h")

curve(-x*zeta(1-x), -.2, 15.03, n=300,
      main = expression(-x \%.\% zeta(1-x)))
legend("top", paste(c("even","odd  "), "Bernoulli numbers"),
       pch=c(1,3), col=2, pt.cex=2, inset=1/64)
abline(h=0,v=0, lty=3, col="gray")
k <- 0:15; k[1] <- 1e-4
points(k, -k*zeta(1-k), col=2, cex=2, pch=1+2*(k\%\%2))

## They pretty much explode for larger k :
k2 <- 2*(1:120)
plot(k2, abs(as.numeric(Bernoulli(k2))), log = "y")
title("Bernoulli numbers exponential growth")

Bernoulli(10000)# - 9.0494239636 * 10^27677
}
\keyword{arith}
