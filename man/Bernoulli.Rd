\name{Bernoulli}
\alias{Bernoulli}
\title{Bernoulli Numbers in Arbitrary Precision}
\description{
  Computes the Bernoulli numbers in the desired (binary) precision.
  The computation happens via the \code{\link{zeta}} function and the
  formula
  \deqn{B_k = -k \zeta(1 - k).}
}
\usage{
Bernoulli(k, precBits = 128)
}
\arguments{
  \item{k}{positive integer vector}
  \item{precBits}{the precision in \emph{bits} desired.}
}
\value{
  an \code{\linkS4class{mpfr}} class vector of the same length as
  \code{k}, with i-th component the \code{k[i]}-th Bernoulli number.
}
\references{\url{http://en.wikipedia.org/wiki/Bernoulli_number}
}
\author{Martin Maechler}
\seealso{\code{\link{zeta}} is used to compute them.
}
\examples{
Bernoulli(1:10)
plot(as.numeric(Bernoulli(1:15)), type = "h")

## They pretty much explode for larger k :
k2 <- 2*(1:120)
plot(k2, abs(as.numeric(Bernoulli(k2))), log = "y")
title("Bernoulli numbers exponential growth")

Bernoulli(10000)# - 9.0494239636 * 10^27677
}
\keyword{arith}
