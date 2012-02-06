\name{sumBinomMpfr}
\alias{sumBinomMpfr}
\title{(Alternating) Binomial Sums via Rmpfr}
\concept{Rice integral}
\concept{Forward Difference}
\description{
  Compute (alternating) binomial sums via high-precision arithmetic.
  Such sums appear in different contexts and are typically challenging,
  i.e., currently impossible, to evaluate reliably as soon as \eqn{n} is
  larger than around \eqn{50--70}.

  The alternating binomial sum \eqn{sB(f,n) := sumBinom(n, f, n0=0)} is
  (up to sign) equal to the \eqn{n}-th forward difference operator \eqn{\Delta^n f},
  \deqn{sB(f,n) = (-1)^n \Delta^n f,} where
  \deqn{\Delta^n f = \sum_{k=0}^{n} (-1)^{n-k}{n \choose k}\cdot f(k),}{%
    Delta^n f = sum(k = n0:n ; (-1)^(n-k) choose(n,k) * f(k)),}
  is the \eqn{n}-fold iterated forward difference
  \eqn{\Delta f(x) = f(x+1) - f(x)} (for \eqn{x = 0}).
}
\usage{
sumBinomMpfr(n, f, n0 = 0, alternating = TRUE, precBits = 256)
}
\arguments{
  \item{n}{upper summation index (integer).}
  \item{f}{\code{\link{function}} to be evaluated at \eqn{k} for
    \code{k in n0:n}.}
  \item{n0}{lower summation index, typically \code{0} (= default) or \code{1}.}
  \item{alternating}{logical indicating if the sum is alternating, see below.}
  \item{precBits}{the number of bits for MPFR precision, see \code{\link{mpfr}}.}
}
\details{
  The current implementation is not particularly fast and might be
  improved in the future, notably for the case where
  \eqn{sB(f,n)=}\code{sumBinomMpfr(n, f, *)} is to be computed for a whole sequence
  \eqn{n = 1,\dots,N}{n = 1,...,N}.
}
\value{
  an \code{\linkS4class{mpfr}} number of precision \code{precBits}, say
  \eqn{s}. If \code{alternating} is true (as per default),
  \deqn{s = \sum_{k = n0}^n (-1)^k {n \choose k}\cdot f(k),}{%
    s = sum(k = n0:n ; (-1)^k choose(n,k) * f(k)),}
  if \code{alternating} is false, the \eqn{(-1)^k} factor is dropped (or
  replaced by \eqn{1}) above.
}
\references{
  Wikipedia (2012)
  The N\"orlund-Rice integral,
  \url{http://en.wikipedia.org/wiki/Rice_integral}

  Flajolet, P. and Sedgewick, R. (1995)
  Mellin Transforms and Asymptotics: Finite Differences and Rice's Integrals,
  \emph{Theoretical Computer Science} \bold{144}, 101--124.
}
\author{Martin Maechler, after conversations with Christophe Dutang.}
\seealso{
  \code{\link{chooseMpfr}}
}
\examples{
## "naive" R implementation:
sumBinom <- function(n, f, n0=0, ...) {
  k <- n0:n
  sum( choose(n, k) * (-1)^k * f(k, ...))
}

## compute  sumBinomMpfr(.) for a whole set of 'n' values:
sumBin.all <- function(n, f, n0=0, precBits = 256, ...)
{
  N <- length(n)
  precBits <- rep(precBits, length = N)
  ll <- lapply(seq_len(N), function(i)
           sumBinomMpfr(n[i], f, n0=n0, precBits=precBits[i], ...))
  sapply(ll, as, "double")
}
sumBin.all.R <- function(n, f, n0=0, ...)
   sapply(n, sumBinom, f=f, n0=n0, ...)

n.set <- 5:80
system.time(res.R   <- sumBin.all.R(n.set, f = sqrt)) ## instantaneous..
system.time(resMpfr <- sumBin.all  (n.set, f = sqrt)) ## takes ~2 seconds

matplot(n.set, cbind(res.R, resMpfr), type = "l", lty=1,
        ylim = extendrange(resMpfr, f = 0.25), xlab = "n",
        main = "sumBinomMpfr(n, f = sqrt)  vs.  R double precision")
legend("topleft", leg=c("double prec.", "mpfr"), lty=1, col=1:2, bty = "n")
}
\keyword{arith}
