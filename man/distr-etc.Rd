\name{mpfr-distr-etc}
\title{Distribution Functions etc (MPFR)}
\alias{pnorm}
\alias{dnorm}
\alias{dbinom}
\alias{dpois}
\alias{mpfr-distr}% <- alternative for \link
\usage{
dpois (x, lambda, log = FALSE)
dbinom(x, size, prob, log = FALSE)
dnorm (x, mean = 0, sd = 1, log = FALSE)

pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
}
\description{
  For some \R standard (probability) density, distribution or quantile
  functions, we provide MPFR versions.
}
\details{
  \code{pnorm()} is based on \code{\link{erf}()} and \code{erfc()} which
  have direct MPFR counter parts and are both reparametrizations
  of \code{pnorm}, \code{erf(x) = 2*pnorm(sqrt(2)*x)} and
  \code{erfc(x) = 2* pnorm(sqrt(2)*x, lower=FALSE)}.
}
\arguments{
  \item{x,q, lambda, size,prob, mean,sd}{\code{\link{numeric}} or
    \code{\linkS4class{mpfr}} vectors.  All of these are
    \dQuote{recycled} to the length of the longest one.}% for the function call.
  \item{log, log.p, lower.tail}{logical, see
    \code{\link[stats:Normal]{pnorm}},
    \code{\link[stats:Poisson]{dpois}}, etc.}
}
\value{
  A vector of the same length as the longest of \code{x,q, ...},
  of class \code{\linkS4class{mpfr}} with the high accuracy results of
  the corresponding standard \R function.
}
\seealso{\code{\link[stats:Normal]{pnorm}},
  \code{\link[stats]{dbinom}},
  \code{\link[stats]{dpois}} in standard package \pkg{stats}.

  \code{\link{pbetaI}(x, a,b)} is a \code{\link{mpfr}} version of
  \code{\link{pbeta}} only for \emph{integer} \code{a} and \code{b}.
}
\examples{
x <- 1400+ 0:10
print(dpois(x, 1000), digits =18) ## standard R's double precision
dpois(mpfr(x, 120), 1000)## more accuracy for the same
}
\keyword{distribution}

