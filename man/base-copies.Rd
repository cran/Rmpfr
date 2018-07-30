\name{Rmpfr-workarounds}
\alias{outer}
\title{Base Functions etc, as an Rmpfr version}
\description{
  Functions from \pkg{base} etc which need a \emph{copy} in the
  \pkg{Rmpfr} namespace so they correctly dispatch.
}
\usage{
outer(X, Y, FUN = "*", ...)
}
\arguments{
  \item{X, Y, FUN, ...}{See \pkg{base} package help: \code{\link[base]{outer}}.}
}
\seealso{
 \code{\link[base]{outer}}.
}
\examples{
outer(1/mpfr(1:10, 70), 0:2)
}
\keyword{misc}
