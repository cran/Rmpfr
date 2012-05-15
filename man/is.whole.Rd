\name{is.whole}
\alias{is.whole.mpfr}
\title{Whole ("Integer") Numbers}
\description{
  Check which elements of \code{x[]} are integer valued aka
  \dQuote{whole} numbers,including MPFR
  numbers (class \code{\linkS4class{mpfr}}).
}
\usage{
\S3method{is.whole}{mpfr}(x)
}
\arguments{
  \item{x}{any \R vector, here of \code{\link{class}} \code{\linkS4class{mpfr}}.}
}
\value{
  logical vector of the same length as \code{x}, indicating where
  \code{x[.]} is integer valued.
}
\author{Martin Maechler}
\seealso{
 \code{\link{is.integer}(x)} (\pkg{base} package) checks for the
 \emph{internal} mode or class, not if \code{x[i]} are integer valued.

 The \code{\link[gmp]{is.whole}()} methods in package \pkg{gmp}.
}
\examples{
 is.integer(3) # FALSE, it's internally a double
 is.whole(3)   # TRUE
 x <- c(as(2,"mpfr") ^ 100, 3, 3.2, 1000000, 2^40)
 is.whole(x) # one FALSE, only
\dontshow{
xt <- c(x, as.bigz(2)^120, as.bigq(3,1:3))
stopifnot(identical(is.whole(x), (1:5)!=3), is(xt, "mpfr"),
          identical(is.whole(xt), is.na(match(1:9, c(3,8)))))
}%dont
}
\keyword{math}
