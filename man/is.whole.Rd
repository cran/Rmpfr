\name{is.whole}
\Rdversion{1.1}
\alias{is.whole}
\title{Whole ("Integer") Numbers}
\description{
  Check which elements are integer valued numbers, including MPFR
  numbers (class \code{\linkS4class{mpfr}}).
}
\usage{
is.whole(x)
}
\arguments{
  \item{x}{any \R vector}
}
\value{
  logical vector of the same length as \code{x}, indicating where
  \code{x[.]} is integer valued.
}
\author{Martin Maechler}
\seealso{
 \code{\link{is.integer}(x)} (\pkg{base} package) checks for the
 \emph{internal} mode or class; not if \code{x[i]} are integer valued.
}
\examples{
 is.integer(3) # FALSE, it's internally a double
 is.whole(3)   # TRUE
 ## integer valued complex numbers  (two FALSE) :
 is.whole(c(7, 1 + 1i, 1.2, 3.4i, 7i))
 x <- c(as(2,"mpfr") ^ 100, 3, 3.2, 1000000, 2^40)
 is.whole(x) # one FALSE, only
}
\keyword{math}
