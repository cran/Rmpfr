\name{mpfrArray}
\alias{mpfrArray}
\title{Construct "mpfrArray" almost as by 'array()'}
\description{
  Utility to construct an \R object of class
  \code{\linkS4class{mpfrArray}}, very analogously to the numeric
  \code{\link{array}} function.
}
\usage{
mpfrArray(x, precBits, dim = length(x), dimnames = NULL)
}
\arguments{
  \item{x}{numeric(like) vector, typically of length \code{prod(dim)} or
    shorter in which case it is recycled.}
  \item{precBits}{ ~~Describe \code{precBits} here~~ }
  \item{dim}{the dimension of the array to be created, that is a
    vector of length one or more giving the maximal indices in
    each dimension.}
  \item{dimnames}{either \code{NULL} or the names for the
    dimensions.  This is a list with one component for each dimension,
    either \code{NULL} or a character vector of the length given by
    \code{dim} for that dimension.}
  %% ?array has more, about named dimnames etc...
}
\value{
  an object of class \code{"\linkS4class{mpfrArray}"}, specifically
 \code{"\linkS4class{mpfrMatrix}"} when \code{length(dim) == 2}.
}
\seealso{\code{\link{mpfr}}, \code{\link{array}}.
}
\examples{
## preallocating is possible here too
ma <- mpfrArray(NA,   prec = 80, dim = 2:4)
validObject(a2 <- mpfrArray(1:24, prec = 64, dim = 2:4))

## recycles, gives an "mpfrMatrix" and dimnames :
mat <- mpfrArray(1:5, 64, dim = c(5,3), dimnames=list(NULL, letters[1:3]))
mat
}
\keyword{array}
