\name{mpfrArray}
\alias{mpfrArray}
\title{Construct "mpfrArray" almost as by 'array()'}
\description{
  Utility to construct an \R object of class
  \code{\linkS4class{mpfrArray}}, very analogously to the numeric
  \code{\link{array}} function.
}
\usage{
mpfrArray(x, precBits, dim = length(x), dimnames = NULL,
          rnd.mode = c("N","D","U","Z","A"))
}
\arguments{
  \item{x}{numeric(like) vector, typically of length \code{prod(dim)} or
    shorter in which case it is recycled.}
  \item{precBits}{a number, the maximal precision to be used, in
    \bold{\emph{bits}}; i.e., \code{53} corresponds to double precision.
    Must be at least 2.}
  \item{dim}{the dimension of the array to be created, that is a
    vector of length one or more giving the maximal indices in
    each dimension.}
  \item{dimnames}{either \code{NULL} or the names for the
    dimensions.  This is a list with one component for each dimension,
    either \code{NULL} or a character vector of the length given by
    \code{dim} for that dimension.}
  %% ?array has more, about named dimnames etc...
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see details of
    \code{\link{mpfr}}.}
}
\value{
  an object of class \code{"\linkS4class{mpfrArray}"}, specifically
 \code{"\linkS4class{mpfrMatrix}"} when \code{length(dim) == 2}.
}
\seealso{\code{\link{mpfr}}, \code{\link{array}};
  \code{\link{asNumeric}()} as \dQuote{inverse} of \code{mpfrArray()},
  to get back a numeric array.

  \code{\link{mpfr2array}(x)} is for \code{"mpfr"} classed \code{x},
  only, whereas \code{mpfrArray(x)} is for numeric (\dQuote{non-mpfr})
  \code{x}.
}
\examples{
## preallocating is possible here too
ma <- mpfrArray(NA, prec = 80, dim = 2:4)
validObject(A2 <- mpfrArray(1:24, prec = 64, dim = 2:4))

## recycles, gives an "mpfrMatrix" and dimnames :
mat <- mpfrArray(1:5, 64, dim = c(5,3), dimnames=list(NULL, letters[1:3]))
mat
asNumeric(mat)
stopifnot(identical(asNumeric(mat),
                    matrix(1:5 +0, 5,3, dimnames=dimnames(mat))))

## Testing the apply() method :
apply(mat, 2, range)
apply(A2, 1:2, range)
apply(A2, 2:3, max)
(fA2 <- apply(A2, 2, fivenum))
a2 <- as(A2, "array")
stopifnot(as(apply(A2, 2, range), "matrix") ==
             apply(a2, 2, range)
        , all.equal(fA2, apply(a2, 2, fivenum))
        , all.equal(apply(A2, 2, quantile),
                    apply(a2, 2, quantile))
        , all.equal(A2, apply(A2, 2:3, identity) -> aA2, check.attributes=FALSE)
        , dim(A2) == dim(aA2)
)
\dontshow{
 for(nf in c("colSums", "colMeans", "rowSums", "rowMeans")) {
   FUN <- getFunction(nf)
   for(di in c(1,2)) {
     r <- FUN(a2, dims = di)
     R <- FUN(A2, dims = di)
     stopifnot(dim(r) == dim(R),
               all.equal(as(R, if(is.array(r)) "array" else "numeric"),
                         unname(r), tol = 1e-15))
   }
 }
## with non-trivial dimnames:
ma2 <- mat
dimnames(ma2) <- list(row=paste0("r",1:5), col=colnames(mat))
stopifnot(identical(ma2,   apply(ma2, 2, identity)),
          identical(ma2, t(apply(ma2, 1, identity))),
          identical(names(apply(ma2,2,sum)), colnames(ma2))
)
}%end{ dontshow }
}
\keyword{array}
