\name{sapplyMpfr}
\alias{sapplyMpfr}
\title{Apply a Function over a "mpfr" Vector}
\description{
  Users may be disappointed to note that \code{\link{sapply}()} or
  \code{\link{vapply}()} typically do not work with \code{"mpfr"}
  numbers.

  This is a simple (but strong) approach to work around the problem,
  based on \code{\link{lapply}()}.
}
\usage{
sapplyMpfr(X, FUN, \dots, drop_1_ = TRUE)
}
\arguments{
  \item{X}{a vector, possibly of class \code{"\linkS4class{mpfr}"}.}
  \item{FUN}{a \code{\link{function}} returning an \code{"\linkS4class{mpfr}"}
    vector or even an \code{"\linkS4class{mpfrArray}"}.  May also be a
    function returning a \code{\link{numeric}} vector or array for numeric
    \code{X}, \emph{and} which returns \code{"mpfr(Array)"} for an \code{X}
    argument inheriting from \code{"\linkS4class{mpfr}"}.}
  \item{\dots}{further arguments passed to \code{\link{lapply}},
    typically further arguments to \code{FUN}.}
  \item{drop_1_}{logical (with unusual name on purpose!) indicating if
    1-column matrices (\code{"mpfrMatrix"}) should be \dQuote{dropped} to
    vectors (\code{"mpfr"}), the same as in base \R's own \code{\link{sapply}}.
    This has been implicitly \code{FALSE} in \pkg{Rmpfr} versions 0.8-5 to
    0.8-9 (Oct 2021 to June 2022), accidentally.
    Since \pkg{Rmpfr} 0.9-0, this has been made an argument with default
    \code{TRUE} to be compatible by default with \R's \code{\link{sapply}}.
  }
}
\details{
  In the case \code{FUN(<length-1>)} returns an \code{\link{array}} or \code{"mpfrArray"}, i.e.,
  with two or more dimensions, \code{sapplyMpfr()} returns an
  \code{"mpfrArray"}; this is analogous to \code{\link{sapply}(X, FUN, simplify = "array")}
  (rather than the default \code{sapply()} behaviour which returns a
   \code{matrix} also when a higher array would be more \dQuote{logical}.)
}
\value{
  an  \code{"\linkS4class{mpfr}"} vector, typically of the same length
  as \code{X}.
}
\author{Martin Maechler}
\note{
  This may still not always work as well as \code{\link{sapply}()} does for
  atomic vectors.  The examples seem to indicate that it typically does
  work as desired, since \pkg{Rmpfr} version 0.9-0.

  If you want to transform back to regular numbers anyway, it maybe simpler
  and more efficient to use \preformatted{
    res <- lapply(....)
    sapply(res, asNumeric, simplify = "array")
  }
  instead of \code{sapplyMpfr()}.
}
\seealso{
  \code{\link{sapply}}, \code{\link{lapply}}, etc.
}
\examples{
sapplyMpfr0 <- ## Originally, the function was simply defined as
  function (X, FUN, ...) new("mpfr", unlist(lapply(X, FUN, ...), recursive = FALSE))

(m1 <- sapply    (     3,      function(k) (1:3)^k)) # 3 x 1  matrix (numeric)
(p1 <- sapplyMpfr(mpfr(3, 64), function(k) (1:3)^k))
stopifnot(m1 == p1, is(p1, "mpfrMatrix"), dim(p1) == c(3,1), dim(p1) == dim(m1))
k.s <- c(2, 5, 10, 20)
(mk <- sapply    (     k.s,      function(k) (1:3)^k)) # 3 x 4    "       "
(pm <- sapplyMpfr(mpfr(k.s, 64), function(k) (1:3)^k))
stopifnot(mk == pm, is(pm, "mpfrMatrix"), dim(pm) == 3:4, 3:4 == dim(mk))
## was *wrongly* 4x3  in Rmpfr 0.8-x
f5k  <- function(k) outer(1:5, k+0:2, `^`)# matrix-valued
(mk5 <- sapply    (     k.s,      f5k))  # sapply()'s default; not "ideal"
(ak5 <- sapply    (     k.s,      f5k, simplify = "array")) # what we want
(pm5 <- sapplyMpfr(mpfr(k.s, 64), f5k))
stopifnot(c(mk5) == c(ak5), ak5 == pm5, is(pm5, "mpfrArray"), is.array(ak5),
          dim(pm5) == dim(ak5), dim(pm5) == c(5,3, 4))
if(require("Bessel")) { # here X, is simple
  bI1 <- function(k) besselI.nuAsym(mpfr(1.31e9, 128), 10, expon.scaled=TRUE, k.max=k)
  bImp1 <- sapplyMpfr (0:4, bI1, drop_1_ = FALSE) # 1x5 mpfrMatrix -- as in DPQ 0.8-8
  bImp  <- sapplyMpfr (0:4, bI1, drop_1_ = TRUE ) # 5 "mpfr" vector {by default}
  bImp0 <- sapplyMpfr0(0:4, bI1) # 5-vector
  stopifnot(identical(bImp, bImp0), bImp == bImp1,
            is(bImp, "mpfr"), is(bImp1, "mpfrMatrix"), dim(bImp1) == c(1, 5))
}# {Bessel}
}
\keyword{manip}
