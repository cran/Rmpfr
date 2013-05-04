\name{asNumeric-methods}
\docType{methods}
%% FIXME??  Should this migrate to base ?? --> have implicit generic and S3 generic there!
%% ----- there, with semantic  storage.mode(x) <- "numeric"
\alias{asNumeric-methods}
\alias{asNumeric,mpfr-method}
\alias{asNumeric,mpfrArray-method}
\title{Methods for \code{asNumeric(<mpfr>)}}
\description{% -->  ~/R/Pkgs/gmp/man/asNumeric.Rd
  Methods for function \code{\link[gmp]{asNumeric}} (in package \pkg{gmp}).
}
\usage{
%\S4method{asNumeric}{mpfr}(x)
\S4method{asNumeric}{mpfrArray}(x)
}
\arguments{
  \item{x}{a \dQuote{number-like} object, here, a
    \code{\linkS4class{mpfr}} or typically
    \code{\linkS4class{mpfrArray}}one.}
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "mpfrArray")}}{this method also dispatches
      for \code{\linkS4class{mpfrMatrix}} and returns a numeric array.}
    \item{\code{signature(x = "mpfr")}}{for non-array/matrix,
      \code{asNumeric(x)} is basically the same as \code{as.numeric(x)}.}
  }
}
\value{
  an \R object of type (\code{\link{typeof}}) \code{"numeric"}, a \code{\link{matrix}}
  or \code{\link{array}} if \code{x} had non-NULL dimension \code{\link{dim}()}.
}
\author{Martin Maechler}
\seealso{
  our lower level (non-generic) \code{\link{toNum}()}. Further,
  \code{\link[gmp]{asNumeric}} (package \pkg{gmp}),
  standard \R's \code{\link{as.numeric}()}.
}
\examples{
x <- (0:7)/8 # (exact)
X <- mpfr(x, 99)
stopifnot(identical(asNumeric(x), x),
	  identical(asNumeric(X), x))

m <- matrix(1:6, 3,2)
(M <- mpfr(m, 99) / 5) ##-> "mpfrMatrix"
asNumeric(M) # numeric matrix
stopifnot(all.equal(asNumeric(M), m/5),
          identical(asNumeric(m), m))# remains matrix
}
\keyword{methods}

