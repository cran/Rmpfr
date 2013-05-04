\name{mpfrMatrix}
\title{Classes "mpfrMatrix" and "mpfrArray"}
\docType{class}
\alias{mpfrMatrix-class}
\alias{mpfrArray-class}
%
\alias{Arith,mpfrArray,mpfr-method}
\alias{Arith,mpfrArray,mpfrArray-method}
\alias{Arith,mpfr,mpfrArray-method}
\alias{Arith,mpfrArray,numeric-method}
\alias{Arith,numeric,mpfrArray-method}
\alias{Compare,mpfrArray,mpfr-method}
\alias{Compare,mpfrArray,numeric-method}
\alias{Compare,mpfr,mpfrArray-method}
\alias{Compare,numeric,mpfrArray-method}
%
\alias{apply,mpfrArray-method}
\alias{colSums,mpfrArray-method}
\alias{colMeans,mpfrArray-method}
\alias{rowSums,mpfrArray-method}
\alias{rowMeans,mpfrArray-method}
\alias{as.vector,mpfrArray,missing-method}
\alias{coerce,mpfrArray,array-method}
\alias{coerce,mpfrArray,matrix-method}
\alias{coerce,mpfrMatrix,matrix-method}
%% the following two from the same setAs(); the 2nd one is auto-produced
%% *and* needed by codoc().  Behavior seems a bit bogous :
\alias{coerce,mpfrArray,vector-method}
\alias{coerce<-,mpfrArray,vector-method}
\alias{diag,mpfrMatrix-method}
\alias{diag<-,mpfrMatrix-method}
%
\alias{[<-,mpfrArray,ANY,ANY,ANY-method}
\alias{[<-,mpfrArray,ANY,ANY,mpfr-method}
\alias{[<-,mpfrArray,ANY,missing,ANY-method}
\alias{[<-,mpfrArray,ANY,missing,mpfr-method}
\alias{[<-,mpfrArray,missing,ANY,ANY-method}
\alias{[<-,mpfrArray,missing,ANY,mpfr-method}
\alias{[<-,mpfrArray,missing,missing,ANY-method}
\alias{[<-,mpfrArray,missing,missing,mpfr-method}
\alias{[<-,mpfrArray,matrix,missing,ANY-method}
\alias{[<-,mpfrArray,matrix,missing,mpfr-method}
%
\alias{[,mpfrArray,ANY,ANY,ANY-method}
\alias{[,mpfrArray,ANY,missing,missing-method}
\alias{[,mpfrArray,matrix,missing,missing-method}

\alias{\%*\%,Mnumber,mpfr-method}
\alias{\%*\%,mpfrMatrix,mpfr-method}
\alias{\%*\%,mpfrMatrix,mpfrMatrix-method}
\alias{\%*\%,mpfr,Mnumber-method}
\alias{\%*\%,mpfr,mpfr-method}
\alias{\%*\%,mpfr,mpfrMatrix-method}
\alias{crossprod,mpfr,missing-method}
\alias{crossprod,mpfrMatrix,mpfr-method}
\alias{crossprod,mpfrMatrix,mpfrMatrix-method}
\alias{crossprod,mpfr,Mnumber-method}
\alias{crossprod,mpfr,mpfr-method}
\alias{crossprod,mpfr,mpfrMatrix-method}
\alias{crossprod,Mnumber,mpfr-method}
\alias{tcrossprod,mpfr,missing-method}
\alias{tcrossprod,mpfrMatrix,mpfr-method}
\alias{tcrossprod,mpfrMatrix,mpfrMatrix-method}
\alias{tcrossprod,mpfr,Mnumber-method}
\alias{tcrossprod,mpfr,mpfr-method}
\alias{tcrossprod,mpfr,mpfrMatrix-method}
\alias{tcrossprod,Mnumber,mpfr-method}

\alias{dim,mpfrArray-method}
\alias{dimnames<-,mpfrArray-method}
\alias{dimnames,mpfrArray-method}

\alias{show,mpfrArray-method}
\alias{sign,mpfrArray-method}
\alias{t,mpfrMatrix-method}
\alias{aperm,mpfrArray-method}
%
\description{The classes \code{"mpfrMatrix"} and \code{"mpfrArray"} are,
  analogously to the \pkg{base} \code{\link{matrix}} and \code{\link{array}}
  functions and classes simply \dQuote{numbers} of class
  \code{\linkS4class{mpfr}} with an additional \code{Dim} and
  \code{Dimnames} slot.
}
\section{Objects from the Class}{
  Objects should typically be created by \code{\link{mpfrArray}()}, but
  can also be created by
  \code{new("mpfrMatrix", ...)} or \code{new("mpfrArray", ...)}, or also
  by \code{t(x)}, \code{dim(x) <- dd}, or \code{\link{mpfr2array}(x,
    dim=dd)} where \code{x} is a an \code{\linkS4class{mpfr}}
  \dQuote{number vector}.

  A (slightly more flexible) alternative to \code{dim(x) <- dd} is
  \code{\link{mpfr2array}(x, dd, dimnames)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{as for the \code{\linkS4class{mpfr}} class, a
      \code{"list"} of  \code{\linkS4class{mpfr1}} numbers.}
    \item{\code{Dim}:}{of class \code{"integer"}, specifying the array
      dimension.}
    \item{\code{Dimnames}:}{of class \code{"list"} and the same length
      as \code{Dim}, each list component either \code{\link{NULL}} or a
      \code{\link{character}} vector of length \code{Dim[j]}.}
  }
}
\section{Extends}{
Class \code{"mpfrMatrix"} extends \code{"mpfrArray"}, directly.

Class \code{"mpfrArray"} extends
class \code{"\linkS4class{mpfr}"}, by class "mpfrArray", distance 2;
class \code{"\linkS4class{list}"}, by class "mpfrArray", distance 3;
class \code{"\linkS4class{vector}"}, by class "mpfrArray", distance 4.
}
\section{Methods}{
  \describe{
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "mpfrArray")}: ... }
    \item{Arith}{\code{signature(e1 = "numeric", e2 = "mpfrArray")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfrArray", e2 = "mpfrArray")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfrArray", e2 = "mpfr")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfrArray", e2 = "numeric")}: ... }
    \item{as.vector}{\code{signature(x = "mpfrArray", mode =
	"missing")}: drops the dimension \sQuote{attribute}, i.e.,
      transforms \code{x} into a simple \code{\linkS4class{mpfr}}
      vector.  This is an inverse of \code{t(.)} or \code{dim(.) <- *}
      on such a vector.}

    \item{atan2}{\code{signature(y = "ANY", x = "mpfrArray")}: ... }
    \item{atan2}{\code{signature(y = "mpfrArray", x = "mpfrArray")}: ... }
    \item{atan2}{\code{signature(y = "mpfrArray", x = "ANY")}: ... }
    \item{[<-}{\code{signature(x = "mpfrArray", i = "ANY", j = "ANY", value = "ANY")}: ... }
    \item{[}{\code{signature(x = "mpfrArray", i = "ANY", j = "ANY", drop = "ANY")}: ... }
    \item{[}{\code{signature(x = "mpfrArray", i = "ANY", j = "missing", drop = "missing")}:
      \code{"mpfrArray"}s can be subset (\dQuote{indexed}) as regular \R
      \code{\link{array}}s.}
    \item{\%*\%}{\code{signature(x = "mpfr", y = "mpfrMatrix")}: Compute
      the matrix/vector product \eqn{x y} when the dimensions
      (\code{\link{dim}}) of \code{x} and \code{y} match.  If \code{x}
      is not a matrix, it is treated as a 1-row or 1-column matrix (aka
      \dQuote{row vector} or \dQuote{column vector}) depending on which
      one makes sense, see the documentation of the \pkg{base}
      function \code{\link[base:matmult]{\%*\%}}.}
    \item{\%*\%}{\code{signature(x = "mpfr", y = "Mnumber")}: method
      definition for cases with one \code{\linkS4class{mpfr}} and any
      \dQuote{number-like} argument are to use MPFR arithmetic as well.}
    \item{\%*\%}{\code{signature(x = "mpfrMatrix", y = "mpfrMatrix")},}
    \item{\%*\%}{\code{signature(x = "mpfrMatrix", y = "mpfr")}, etc.
      Further method definitions with identical semantic.}

    \item{crossprod}{\code{signature(x = "mpfr", y = "missing")}:
      Computes \eqn{x'x}, i.e., \code{t(x) \%*\% x}, typically more efficiently.}
    \item{crossprod}{\code{signature(x = "mpfr", y = "mpfrMatrix")}:
      Computes \eqn{x'y}, i.e., \code{t(x) \%*\% y}, typically more efficiently.}
    \item{crossprod}{\code{signature(x = "mpfrMatrix", y = "mpfrMatrix")}: ... }
    \item{crossprod}{\code{signature(x = "mpfrMatrix", y = "mpfr")}: ... }

    \item{tcrossprod}{\code{signature(x = "mpfr", y = "missing")}:
      Computes \eqn{xx'}, i.e., \code{x \%*\% t(x)}, typically more efficiently.}
    \item{tcrossprod}{\code{signature(x = "mpfrMatrix", y = "mpfrMatrix")}:
      Computes \eqn{xy'}, i.e., \code{x \%*\% t(y)}, typically more efficiently.}
    \item{tcrossprod}{\code{signature(x = "mpfrMatrix", y = "mpfr")}: ... }
    \item{tcrossprod}{\code{signature(x = "mpfr", y = "mpfrMatrix")}: ... }

    \item{coerce}{\code{signature(from = "mpfrArray", to = "array")}:
      coerces \code{from} to a \emph{numeric} array of the same dimension.}
    \item{coerce}{\code{signature(from = "mpfrArray", to = "vector")}:
      as for standard \code{\link{array}}s, this \dQuote{drops} the
      \code{dim} (and \code{dimnames}), i.e., returns an
      \code{\linkS4class{mpfr}} vector.}

    \item{Compare}{\code{signature(e1 = "mpfr", e2 = "mpfrArray")}: ... }
    \item{Compare}{\code{signature(e1 = "numeric", e2 = "mpfrArray")}: ... }
    \item{Compare}{\code{signature(e1 = "mpfrArray", e2 = "mpfr")}: ... }
    \item{Compare}{\code{signature(e1 = "mpfrArray", e2 = "numeric")}: ... }
    \item{dim}{\code{signature(x = "mpfrArray")}: ... }
    \item{dimnames<-}{\code{signature(x = "mpfrArray")}: ... }
    \item{dimnames}{\code{signature(x = "mpfrArray")}: ... }
    \item{show}{\code{signature(object = "mpfrArray")}: ... }
    \item{sign}{\code{signature(x = "mpfrArray")}: ... }

    \item{t}{\code{signature(x = "mpfrMatrix")}: tranpose the mpfrMatrix.}
    \item{aperm}{\code{signature(a = "mpfrArray")}: \code{aperm(a,
	perm)} is a generalization of \code{t(.)} to \emph{perm}ute the
      dimensions of an mpfrArray; it has the same semantics as the
      standard \code{\link{aperm}()} method for simple \R \code{\link{array}}s.}
  }
}
\author{Martin Maechler}
\seealso{
 \code{\link{mpfrArray}}, also for more examples.
}
\examples{
showClass("mpfrMatrix")

validObject(mm <- new("mpfrMatrix"))
validObject(aa <- new("mpfrArray"))

v6 <- mpfr(1:6, 128)
m6 <- new("mpfrMatrix", v6, Dim = c(2L, 3L))
validObject(m6)
m6
which(m6 == 3, arr.ind = TRUE) # |-->  (1, 2)
## Coercion back to "vector": Both of these work:
stopifnot(identical(as(m6, "mpfr"), v6),
          identical(as.vector(m6), v6)) # < but this is a "coincidence"

S2 <- m6[,-3] # 2 x 2
S3 <- rbind(m6, c(1:2,10)) ; s3 <- asNumeric(S3)
det(S2)
str(determinant(S2))
det(S3)
stopifnot(all.equal(det(S2), det(asNumeric(S2)), tol=1e-15),
	  all.equal(det(S3), det(s3), tol=1e-15))

## 2-column matrix indexing and replacement:
(sS <- S3[i2 <- cbind(1:2, 2:3)])
stopifnot(identical(asNumeric(sS), s3[i2]))
C3 <- S3; c3 <- s3
C3[i2] <- 10:11
c3[i2] <- 10:11
stopifnot(identical(asNumeric(C3), c3))

AA <- new("mpfrArray", as.vector(cbind(S3, -S3)), Dim=c(3L,3:2))
stopifnot(identical(AA[,,1] , S3), identical(AA[,,2] , -S3))
aa <- asNumeric(AA)

i3 <- cbind(3:1, 1:3, c(2L, 1:2))
ii3 <- Rmpfr:::.mat2ind(i3, dim(AA), dimnames(AA))
stopifnot(aa[i3] == new("mpfr", getD(AA)[ii3]))
stopifnot(identical(aa[i3], asNumeric(AA[i3])))
CA <- AA; ca <- aa
ca[i3] <- ca[i3] ^ 3
CA[i3] <- CA[i3] ^ 3
}
\keyword{classes}
