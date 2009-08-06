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
\alias{coerce,mpfrArray,array-method}
\alias{coerce,mpfrMatrix,matrix-method}
%
\alias{as.vector,mpfrArray,missing-method}
%
\alias{[<-,mpfrArray,ANY,ANY,ANY-method}
\alias{[<-,mpfrArray,ANY,ANY,mpfr-method}
\alias{[<-,mpfrArray,ANY,missing,mpfr-method}
\alias{[<-,mpfrArray,missing,ANY,mpfr-method}
%
\alias{[,mpfrArray,ANY,ANY,ANY-method}
\alias{[,mpfrArray,ANY,missing,missing-method}

\alias{\%*\%,mpfrMatrix,mpfrMatrix-method}
\alias{\%*\%,mpfrMatrix,array_or_vector-method}
\alias{\%*\%,array_or_vector,mpfrMatrix-method}
\alias{crossprod,mpfrMatrix,mpfrMatrix-method}
\alias{tcrossprod,mpfrMatrix,mpfrMatrix-method}

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
  Objects can be created by calls of the form
  \code{new("mpfrMatrix", ...)} or \code{new("mpfrArray", ...)}, or also
  by \code{dim(x) <- dd} or \code{t(x)} where \code{x} is a an
  \code{\linkS4class{mpfr}} \dQuote{number vector}.
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
    \item{coerce}{\code{signature(from = "mpfrArray", to = "array")}: ... }
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
% \note{ ~~further notes~~ }
% \seealso{
% 	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
% 	or \code{\linkS4class{CLASSNAME}} for links to other classes
% }
\examples{
showClass("mpfrMatrix")

validObject(mm <- new("mpfrMatrix"))
validObject(aa <- new("mpfrArray"))

m6 <- new("mpfrMatrix", mpfr(1:6, 128), Dim = c(2L, 3L))
validObject(m6)
m6
which(m6 == 3, arr.ind = TRUE) # |-->  (1, 2)
}
\keyword{classes}
