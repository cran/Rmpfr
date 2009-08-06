\name{array_or_vector-class}
\docType{class}
\alias{array_or_vector-class}

\title{Auxiliary Class "array\_or\_vector"}
\description{
  \code{"array_or_vector"} is the class union of
  \code{c("array", "matrix", "vector")} and exists for its use in
  signatures of method definitions.
}
\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "mpfrMatrix", y = "array_or_vector")}: ... }
    \item{crossprod}{\code{signature(x = "mpfr", y = "array_or_vector")}: ... }
    \item{tcrossprod}{\code{signature(x = "array_or_vector", y =
	"mpfr")}: ...}
  }
  etc.  These are documented with the classes \code{\linkS4class{mpfr}}
  and or \code{\linkS4class{mpfrMatrix}}.
}
\details{
  Using \code{"array_or_vector"} instead of just \code{"vector"} in a
  signature makes an important difference: E.g., if we had
  \code{setMethod(crossprod, c(x="mpfr", y="vector"), function(x,y) CPR(x,y))},
  a call \code{crossprod(x, matrix(1:6, 2,3))} would extend into a call
  of \code{CPR(x, as(y, "vector"))} such that \code{CPR()}'s second
  argument would simply be a vector instead of the desired
  \eqn{2\times3}{2 x 3} matrix.
}
%\author{Martin Maechler}
% \examples{
% showClass("array_or_vector")
% }
\keyword{classes}
