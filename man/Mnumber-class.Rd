\name{Mnumber-class}
\Rdversion{1.1}
\docType{class}
\alias{Mnumber-class}
\title{Class "Mnumber" of "mpfr" and regular numbers and arrays from them}
\description{
  Class \code{"Mnumber"} is a class union of \code{"\linkS4class{mpfr}"}
  and regular numbers and arrays from them.\cr
  Its purpose is for method dispatch, notably defining a
  \code{cbind(...)} method where \code{...} contains objects of one of
  the member classes of \code{"Mnumber"}.
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "mpfrMatrix", y = "Mnumber")}: ... }
    \item{crossprod}{\code{signature(x = "mpfr", y = "Mnumber")}: ... }
    \item{tcrossprod}{\code{signature(x = "Mnumber", y = "mpfr")}: ...}
  }
  etc.  These are documented with the classes \code{\linkS4class{mpfr}}
  and or \code{\linkS4class{mpfrMatrix}}.
}
\seealso{
  the \code{\linkS4class{array_or_vector}} sub class;
  \code{\link{cbind-methods}}.
}
\examples{
## "Mnumber" encompasses (i.e., "extends") quite a few
##  "vector / array - like" classes:
showClass("Mnumber")
stopifnot(extends("mpfrMatrix", "Mnumber"),
          extends("array",      "Mnumber"))
}
\keyword{classes}
