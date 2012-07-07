\name{Mnumber-class}
\Rdversion{1.1}
\docType{class}
\alias{Mnumber-class}
\alias{mNumber-class}
\alias{numericVector-class}
\title{Class "Mnumber" and "mNumber" of "mpfr" and regular numbers and arrays from them}
\description{
  Classes \code{"Mnumber"}  \code{"mNumber"} are class unions of \code{"\linkS4class{mpfr}"}
  and regular numbers and arrays from them.\cr
  Its purpose is for method dispatch, notably defining a
  \code{cbind(...)} method where \code{...} contains objects of one of
  the member classes of \code{"Mnumber"}.

  Classes \code{"mNumber"} is considerably smaller is it does \emph{not}
  contain \code{"matrix"} and \code{"array"} since these also extend \code{"character"}
  which is not really desirable for generalized numbers.
  It extends the simple \code{"numericVector"} class by \code{mpfr*} classes.
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

Mnsub <- names(getClass("Mnumber")@subclasses)
(mNsub <- names(getClass("mNumber")@subclasses))
## mNumber has *one* subclass which is not in Mnumber:
setdiff(mNsub, Mnsub)# namely "numericVector"
## The following are only subclasses of "Mnumber", but not of "mNumber":
setdiff(Mnsub, mNsub)
}
\keyword{classes}
