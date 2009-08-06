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
\seealso{
  \code{\link{cbind-methods}}
}
\examples{
showClass("Mnumber")
}
\keyword{classes}
