\name{pmax}
\Rdversion{1.1}
\alias{pmax}
\alias{pmin}
\alias{pmax-methods}
\alias{pmin-methods}
\alias{pmax,ANY-method}
\alias{pmax,mNumber-method}
\alias{pmin,ANY-method}
\alias{pmin,mNumber-method}
\title{Parallel Maxima and Minima}
\description{
  Returns the parallel maxima and minima of the input values.

  The functions \code{pmin} and \code{pmax} have been made S4 generics,
  and this page documents the \dQuote{\code{...} method for class
    \code{"mNumber"}}, i.e., for arguments that are numeric or from
  \code{\link{class} "\linkS4class{mpfr}"}.
}
\usage{
pmax(\dots, na.rm = FALSE)
pmin(\dots, na.rm = FALSE)
}
\arguments{
  \item{\dots}{numeric or arbitrary precision numbers (class
  \code{\linkS4class{mpfr}}).}
  \item{na.rm}{a logical indicating whether missing values should be
    removed.}
}
\details{
  See \code{\link[base:Extremes]{pmax}}, the documentation of the base
  functions, i.e., default methods.
}
\value{
  vector-like, of length the longest of the input vectors; typically of
  class \code{\linkS4class{mpfr}}, for the methods here.
}
\section{Methods}{
  \describe{
    \item{... = "ANY"}{the default method, really just
      \code{\link[base:Extremes]{base::pmin}} or \code{base::pmax},
      respectively.}
    \item{... = "mNumber"}{the method for \code{\linkS4class{mpfr}}
      arguments, mixed with numbers; designed to follow the same semantic as
      the default method.}
  }
}
\seealso{
  The documentation of the \pkg{base} functions,
  \code{\link[base:Extremes]{pmin}} and \code{pmax}; also
  \code{\link{min}} and \code{max}; further,

  \code{\link{range}} (\emph{both} min and max).
}
\examples{
 (pm <- pmin(1.35, mpfr(0:10, 77)))
 stopifnot(pm == pmin(1.35, 0:10))
}
\keyword{methods}
\keyword{univar}
\keyword{arith}
