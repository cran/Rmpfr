\name{bind-methods}
\title{"mpfr" '...' - Methods for Functions cbind(), rbind()}
\docType{methods}
\alias{cbind}
\alias{rbind}
\alias{cbind-methods}
\alias{rbind-methods}
\alias{cbind,ANY-method}
\alias{cbind,Mnumber-method}
\alias{rbind,ANY-method}
\alias{rbind,Mnumber-method}
\description{
  \code{\link{cbind}} and \code{\link{rbind}} methods for signature
  \code{...} (see \code{\link[methods]{dotsMethods}}
  are provided for class \code{\linkS4class{Mnumber}}, i.e., for binding
  numeric vectors and class \code{"\linkS4class{mpfr}"} vectors and
  matrices (\code{"\linkS4class{mpfrMatrix}"}) together.
}
\usage{
cbind(\dots, deparse.level = 1)
rbind(\dots, deparse.level = 1)
}
\arguments{
  \item{\dots}{matrix-/vector-like \R objects to be bound together, see
    the \pkg{base} documentation, \code{\link[base:cbind]{cbind}}.
  }
  \item{deparse.level}{integer determining under which circumstances
    column and row names are built from the actual arguments'
    \sQuote{expression}, see \code{\link{cbind}}.}
}
\section{Methods}{
  \describe{
    \item{... = "Mnumber"}{is used to (c|r)bind multiprecision \dQuote{numbers}
      (inheriting from class \code{"\linkS4class{mpfr}"})
      together, maybe combined with simple numeric vectors.}
    \item{... = "ANY"}{reverts to \code{\link[base]{cbind}} and \code{rbind}
      from package \pkg{base}.}
  }
}
\value{
  typically a \sQuote{matrix-like} object, here typically of
  \code{\link{class} "\linkS4class{mpfrMatrix}"}.
}
% \references{ ~put references to the literature/web site here ~ }
\author{Martin Maechler}
\seealso{\code{\link{cbind2}}, \code{\link{cbind}},
  \code{\link[methods]{Methods}}.
}
\examples{
  cbind(1, mpfr(6:3, 70)/7, 3:0)
}
\keyword{methods}

