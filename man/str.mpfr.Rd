\name{str.mpfr}
\alias{str.mpfr}
\title{Compactly Show STRucture of Rmpfr Number Object}
\description{
  The \code{\link{str}} method for objects of class
  \code{\linkS4class{mpfr}} produces a bit more useful output than
  the default method \code{\link{str.default}}.
}
\usage{
\method{str}{mpfr}(object, nest.lev, give.head=TRUE, \dots)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{mpfr}}.}
  \item{nest.lev}{for \code{\link{str}()}, typically only used when
    called by a higher level \code{str()}.}
  \item{give.head}{logical indicating if the \dQuote{header} should be printed.}
  \item{\dots}{further arguments passed to \code{\link{str}}.}
}
\seealso{
  \code{\link{mpfr.is.0}} for many more utilities.
}
\examples{
(x <- c(Const("pi", 64), mpfr(-2:2, 64)))
str(x)
str(list(pi = pi, x.mpfr = x))
}
\keyword{utilities}
