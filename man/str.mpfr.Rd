\name{str.mpfr}
\alias{str.mpfr}
\title{Compactly Show STRucture of Rmpfr Number Object}
\description{
  The \code{\link{str}} method for objects of class
  \code{\linkS4class{mpfr}} produces a bit more useful output than
  the default method \code{\link{str.default}}.
}
\usage{
\method{str}{mpfr}(object, nest.lev, internal = FALSE,
    give.head = TRUE, digits.d = 12, vec.len = NULL, drop0trailing=TRUE,
    width = getOption("width"), \dots)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{mpfr}}.}
  \item{nest.lev}{for \code{\link{str}()}, typically only used when
    called by a higher level \code{str()}.}
  \item{internal}{logical indicating if the low-level internal structure
    should be shown; if true (not by default), uses \code{str(object@.Data)}.}
  \item{give.head}{logical indicating if the \dQuote{header} should be
    printed.}
  \item{digits.d}{the number of digits to be used, will be passed
    \code{\link{formatMpfr}()} and hence \code{NULL} will use \dQuote{as
      many as needed}, i.e. often too many.  If this is a number, as per
    default, less digits will be used in case the precision
    (\code{\link{getPrec}(object)}) is smaller.}
  \item{vec.len}{the number of \emph{elements} that will be shown.  The
    default depends on the precision of \code{object} and \code{width}
    (since \pkg{Rmpfr} 0.6-0, it was \code{3} previously).}
  \item{drop0trailing}{logical, passed to \code{\link{formatMpfr}()} (with a
    different default here).}
  \item{width}{the (approximately) desired width of output, see
    \code{\link{options}(width = .)}.}
  \item{\dots}{further arguments, passed to \code{\link{formatMpfr}()}.}
}
\seealso{
  \code{\link{mpfrIs0}} for many more utilities.
}
\examples{
(x <- c(Const("pi", 64), mpfr(-2:2, 64)))
str(x)
str(list(pi = pi, x.mpfr = x))
str(x ^ 1000)
str(x ^ -1e4, digits=NULL) # full precision

str(x, internal = TRUE) # internal low-level (for experts)

uu <- Const("pi", 16)# unaccurate
str(uu) # not using default of 12 digits
}
\keyword{utilities}
