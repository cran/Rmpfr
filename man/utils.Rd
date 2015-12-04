\name{mpfr.utils}
\title{MPFR Number Utilities}
\alias{c.mpfr}
\alias{diff.mpfr}
\alias{mpfrIs0}
\alias{.mpfr.is.whole}
\alias{mpfrVersion}
%% next two are deprecated since Aug.2015 (for Rmpfr 0.6.0):
\alias{mpfr.is.0}
\alias{mpfr.is.integer}
%
\description{
  \code{mpfrVersion()} returns the version of the MPFR library which
  \pkg{Rmpfr} is currently linked to.

  \code{\link{c}(x,y,...)} can be used to combine MPFR numbers in the
  same way as regular numbers \bold{IFF} the first argument \code{x} is
  of class \code{\linkS4class{mpfr}}.

  \code{mpfrIs0(.)} uses the MPFR library in the documented way to
  check if (a vector of) MPFR numbers are zero.  It was called
  \code{mpfr.is.0} which is strongly deprecated now.

  \code{.mpfr.is.whole(x)} uses the MPFR library in the documented way to
  check if (a vector of) MPFR numbers is integer \emph{valued}.  This is
  equivalent to \code{x == round(x)}, but \emph{not} at all to
  \code{is.integer(as(x, "numeric"))}.\cr
  You should typically rather use (the \code{"mpfr"} method of the
  generic function) \code{\link{is.whole}(x)} instead.
  The former name \code{mpfr.is.integer} is deprecated now.
}
\usage{
mpfrVersion()
mpfrIs0(x)
%% .mpfr.is.whole(x)
\method{c}{mpfr}(\dots)
\method{diff}{mpfr}(x, lag = 1L, differences = 1L, \dots)
}
\arguments{
  \item{x}{an object of class \code{\linkS4class{mpfr}}.}
  \item{\dots}{for \code{diff}, further \code{\linkS4class{mpfr}} class
    objects or simple numbers (\code{\link{numeric}} vectors) which are coerced to
    \code{mpfr} with default precision of 128 bits.}
  \item{lag, differences}{for \code{diff()}: exact same meaning as in
    \code{\link{diff}()}'s default method, \code{\link{diff.default}}.}
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\value{
  \code{mpfrIs0} returns a logical vector of length \code{length(x)}
  with values \code{TRUE} iff the corresponding \code{x[i]} is an MPFR
  representation of zero (\code{0}).
  \cr Similarly, \code{.mpfr.is.whole} and \code{is.whole} return a
  logical vector of length \code{length(x)}.

  \code{mpfrVersion} returns an object of S3 class
  \code{"\link{numeric_version}"}, so it can be used in comparisons.

  The other functions return MPFR number (vectors), i.e., extending
  class \code{\linkS4class{mpfr}}.
}
\seealso{
  \code{\link{str.mpfr}} for the \code{\link{str}} method.
  \code{\link{erf}} for special mathematical functions on MPFR.

  The class description \code{\linkS4class{mpfr}} page mentions many
  generic arithmetic and mathematical functions for which \code{"mpfr"}
  methods are available.
}
\examples{
mpfrVersion()

(x <- c(Const("pi", 64), mpfr(-2:2, 64)))
mpfrIs0(x)    # one of them is
x[mpfrIs0(x)] # but it may not have been obvious..
str(x)

x <- rep(-2:2, 5)
stopifnot(is.whole(mpfr(2, 500) ^ (1:200)),
          all.equal(diff(x), diff(as.numeric(x))))
}
\keyword{arith}
