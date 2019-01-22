\name{sapplyMpfr}
\alias{sapplyMpfr}
\title{Apply a Function over a "mpfr" Vector}
\description{
  Users may be disappointed to note that \code{\link{sapply}()} or
  \code{\link{vapply}()} typically do not work with \code{"mpfr"}
  numbers.

  This is a simple (but strong) approach to work around the problem,
  based on \code{\link{lapply}()}.

  Note that this is not yet as flexible as \code{\link{sapply}()} for
  atomic vectors.
}
\usage{
sapplyMpfr(X, FUN, \dots)
}
\arguments{
  \item{X}{a vector, possibly of class \code{"\linkS4class{mpfr}"}.}
  \item{FUN}{a function returning an \code{"\linkS4class{mpfr}"}
    number. (TODO: A function returning a \emph{vector} of mpfr numbers
    or even \code{"\linkS4class{mpfrArray}"}).}
  \item{\dots}{further arguments passed to \code{\link{lapply}},
    typically further arguments to \code{FUN}.}
}
%% \details{
%% }
\value{
  an  \code{"\linkS4class{mpfr}"} vector, typically of the same length
  as \code{X}.
}
\author{Martin Maechler}
\note{
  Another workaround could be to use \preformatted{
    res <- lapply(....)
    sapply(res, asNumeric)
  }
}
%% \seealso{
%%   \code{\link{..}},
%% }
\examples{
## The function is simply defined as
function (X, FUN, ...)
  new("mpfr", unlist(lapply(X, FUN, ...), recursive = FALSE))

if(require("Bessel")) # here X, is simple
  bImp <- sapplyMpfr(0:4, function(k)
                     besselI.nuAsym(mpfr(1.31e9, 128), 10, expon.scaled=TRUE, k.max=k))
}
\keyword{manip}
