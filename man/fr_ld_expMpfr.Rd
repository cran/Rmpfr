\name{frexpMpfr}
\alias{frexpMpfr}
\alias{ldexpMpfr}
\title{Base-2 Representation and Multiplication of Mpfr Numbers}% ../DPQ/man/fr_ld_exp.Rd
\description{
  MPFR - versions of the C99 (and POSIX) standard C (and C++) mathlib
  functions \code{frexp()} and \code{ldexp()}.

  \code{frexpMpfr(x)} computes base-2 exponent \code{e} and \dQuote{mantissa},
  or \emph{fraction} \code{r}, such that \eqn{x = r * 2^e}, where \eqn{r \in
    [0.5, 1)} (unless when \code{x} is in \code{c(0, -Inf, Inf, NaN)}
  where \code{r == x} and \code{e} is 0),
  and \eqn{e} is integer valued.

  \code{ldexpMpfr(f, E)} is the \emph{inverse} of \code{frexpMpfr()}: Given
  fraction or mantissa \code{f} and integer exponent \code{E}, it returns
  \eqn{x = f * 2^E}.
  Viewed differently, it's the fastest way to multiply or divide MPFR
  numbers with \eqn{2^E}.
}
\usage{
frexpMpfr(x,    rnd.mode = c("N", "D", "U", "Z", "A"))
ldexpMpfr(f, E, rnd.mode = c("N", "D", "U", "Z", "A"))
}
\arguments{
  \item{x}{numeric (coerced to \code{double}) vector.}
  \item{f}{numeric fraction (vector), in \eqn{[0.5, 1)}.}
  \item{E}{integer valued, exponent of \code{2}, i.e., typically in
    \code{(-1024-50):1024}, otherwise the result will underflow to 0 or
    overflow to \code{+/- Inf}.}
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see \code{\link{mpfr}}.}
}
%% \details{

%% }
\value{
  \code{frexpMpfr} returns a \code{\link{list}} with named components \code{r}
  (of class \code{mpfr}) and \code{e} (integer valued, of type
  \code{integer} is small enough, \code{"double"} otherwise).
}
\references{
  On unix-alikes, typically \command{man frexp} and \command{man ldexp}
}
\author{Martin Maechler}
%% \note{

%% }
\seealso{
  Somewhat related, \code{\link{.mpfr2exp}()}.
  %_TODO_: that DPQ version is not yet on CRAN:
  %\code{\link[DPQ]{frexp}()} and \code{ldexp()} in package \CRANpkg{DPQ}.
}
\examples{
set.seed(47)
x <- c(0, 2^(-3:3), (-1:1)/0,
       sort(rlnorm(2^12, 10, 20) * sample(c(-1,1), 512, replace=TRUE)))
head(xM <- mpfr(x, 128), 11)
str(rFM <- frexpMpfr(xM))
d.fr <- with(rFM, data.frame(x=x, r=asNumeric(r), e=e))
head(d.fr , 16)
tail(d.fr)
ar <- abs(rFM$r)
stopifnot(0.5 <= ar[is.finite(x) & x != 0], ar[is.finite(x)] < 1,
          is.integer(rFM$e))
ldx <- mpfr(with(rFM, ldexpMpfr(r, e)))
(iN <- which(is.na(x))) # 10
stopifnot(exprs = {
  all.equal(xM, ldx, tol = 2^-124) # allow 4 bits loss, but apart from the NA, even:
  identical(xM[-iN], ldx[-iN])
  is.na(xM [iN])
  is.na(ldx[iN])
})
}
\keyword{arithmetic}
\keyword{utilities}
