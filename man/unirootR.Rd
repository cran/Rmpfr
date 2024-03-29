\name{unirootR}
\alias{unirootR}
\title{One Dimensional Root (Zero) Finding -- in pure \R}
\usage{
unirootR(f, interval, \dots,
        lower = min(interval), upper = max(interval),
        f.lower = f(lower, ...), f.upper = f(upper, ...),
        extendInt = c("no", "yes", "downX", "upX"),
        trace = 0, verbose = as.logical(trace),
        verbDigits = max(3, min(20, -log10(tol)/2)),
        tol = .Machine$double.eps^0.25, maxiter = 1000L,
        check.conv = FALSE,
        warn.no.convergence = !check.conv,
        epsC = NULL)
}
\arguments{% base R : ~/R/D/r-devel/R/src/library/stats/man/uniroot.Rd
  \item{f}{the function for which the root is sought.}
  \item{interval}{a vector containing the end-points of the interval
    to be searched for the root.}
  \item{\dots}{additional named or unnamed arguments to be passed
    to \code{f}}
  \item{lower, upper}{the lower and upper end points of the interval to
    be searched.}
  \item{f.lower, f.upper}{the same as \code{f(upper)} and
    \code{f(lower)}, respectively.  Passing these values from the caller
    where they are often known is more economical as soon as \code{f()}
    contains non-trivial computations.}
  \item{extendInt}{character string specifying if the interval
    \code{c(lower,upper)} should be extended or directly produce an error
    when \code{f()} does not have differing signs at the endpoints.  The
    default, \code{"no"}, keeps the search interval and hence produces
    an error.  Can be abbreviated.}
  \item{trace}{integer number; if positive, tracing information is
    produced.  Higher values giving more details.}
  \item{verbose}{logical (or integer) indicating if (and how much)
    verbose output should be produced during the iterations.}
  \item{verbDigits}{used only if \code{verbose} is true, indicates the
    number of digits numbers should be printed with, using
    \code{\link{format}(., digits=verbDigits)}.}
  \item{tol}{the desired accuracy (convergence tolerance).}
  \item{maxiter}{the maximum number of iterations.}
  \item{check.conv}{logical indicating whether non convergence
    should be caught as an error, notably non-convergence in \code{maxiter}
    iterations should be an error instead of a warning.}
  \item{warn.no.convergence}{if set to \code{FALSE} there's no warning
    about non-convergence.  Useful to just run a few iterations.}
  \item{epsC}{positive number or \code{NULL} in which case a smart
    default is sought.  This should specify the \dQuote{achievable
      machine precision} \emph{for} the given numbers and their
    arithmetic.

    The default will set this to \code{\link{.Machine}$double.eps} for
    double precision numbers, and will basically use
    \code{2 ^ - min(getPrec(f.lower), getPrec(f.upper))} when that works
    (as, e.g., for \code{\linkS4class{mpfr}}-numbers) otherwise.

    This is factually a lower bound for the achievable lower bound, and
    hence, setting \code{tol} smaller than \code{epsC} is typically
    non-sensical and produces a warning.
  }
}
\description{
  The function \code{unirootR} searches the interval from \code{lower}
  to \code{upper} for a root (i.e., zero) of the function \code{f} with
  respect to its first argument.

  \code{unirootR()} is \dQuote{clone} of \code{\link{uniroot}()},
  written entirely in \R, in a way that it works with
  \code{\linkS4class{mpfr}}-numbers as well.
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  Either \code{interval} or both \code{lower} and \code{upper} must be
  specified: the upper endpoint must be strictly larger than the lower
  endpoint.  The function values at the endpoints must be of opposite
  signs (or zero), for \code{extendInt="no"}, the default.  Otherwise, if
  \code{extendInt="yes"}, the interval is extended on both sides, in
  search of a sign change, i.e., until the search interval \eqn{[l,u]}
  satisfies \eqn{f(l) \cdot f(u) \le 0}{f(l) * f(u) <= 0}.

  If it is \emph{known how} \eqn{f} changes sign at the root
  \eqn{x_0}{x0}, that is, if the function is increasing or decreasing there,
  \code{extendInt} can (and typically should) be specified as
  \code{"upX"} (for \dQuote{upward crossing}) or \code{"downX"},
  respectively.  Equivalently, define \eqn{S := \pm 1}{S:= +/- 1}, to
  require \eqn{S = \mathrm{sign}(f(x_0 + \epsilon))}{S = sign(f(x0 +
    eps))} at the solution.  In that case, the search interval \eqn{[l,u]}
  possibly is extended to be such that \eqn{S\cdot f(l)\le 0}{%
    S * f(l) <= 0} and \eqn{S \cdot f(u) \ge 0}{S * f(u) >= 0}.


  The function only uses \R code with basic arithmetic, such that it
  should also work with \dQuote{generalized} numbers (such as
  \code{\linkS4class{mpfr}}-numbers) as long the necessary
  \code{\link{Ops}} methods are defined for those.

  The underlying algorithm assumes a continuous function (which then is
  known to have at least one root in the interval).

  Convergence is declared either if \code{f(x) == 0} or the change in
  \code{x} for one step of the algorithm is less than \code{tol} (plus an
  allowance for representation error in \code{x}).

  If the algorithm does not converge in \code{maxiter} steps, a warning
  is printed and the current approximation is returned.

  \code{f} will be called as \code{f(\var{x}, ...)} for a (generalized)
  numeric value of \var{x}.
}
\value{
  A list with four components: \code{root} and \code{f.root} give the
  location of the root and the value of the function evaluated at that
  point. \code{iter} and \code{estim.prec} give the number of iterations
  used and an approximate estimated precision for \code{root}.  (If the
  root occurs at one of the endpoints, the estimated precision is
  \code{NA}.)
}
\source{
  Based on \code{zeroin()} (in package \pkg{rootoned}) by John Nash who
  manually translated the C code in \R's \code{zeroin.c} and on
  \code{\link{uniroot}()} in \R's sources.
}
\references{
  Brent, R. (1973),   see \code{\link{uniroot}}.
}
\seealso{
  \R's own (\pkg{stats} package) \code{\link{uniroot}}.
  \code{\link{polyroot}} for all complex roots of a polynomial;
  \code{\link{optimize}}, \code{\link{nlm}}.
}
\examples{
require(utils) # for str

## some platforms hit zero exactly on the first step:
## if so the estimated precision is 2/3.
f <- function (x,a) x - a
str(xmin <- unirootR(f, c(0, 1), tol = 0.0001, a = 1/3))

## handheld calculator example: fixpoint of cos(.):
rc <- unirootR(function(x) cos(x) - x, lower=-pi, upper=pi, tol = 1e-9)
rc$root

## the same with much higher precision:
rcM <- unirootR(function(x) cos(x) - x,
                 interval= mpfr(c(-3,3), 300), tol = 1e-40)
rcM
x0 <- rcM$root
stopifnot(all.equal(cos(x0), x0,
                    tol = 1e-40))## 40 digits accurate!

str(unirootR(function(x) x*(x^2-1) + .5, lower = -2, upper = 2,
            tol = 0.0001), digits.d = 10)
str(unirootR(function(x) x*(x^2-1) + .5, lower = -2, upper = 2,
            tol = 1e-10 ), digits.d = 10)

## A sign change of f(.), but not a zero but rather a "pole":
tan. <- function(x) tan(x * (Const("pi",200)/180))# == tan( <angle> )
(rtan <- unirootR(tan., interval = mpfr(c(80,100), 200), tol = 1e-40))
## finds 90 {"ok"}, and now gives a warning
\dontshow{stopifnot(all.equal(rtan$root, 90, tolerance = 1e-38))}

## Find the smallest value x for which exp(x) > 0 (numerically):
r <- unirootR(function(x) 1e80*exp(x)-1e-300, c(-1000,0), tol = 1e-15)
str(r, digits.d = 15) ##> around -745, depending on the platform.

exp(r$root)     # = 0, but not for r$root * 0.999...
minexp <- r$root * (1 - 10*.Machine$double.eps)
exp(minexp)     # typically denormalized

## --- using mpfr-numbers :

## Find the smallest value x for which exp(x) > 0 ("numerically");
## Note that mpfr-numbers underflow *MUCH* later than doubles:
## one of the smallest mpfr-numbers {see also ?mpfr-class } :
(ep.M <- mpfr(2, 55) ^ - ((2^30 + 1) * (1 - 1e-15)))
r <- unirootR(function(x) 1e99* exp(x) - ep.M, mpfr(c(-1e20, 0), 200))
r # 97 iterations; f.root is very similar to ep.M

## interval extension 'extendInt'  --------------

f1 <- function(x) (121 - x^2)/(x^2+1)
f2 <- function(x) exp(-x)*(x - 12)
tools::assertError(unirootR(f1, c(0,10)), verbose=TRUE)
##--> error: f() .. end points not of opposite sign

## where as  'extendInt="yes"'  simply first enlarges the search interval:
u1 <- unirootR(f1, c(0,10),extendInt="yes", trace=1)
u2 <- unirootR(f2, mpfr(c(0,2), 128), extendInt="yes", trace=2, verbose=FALSE, tol = 1e-25)
stopifnot(all.equal(u1$root, 11, tolerance = 1e-5),
          all.equal(u2$root, 12, tolerance = 1e-23))

## The *danger* of interval extension:
## No way to find a zero of a positive function, but
## numerically, f(-|M|) becomes zero :
u3 <- unirootR(exp, c(0,2), extendInt="yes", trace=TRUE)

## Nonsense example (must give an error):
tools::assertCondition( unirootR(function(x) 1, 0:1, extendInt="yes"),
                       "error", verbose=TRUE)
}
\keyword{optimize}
