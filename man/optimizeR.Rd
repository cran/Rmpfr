\name{optimizeR}
\alias{optimizeR}
\title{High Precisione One-Dimensional Optimization}
%% Compare with ~/R/D/r-devel/R/src/library/stats/man/optimize.Rd
\description{
  \code{optimizeR} searches the intervalfrom
  \code{lower} to \code{upper} for a minimum % TODO ? or maximum
  of the function \code{f} with respect to its first argument.
}
\usage{
optimizeR(f, lower, upper, ..., tol = 1e-20,
          method = c("Brent", "GoldenRatio"),
          maximum = FALSE,
          precFactor = 2.0, precBits = -log2(tol) * precFactor,
          maxiter = 1000, trace = FALSE)
}
\arguments{
  \item{f}{the function to be optimized. \code{f(x)} must work
    \dQuote{in \pkg{Rmpfr} arithmetic} for \code{optimizer()} to make sense.
    The function is either minimized or maximized over its first argument
    depending on the value of \code{maximum}.}
  \item{\dots}{additional named or unnamed arguments to be passed
    to \code{f}.}
  \item{lower}{the lower end point of the interval to be searched.}
  \item{upper}{the upper end point of the interval to be searched.}
  \item{tol}{the desired accuracy, typically higher than double
    precision, i.e., \code{tol < 2e-16}.}
  \item{method}{\code{\link{character}} string specifying the
    optimization method.}%
  \item{maximum}{logical indicating if \eqn{f()} should be maximized or
    minimized (the default).}
  \item{precFactor}{only for default \code{precBits} construction: a factor
    to multiply with the number of bits directly needed for \code{tol}.}
  \item{precBits}{number of bits to be used for
    \code{\link[Rmpfr]{mpfr}} numbers used internally.}
  \item{maxiter}{maximal number of iterations to be used.}
  \item{trace}{integer or logical indicating if and how iterations
    should be monitored; if an integer \eqn{k}, print every \eqn{k}-th
    iteration.}
}
\details{
  \describe{
    \item{\code{"Brent"}:}{Brent(1973)'s simple and robust algorithm
      is a hybrid, using a combination of the golden ratio and local
      quadratic (\dQuote{parabolic}) interpolation.  This is the same
      algorithm as standard \R's \code{\link{optimize}()}, adapted to
      high precision numbers.

      In smooth cases, the convergence is considerably faster than the golden
      section or Fibonacci ratio algorithms.
    }
    \item{\code{"GoldenRatio"}:}{The golden ratio method works as follows:
      from a given interval containing the solution, it constructs the
      next point in the golden ratio between the interval boundaries.
    }
  }% end{describe}
}
\value{
  A \code{\link{list}} with components \code{minimum} (or \code{maximum})
  and \code{objective} which give the location of the minimum (or maximum)
  and the value of the function at that point;
  \code{iter} specifiying the number of iterations, the logical
  \code{convergence} indicating if the iterations converged and
  \code{estim.prec} which is an estimate or an upper bound of the final
  precision (in \eqn{x}).
  \code{method} the string of the method used.
}
% \references{
% }
\author{
  \code{"GoldenRatio"} is based on Hans W Borchert's
  \code{\link[pracma]{golden_ratio}};
  modifications and \code{"Brent"} by Martin Maechler.
}
% \note{
% }
\seealso{
  \R's standard \code{\link{optimize}}; \pkg{Rmpfr}'s \code{\link{unirootR}}.
}
\examples{
iG5 <- function(x) -exp(-(x-5)^2/2)
curve(iG5, 0, 10, 200)
o.dp  <- optimize (iG5, c(0, 10)) #->  5 of course
oM.gs <- optimizeR(iG5, 0, 10, method="Golden")
oM.Br <- optimizeR(iG5, 0, 10, method="Brent", trace=TRUE)
oM.gs$min ; oM.gs$iter
oM.Br$min ; oM.Br$iter
(doExtras <- Rmpfr:::doExtras())
if(doExtras) {## more accuracy {takes a few seconds}
 oM.gs <- optimizeR(iG5, 0, 10, method="Golden", tol = 1e-70)
 oM.Br <- optimizeR(iG5, 0, 10,                  tol = 1e-70)
}
rbind(Golden = c(err = as.numeric(oM.gs$min -5), iter = oM.gs$iter),
      Brent  = c(err = as.numeric(oM.Br$min -5), iter = oM.Br$iter))

## ==> Brent is orders of magnitude more efficient !

## Testing on the sine curve with 40 correct digits:
sol <- optimizeR(sin, 2, 6, tol = 1e-40)
str(sol)
sol <- optimizeR(sin, 2, 6, tol = 1e-50,
                 precFactor = 3.0, trace = TRUE)
pi.. <- 2*sol$min/3
print(pi.., digits=51)
stopifnot(all.equal(pi.., Const("pi", 256), tol = 10*1e-50))

## a harder one:
f.sq <- function(x) sin(x-2)^4 + sqrt(pmax(0,(x-1)*(x-4)))*(x-2)^2
curve(f.sq, 0, 4.5, n=1000)
msq <- optimizeR(f.sq, 0, 5, tol = 1e-50, trace=5)
str(msq) # ok
stopifnot(abs(msq$minimum - 2) < 1e-49)

if(doExtras) { # considerably more expensive
## find the other local minimum: -- non-smooth ==> Golden-section is used
msq2 <- optimizeR(f.sq, 3.5, 5, tol = 1e-50, trace=10)
stopifnot(abs(msq2$minimum - 4) < 1e-49)
}

## and a local maximum:
msq3 <- optimizeR(f.sq, 3, 4, maximum=TRUE, trace=2)
stopifnot(abs(msq3$maximum - 3.57) < 1e-2)




##----- "impossible" one to get precisely ------------------------

ff <- function(x) exp(-1/(x-8)^2)
curve(exp(-1/(x-8)^2), -3, 13, n=1001)
(opt. <- optimizeR(function(x) exp(-1/(x-8)^2), -3, 13, trace = 5))
## -> close to 8 {but not very close!}
ff(opt.$minimum) # gives 0
if(doExtras) {
 ## try harder ... in vain ..
 str(opt1 <- optimizeR(ff, -3,13, tol = 1e-60, precFactor = 4))
 print(opt1$minimum, digits=20)
 ## still just  7.99998038 or 8.000036655 {depending on method}
}
}
\keyword{optimize}
