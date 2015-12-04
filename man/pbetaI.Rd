\name{pbetaI}
\alias{pbetaI}
\title{Accurate Incomplete Beta / Beta Probabilities For Integer Shapes}
\description{
  For integers \eqn{a}, \eqn{b}, \eqn{I_x(a,b)}{I(x; a,b)} aka
  \code{pbeta(x, a,b)} is a polynomial in x with rational coefficients,
  and hence arbitarily accurately computable.
}
\usage{
pbetaI(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE,
       precBits = NULL, rnd.mode = c("N","D","U","Z","A"))
}
\arguments{
  \item{q}{called \eqn{x}, above; vector of quantiles, in \eqn{[0,1]}.}
  \item{shape1, shape2}{the positive Beta \dQuote{shape} parameters,
    called \eqn{a, b}, above. \bold{Must} be integer valued for this
    function.}
  \item{ncp}{unused, only for compatibility with \code{\link{pbeta}},
    must be kept at its default, 0.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{precBits}{the precision (in number of bits) to be used in
    \code{\link{sumBinomMpfr}()}.}
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see \code{\link{mpfr}}.}
}
\value{
  an \code{"\linkS4class{mpfr}"} vector of the same length as \code{q}.
}
\note{For upper tail probabilities, i.e., when \code{lower.tail=FALSE},
  we may need large \code{precBits}, because the implicit or explicit
  \eqn{1 - P} computation suffers from severe cancellation.
  %% FIXME: Indeed, we may want to
  %% 1) determine a smart default  precBits        there
  %% 2) produce a  warning  when  cancellation is in the same order as precBits
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{Martin Maechler}
\seealso{
  \code{\link{pbeta}}, \code{\link{sumBinomMpfr}} \code{\link[gmp]{chooseZ}}.
}
\examples{
x <- (0:12)/16 # not all the way up ..
a <- 7; b <- 788

p.  <- pbetaI(x, a, b) ## still slow: %% TOO slow -- FIXME
pp  <- pbetaI(x, a, b, precBits = 2048)
## Currently, the lower.tail=FALSE  are computed "badly":
lp  <- log(pp)    ## = pbetaI(x, a, b, log.p=TRUE)
lIp <- log1p(-pp) ## = pbetaI(x, a, b, lower.tail=FALSE, log.p=TRUE)
 Ip <- 1 - pp     ## = pbetaI(x, a, b, lower.tail=FALSE)

if(Rmpfr:::doExtras()) { ## somewhat slow
   stopifnot(
     all.equal(lp,  pbetaI(x, a, b, precBits = 2048, log.p=TRUE)),
     all.equal(lIp, pbetaI(x, a, b, precBits = 2048, lower.tail=FALSE, log.p=TRUE),
               tol = 1e-230),
     all.equal( Ip, pbetaI(x, a, b, precBits = 2048, lower.tail=FALSE))
   )
}

rErr <- function(approx, true, eps = 1e-200) {
    true <- as.numeric(true) # for "mpfr"
    ifelse(Mod(true) >= eps,
           ## relative error, catching '-Inf' etc :
	   ifelse(true == approx, 0, 1 - approx / true),
           ## else: absolute error (e.g. when true=0)
	   true - approx)
}

rErr(pbeta(x, a, b), pp)
rErr(pbeta(x, a, b, lower=FALSE), Ip)
rErr(pbeta(x, a, b, log = TRUE),  lp)
rErr(pbeta(x, a, b, lower=FALSE, log = TRUE),  lIp)

a.EQ <- function(..., tol=1e-15) all.equal(..., tolerance=tol)
stopifnot(
  a.EQ(pp,  pbeta(x, a, b)),
  a.EQ(lp,  pbeta(x, a, b, log.p=TRUE)),
  a.EQ(lIp, pbeta(x, a, b, lower.tail=FALSE, log.p=TRUE)),
  a.EQ( Ip, pbeta(x, a, b, lower.tail=FALSE))
 )
}
\keyword{arith}
\keyword{distribution}
