\name{qnormI}
\alias{qnormI}
\title{Gaussian / Normal Quantiles \code{qnorm()} via Inversion}
\description{
  Compute Gaussian or Normal Quantiles \code{\link{qnorm}(p, *)} via
  inversion of our \dQuote{mpfr-ified} arbitrary accurate
  \code{\link[Rmpfr]{pnorm}()}, using our \code{\link{unirootR}()} root
  finder.
}
\usage{%--> ../R/unirootR.R
qnormI(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE,
       trace = 0, verbose = as.logical(trace),
       tol,
       useMpfr = any(prec > 53),
       give.full = FALSE,
       \dots)
}
\arguments{% ~/R/D/r-devel/R/src/library/stats/man/Normal.Rd <<<<<
  \item{p}{vector of probabilities.}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.}
  %%--
  \item{trace}{integer passed to \code{\link{unirootR}()}.  If positive,
    information about a search interval extension will be printed to the console.}
  \item{verbose}{logical indicating if progress details should be printed
    to the console.}
  \item{tol}{optionally the desired accuracy (convergence tolerance); if
    missing or not finite, it is computed as \eqn{2^-{pr+2}} where the
    precision \eqn{pr} is basically \code{max(\link{getPrec}(p+mean+sd))}.}
  \item{useMpfr}{logical indicating if \code{\link{mpfr}} arithmetic should
    be used.}%  Its default is derived from \code{tol} when that is specified.}
  \item{give.full}{logical indicating if the \emph{full} result of
    \code{\link{unirootR}()} should be returned (when applicable).}
  \item{\dots}{optional further arguments passed to \code{\link{unirootR}()}
    such as \code{maxiter}, \code{verbDigits}, \code{check.conv},
    \code{warn.no.convergence}, and \code{epsC}.}
}
\value{
  If \code{give.full} is true, return a \code{\link{list}}, say \code{r}, of
  \code{\link{unirootR}(.)} results, with \code{length(r) == length(p)}.

  Otherwise, return a \dQuote{numeric vector} like \code{p}, e.g., of
  \code{class "mpfr"} when \code{p} is.
}
\author{Martin Maechler}

\seealso{
 Standard \R's \code{\link[stats]{qnorm}}.
}
\examples{
doX <- Rmpfr:::doExtras() # slow parts only if(doX)
cat("doExtras: ", doX, "\n")
p  <- (0:32)/32
lp <- -c(1000, 500, 200, 100, 50, 20:1, 2^-(1:8))
if(doX) {
  tol1 <- 2.3e-16
  tolM <- 1e-20
  tolRIlog <- 4e-14
} else { # use one more than a third of the points:
   ip <- c(TRUE,FALSE, rep_len(c(TRUE,FALSE,FALSE), length(p)-2L))
   p <-  p[ip]
  lp <- lp[ip]
  tol1 <- 1e-9
  tolM <- 1e-12
  tolRIlog <- 25*tolM
}

f.all.eq <- function(a,b)
  sub("^Mean relative difference:", '', format(all.equal(a, b, tol=0)))
for(logp  in c(FALSE,TRUE)) {
  pp <- if(logp) lp else p
  mp <- mpfr(pp, precBits = if(doX) 80 else 64) # precBits = 128 gave "the same" as 80
  for(l.tail in c(FALSE,TRUE)) {
      qn <- qnorm (pp, lower.tail = l.tail, log.p = logp)
     qnI <- qnormI(pp, lower.tail = l.tail, log.p = logp, tol = tol1)
     qnM <- qnormI(mp, lower.tail = l.tail, log.p = logp, tol = tolM)
     cat(sprintf("Accuracy of qnorm(*, lower.t=\%-5s, log.p=\%-5s): \%s || qnI: \%s\n",
                 l.tail, logp, f.all.eq(qnM, qn ),
                               f.all.eq(qnM, qnI)))
     stopifnot(exprs = {
        all.equal(qn,  qnI, tol = if(logp) tolRIlog else 4*tol1)
        all.equal(qnM, qnI, tol = tol1)
     })
  }
}

## useMpfr, using mpfr()  :
if(doX) {
  p2 <- 2^-c(1:27, 5*(6:20), 20*(6:15))
  e2 <- 88
} else {
  p2 <- 2^-c(1:2, 7, 77, 177, 307)
  e2 <- 60
}
system.time( pn2 <- pnorm(qnormI(mpfr(p2, e2))) ) # 4.1 or 0.68
           all.equal(p2, pn2, tol = 0) # 5.48e-29 // 5.2e-18
2^-e2
stopifnot(all.equal(p2, pn2, tol = 6 * 2^-e2)) # '4 *' needed


## Boundary -- from limits in mpfr maximal exponent range!
## 1) Use maximal ranges:
(old_eranges <- .mpfr_erange()) # typically -/+ 2^30
(myERng <- (1-2^-52) * .mpfr_erange(c("min.emin","max.emax")))
(doIncr <- !isTRUE(all.equal(unname(myERng), unname(old_eranges))))
if(doIncr) .mpfr_erange_set(value = myERng)

log2(abs(.mpfr_erange()))# 62 62 if the above increase worked
(lrgOK <- all(log2(abs(.mpfr_erange())) >= 62))
                                        # currently(2023-01) FALSE on Windows
## The largest quantile for which our mpfr-ized qnorm() does *NOT* underflow :
cM <- if(doX) { "2528468770.343293436810768159197281514373932815851856314908753969469064"
      } else    "2528468770.34329343681"
##               1 3 5 7 9  1 3 5 7 9 1 3 5 7 9 1 3 5 7 9 1 3 5 7 9 1 3 5 7 9 1 3 5 7 9 1 3
##                       10         20        30        40        50        60        70
(qM <- mpfr(cM))
(pM <- pnorm(-qM)) # precision  233   bits  if(doX)
## 7.64890682545699845135633468495894619457903458325606933043966616334460003e-1388255822130839040
log(pM) # 233 bits: -3196577161300663205.8575919621115614148120323933633827052786873078552904

if(lrgOK) withAutoprint({
  %% FIXME : here, the internal qnInt() gives (-Inf, -Inf) :
  try( qnormI(pM) ) ## Error: lower < upper not fulfilled (evt. TODO)
  ## but this works
  print(qnI <- qnormI(log(pM), log.p=TRUE)) #  -2528468770.343293436
  all.equal(-qM, qnI, tol = 0) # << show how close; seen  1.084202e-19
  stopifnot( all.equal(-qM, qnI, tol = 1e-18) )
})

if(FALSE) # this (*SLOW*) gives 21 x the *same* (wrong) result --- FIXME!
  qnormI(log(pM) * (2:22), log.p=TRUE)
%% 21 'mpfr' numbers of precision 233 bits --- actually see *2* slightly differing:
%%  [1] -2528468770.34329343681076815919728151437393281585185631490875396946906388
%%  [2] -2528468770.34329343681076815919728151437393281585185631490875396946906388
%%  ...............
%% [20] -2528468770.34329343681076815919728151437393281585185631490875396946906419 <<
%% [21] -2528468770.34329343681076815919728151437393281585185631490875396946906388
if(doX) ## Show how bad it is (currently ca. 220 iterations, and then *wrong*)
 str(qnormI(round(log(pM)), log.p=TRUE, trace=1, give.full = TRUE))
if(requireNamespace("DPQ"))
  new("mpfr", as(DPQ::qnormR(pM, trace=1), "mpfr")) # as(*, "mpfr") also works for +/- Inf
  # qnormR1(p=         0, m=0, s=1, l.t.= 1, log= 0): q = -0.5
  #    somewhat close to 0 or 1: r := sqrt(-lp) =  1.7879e+09
  #    r > 5, using rational form R_3(t), for t=1.787897e+09  -- that is *not* accurate
  # [1] -94658744.369295865460462720............

## reset to previous status if needed
if(doIncr) .mpfr_erange_set( , old_eranges)
}
\keyword{distribution}
\keyword{math}
