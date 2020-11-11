\name{igamma}
\alias{igamma}
\title{Incomplete Gamma Function}
\description{% >> ../R/special-fun.R <<<
  For MPFR version >= 3.2.0, the following MPFR library function is provided:
  \code{mpfr_gamma_inc(a,x)}, the \R interface of which is \code{igamma(a,x)}, where
  \code{igamma(a,x)} is the \dQuote{upper} incomplete gamma function
  %% fails in LaTeX (R "bug"): \deqn{γ(a,x) :=: Γ(a) - Γ(a,x),}
  \deqn{\Gamma(a,x) :=: \Gamma(a) - \gamma(a,x),}{Γ(a,x) :=: Γ(a) - γ(a,x),}
  where
  \deqn{\gamma(a,x) := \int_0^x      t^{a-1} e^{-t} dt,}{γ(a,x) := ∫₀ˣ tᵃ⁻¹ e⁻ᵗ dt,}
  and hence
  \deqn{\Gamma(a,x) := \int_x^\infty t^{a-1} e^{-t} dt,}{Γ(a,x) := ∫ (x..∞) tᵃ⁻¹ e⁻ᵗ dt,}
  and \deqn{\Gamma(a) := \gamma(a, \infty).}{Γ(a) := γ(a, ∞).}

  As \R's \code{\link{pgamma}(x,a)} is
  \deqn{\code{pgamma(x, a)} := \gamma(a,x) / \Gamma(a),}{\code{pgamma(x, a)} := γ(a,x) / Γ(a),}
  we get \preformatted{        igamma(a,x) ==  gamma(a) * pgamma(x, a, lower.tail=FALSE)}
}
\usage{
igamma(a, x, rnd.mode = c("N", "D", "U", "Z", "A"))
}
%% MM FIXME: "Copy paste" from hypot() -- which strangely we have in >> ./mpfr-class.Rd <<
\arguments{
  \item{a, x}{an object of class \code{mpfr} or \code{\link{numeric}}.}
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see \code{\link{mpfr}}.}
}
%% \details{

%% }
\value{
  a numeric vector of \dQuote{common length}, recyling along \code{a} and \code{x}.
}
\references{
  NIST Digital Library of Mathematical Functions, section 8.2.
  \url{https://dlmf.nist.gov/8.2.i}

  Wikipedia (2019). \emph{Incomplete gamma function};
  \url{https://en.wikipedia.org/wiki/Incomplete_gamma_function}

  %% .. arXiv  paper TODO % see --> ../TODO
}
\author{\R interface: Martin Maechler}
\seealso{
  \R's \code{\link{gamma}} (function) and \code{\link{pgamma}} (probability distribution).
}
\examples{
## show how close pgamma() is :
x <- c(seq(0,20, by=1/4), 21:50, seq(55, 100, by=5))
if(mpfrVersion() >= "3.2.0") { print(
all.equal(igamma(Const("pi", 80), x),
          pgamma(x, pi, lower.tail=FALSE) * gamma(pi),
          tol=0, formatFUN = function(., ...) format(., digits = 7)) #-> 3.13e-16
)
## and ensure *some* closeness:
stopifnot(exprs = {
   all.equal(igamma(Const("pi", 80), x),
             pgamma(x, pi, lower.tail=FALSE) * gamma(pi),
             tol = 1e-15)
})
} # only if MPFR version >= 3.2.0
}
\keyword{math}
