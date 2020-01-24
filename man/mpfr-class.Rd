\name{mpfr-class}
\docType{class}
\title{Class "mpfr" of Multiple Precision Floating Point Numbers}
%
\alias{mpfr-class}
\alias{mpfr1-class}
\alias{summaryMpfr-class}
%
\alias{[[,mpfr-method}
\alias{[,mpfr,ANY,missing,missing-method}
\alias{[<-,mpfr,missing,missing,ANY-method}
\alias{[<-,mpfr,ANY,missing,ANY-method}
\alias{[<-,mpfr,ANY,missing,mpfr-method}
% FIXME: many of these are just \alias{}ed, but *NOT* documented:
\alias{as.numeric,mpfr-method}
\alias{as.vector,mpfrArray-method}
\alias{as.integer,mpfr-method}
\alias{beta,ANY,mpfr-method}
\alias{beta,ANY,mpfrArray-method}
\alias{beta,mpfr,ANY-method}
\alias{beta,mpfr,mpfr-method}
\alias{beta,mpfr,numeric-method}
\alias{beta,numeric,mpfr-method}
\alias{beta,mpfrArray,ANY-method}
\alias{beta,mpfrArray,mpfrArray-method}
\alias{lbeta,ANY,mpfr-method}
\alias{lbeta,ANY,mpfrArray-method}
\alias{lbeta,mpfr,ANY-method}
\alias{lbeta,mpfr,mpfr-method}
\alias{lbeta,mpfr,numeric-method}
\alias{lbeta,numeric,mpfr-method}
\alias{lbeta,mpfrArray,ANY-method}
\alias{lbeta,mpfrArray,mpfrArray-method}
\alias{atan2,ANY,mpfr-method}
\alias{atan2,ANY,mpfrArray-method}
\alias{atan2,mpfr,ANY-method}
\alias{atan2,mpfr,mpfr-method}
\alias{atan2,mpfr,numeric-method}
\alias{atan2,numeric,mpfr-method}
\alias{atan2,mpfrArray,ANY-method}
\alias{atan2,mpfrArray,mpfrArray-method}
\alias{hypot}
%
\alias{coerce,mpfr,character-method}
\alias{coerce,mpfr,numeric-method}
\alias{coerce,mpfr,bigz-method}
\alias{coerce,mpfr,integer-method}
\alias{coerce,mpfr1,numeric-method}
\alias{coerce,mpfr1,mpfr-method}
\alias{coerce,integer,mpfr-method}
\alias{coerce,logical,mpfr-method}
\alias{coerce,raw,mpfr-method}
\alias{coerce,numeric,mpfr-method}
\alias{coerce,numeric,mpfr1-method}
\alias{coerce,array,mpfr-method}
\alias{coerce,character,mpfr-method}
\alias{coerce,mpfr,mpfr1-method}
%
\alias{Ops,mpfr,ANY-method}
\alias{Ops,ANY,mpfr-method}
\alias{Ops,mpfr,bigq-method}
\alias{Ops,bigq,mpfr-method}
\alias{Ops,mpfr,bigz-method}
\alias{Ops,bigz,mpfr-method}
\alias{Ops,array,mpfr-method}
\alias{Ops,mpfr,array-method}
\alias{Ops,mpfr,vector-method}
\alias{Ops,vector,mpfr-method}
\alias{Arith,mpfr,array-method}
\alias{Arith,mpfr,missing-method}
\alias{Arith,mpfr,mpfr-method}
\alias{Arith,mpfr,integer-method}
\alias{Arith,mpfr,numeric-method}
\alias{Arith,integer,mpfr-method}
\alias{Arith,numeric,mpfr-method}
\alias{Arith,array,mpfr-method}
%
\alias{Compare,mpfr,mpfr-method}
\alias{Compare,mpfr,integer-method}
\alias{Compare,mpfr,numeric-method}
\alias{Compare,integer,mpfr-method}
\alias{Compare,numeric,mpfr-method}
\alias{Compare,mpfr,array-method}
\alias{Compare,array,mpfr-method}
%
\alias{Logic,mpfr,mpfr-method}
\alias{Logic,mpfr,numeric-method}
\alias{Logic,numeric,mpfr-method}
\alias{Summary,mpfr-method}
\alias{Math,mpfr-method}
\alias{Math2,mpfr-method}
\alias{abs,mpfr-method}
\alias{log,mpfr-method}
\alias{factorial,mpfr-method}
\alias{sign,mpfr-method}
% "Complex" (cheap methods for "real"s):
\alias{Re,mpfr-method}
\alias{Im,mpfr-method}
\alias{Mod,mpfr-method}
\alias{Arg,mpfr-method}
\alias{Conj,mpfr-method}

\alias{format,mpfr-method}
\alias{is.finite,mpfr-method}
\alias{is.infinite,mpfr-method}
\alias{is.na,mpfr-method}
\alias{is.nan,mpfr-method}
\alias{is.finite,mpfrArray-method}
\alias{is.infinite,mpfrArray-method}
\alias{is.na,mpfrArray-method}
\alias{is.nan,mpfrArray-method}
\alias{unique,mpfr-method}% <---  wanted        by  R CMD check
\alias{unique,mpfr,ANY-method}% not sufficient for  R CMD ..
\alias{unique.mpfr}
%
\alias{all.equal,mpfr,mpfr-method}
\alias{all.equal,mpfr,ANY-method}
\alias{all.equal,ANY,mpfr-method}
\alias{mean,mpfr-method}
\alias{median,mpfr-method}
\alias{quantile,mpfr-method}
\alias{summary,mpfr-method}

\alias{dim<-,mpfr-method}
\alias{t,mpfr-method}

\alias{\%*\%,array_or_vector,mpfr-method}
\alias{\%*\%,mpfr,array_or_vector-method}
\alias{crossprod,array_or_vector,mpfr-method}
\alias{crossprod,mpfr,array_or_vector-method}
\alias{tcrossprod,array_or_vector,mpfr-method}
\alias{tcrossprod,mpfr,array_or_vector-method}

\alias{which.min,mpfr-method}
\alias{which.max,mpfr-method}

\alias{show,mpfr-method}
\alias{show,mpfr1-method}
\alias{show,summaryMpfr-method}
\alias{print.mpfr1}
\alias{print.summaryMpfr}
%
\description{
  \code{"mpfr"} is the class of \bold{M}ultiple \bold{P}recision
  \bold{F}loatingpoint numbers with \bold{R}eliable arithmetic.

  sFor the high-level user, \code{"mpfr"} objects should behave
  as standard \R's \code{\link{numeric}} \emph{vectors}. They would just
  print differently and use the prespecified (typically high) precision
  instead of the double precision of \sQuote{traditional} \R numbers
  (with \code{\link{class}(.) == "numeric"} and
  \code{\link{typeof}(.) == "double"}).

  \code{hypot(x,y)} computes the hypothenuse length \eqn{z} in a rectangular
  triangle with \dQuote{leg} side lengths \eqn{x} and \eqn{y}, i.e.,
  \deqn{z = hypot(x,y) = \sqrt{x^2 + y^2},}{z = hypot(x,y) = sqrt(x^2 + y^2),}
  in a numerically stable way.
}
\usage{
hypot(x,y, rnd.mode = c("N","D","U","Z","A"))
}
\arguments{
  \item{x,y}{an object of class \code{mpfr}.}
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see \code{\link{mpfr}}.}
}
\section{Objects from the Class}{
  Objects are typically created by \code{\link{mpfr}(<number>, precBits)}.

  \code{summary(<mpfr>)} returns an object of class \code{"summaryMpfr"}
  which contains \code{"mpfr"} but has its own \code{\link{print}} method.
}
\section{Slots}{
  Internally, \code{"mpfr"} objects just contain standard \R
  \code{\link{list}}s where each list element is of class
  \code{"mpfr1"}, representing \emph{one} MPFR number, in a structure
  with four slots, very much parallelizing the C \code{struc} in the
  \code{mpfr} C library to which the \pkg{Rmpfr} package interfaces.

  An object of class \code{"mpfr1"} has slots
  \describe{
    \item{\code{prec}:}{\code{"integer"} specifying the maxmimal
      precision in \bold{bits}.}
    \item{\code{exp}:}{\code{"integer"} specifying the base-\bold{2}
      exponent of the number.}
    \item{\code{sign}:}{\code{"integer"}, typically \code{-1} or
      \code{1}, specifying the sign (i.e. \code{\link{sign}(.)}) of the
      number.}
    \item{\code{d}:}{an \code{"integer"} vector (of 32-bit
      \dQuote{limbs}) which corresponds to the full mantissa of the
      number.}
  }% desc
}
% \section{Extends}{
% Class \code{"\linkS4class{list}"}, from data part.
% Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
% }
\section{Methods}{
  \describe{
    \item{abs}{\code{signature(x = "mpfr")}: ... }

    \item{atan2}{\code{signature(y = "mpfr", x = "ANY")}, and}
    \item{atan2}{\code{signature(x = "ANY", y = "mpfr")}: compute the
     arc-tangent of two arguments: \code{atan2(y, x)} returns the angle
     between the x-axis and the vector from the origin to \eqn{(x, y)},
     i.e., for positive arguments \code{atan2(y, x) == \link{atan}(y/x)}.}

    \item{lbeta}{\code{signature(a = "ANY", b = "mpfrArray")}, is
      \eqn{\log(|B(a,b)|)}{log(abs(B(a,b)))} where \eqn{B(a,b)} is the
      Beta function, \code{beta(a,b)}.}
    \item{beta}{\code{signature(a = "mpfr", b = "ANY")},}
    \item{beta}{\code{signature(a = "mpfr", b = "mpfr")}, \dots, etc:
      Compute the beta function \eqn{B(a,b)}, using high precision,
      building on internal \code{\link{gamma}} or \code{\link{lgamma}}.
      See the help for \R's base function \code{\link[base]{beta}} for
      more.  Currently, there, \eqn{a,b \ge 0}{a,b >= 0} is required.
      Here, we provide (non-\code{\link{NaN}}) for all numeric \code{a, b}.

      When either \eqn{a}, \eqn{b}, or \eqn{a+b} is a negative
      \emph{integer}, \eqn{\Gamma(.)} has a pole there and is undefined
      (\code{NaN}).  However the Beta function can be defined there as
      \dQuote{limit}, in some cases.  Following other software such as
      SAGE, Maple or Mathematica, we provide finite values in these
      cases.  However, note that these are not proper limits
      (two-dimensional in \eqn{(a,b)}), but useful for some
      applications.  E.g., \eqn{B(a,b)} is defined as zero when
      \eqn{a+b} is a negative integer, but neither \eqn{a} nor \eqn{b} is.
      Further, if \eqn{a > b > 0} are integers, \eqn{B(-a,b)= B(b,-a)}
      can be seen as \eqn{(-1)^b * B(a-b+1,b)}.
    }
    \item{dim<-}{\code{signature(x = "mpfr")}: Setting a dimension
      \code{\link{dim}} on an \code{"mpfr"} object makes it into an object
      of class \code{"\linkS4class{mpfrArray}"} or (more specifically)
      \code{"mpfrMatrix"} for a length-2 dimension, see their help page;
      note that \code{t(x)} (below) is a special case of this.}

    \item{Ops}{\code{signature(e1 = "mpfr", e2 = "ANY")}: ... }
    \item{Ops}{\code{signature(e1 = "ANY", e2 = "mpfr")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "missing")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "mpfr")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "integer")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "numeric")}: ... }
    \item{Arith}{\code{signature(e1 = "integer", e2 = "mpfr")}: ... }
    \item{Arith}{\code{signature(e1 = "numeric", e2 = "mpfr")}: ... }
    \item{Compare}{\code{signature(e1 = "mpfr", e2 = "mpfr")}: ... }
    \item{Compare}{\code{signature(e1 = "mpfr", e2 = "integer")}: ... }
    \item{Compare}{\code{signature(e1 = "mpfr", e2 = "numeric")}: ... }
    \item{Compare}{\code{signature(e1 = "integer", e2 = "mpfr")}: ... }
    \item{Compare}{\code{signature(e1 = "numeric", e2 = "mpfr")}: ... }

    \item{Logic}{\code{signature(e1 = "mpfr", e2 = "mpfr")}: ... }

    \item{Summary}{\code{signature(x = "mpfr")}: The S4
      \code{\link[methods:S4groupGeneric]{Summary}} group functions,
      \code{\link{max}}, \code{\link{min}}, \code{\link{range}},
      \code{\link{prod}}, \code{\link{sum}},
      \code{\link{any}}, and \code{\link{all}}
      are all defined for MPFR numbers.  \code{\link{mean}(x, trim)} for
      non-0 \code{trim} works analogously to \code{\link{mean.default}}.}

    \item{median}{\code{signature(x = "mpfr")}: works via}% quantile(*, probs = 0.5) :
    \item{quantile}{\code{signature(x = "mpfr")}: a simple wrapper of
      the \code{quantile.default} method from \pkg{stats}.}
    \item{summary}{\code{signature(object = "mpfr")}: modeled after
      \code{\link{summary.default}}, ensuring to provide the full "mpfr"
      range of numbers.}

    \item{Math}{\code{signature(x = "mpfr")}: All the S4
      \code{\link[methods:S4groupGeneric]{Math}} group functions are
      defined, using multiple precision (MPFR) arithmetic, from
      \code{\link{getGroupMembers}("Math")}, these are (in alphabetical order):

      \code{\link{abs}}, \code{\link{sign}}, \code{\link{sqrt}},
      \code{\link{ceiling}}, \code{\link{floor}}, \code{\link{trunc}},
      \code{\link{cummax}}, \code{\link{cummin}}, \code{\link{cumprod}},
      \code{\link{cumsum}}, \code{\link{exp}}, \code{\link{expm1}},
      \code{\link{log}}, \code{\link{log10}}, \code{\link{log2}},
      \code{\link{log1p}}, \code{\link{cos}}, \code{\link{cosh}},
      \code{\link{sin}}, \code{\link{sinh}}, \code{\link{tan}},
      \code{\link{tanh}}, \code{\link{acos}}, \code{\link{acosh}},
      \code{\link{asin}}, \code{\link{asinh}}, \code{\link{atan}},
      \code{\link{atanh}}, \code{\link{gamma}}, \code{\link{lgamma}},
      \code{\link{digamma}}, and \code{\link{trigamma}}.

      Currently, \code{trigamma} is not provided by
      the MPFR library and hence not yet implemented.
      \cr
      Further, the \code{cum*()} methods are \emph{not yet} implemented.}

    \item{factorial}{\code{signature(x = "mpfr")}: this will
      \code{\link{round}} the result when \code{x} is integer valued.
      Note however that \code{\link{factorialMpfr}(n)} for integer
      \code{n} is slightly more efficient, using the MPFR function
      \samp{mpfr_fac_ui}.}

    \item{Math2}{\code{signature(x = "mpfr")}: \code{\link{round}(x,
	digits)} and \code{\link{signif}(x, digits)} methods.  Note that
	these do not change the formal precision (\code{'prec'} slot),
	and you may often want to apply \code{\link{roundMpfr}()} in
	addition or preference.}

    \item{as.numeric}{\code{signature(x = "mpfr")}: ... }

    \item{as.vector}{\code{signature(x = "mpfrArray")}: as for standard
      \code{\link{array}}s, this \dQuote{drops} the \code{dim} (and
      \code{dimnames}), i.e., transforms \code{x} into an \sQuote{MPFR}
      number vector, i.e., class \code{\linkS4class{mpfr}}.}

    \item{[[}{\code{signature(x = "mpfr", i = "ANY")}, and}
    \item{[}{\code{signature(x = "mpfr", i = "ANY", j = "missing", drop = "missing")}:
      subsetting aka \dQuote{indexing} happens as for numeric vectors.}

    \item{format}{\code{signature(x = "mpfr")}, further arguments
      \code{digits = NULL, scientific = NA}, etc:
      returns \code{\link{character}} vector of same length as \code{x};
      when \code{digits} is \code{NULL}, with \emph{enough} digits to
      recreate \code{x} accurately.  For details, see
      \code{\link{formatMpfr}}.}
    \item{is.finite}{\code{signature(x = "mpfr")}: ... }
    \item{is.infinite}{\code{signature(x = "mpfr")}: ... }
    \item{is.na}{\code{signature(x = "mpfr")}: ... }
    \item{is.nan}{\code{signature(x = "mpfr")}: ... }
    \item{log}{\code{signature(x = "mpfr")}: ... }
    \item{show}{\code{signature(object = "mpfr")}: ... }
    \item{sign}{\code{signature(x = "mpfr")}: ... }

    \item{Re, Im}{\code{signature(z = "mpfr")}: simply return \code{z}
      or \code{0} (as \code{"mpfr"} numbers of correct precision), as mpfr
      numbers are \sQuote{real} numbers.}
    \item{Arg, Mod, Conj}{\code{signature(z = "mpfr")}: these are
      trivial for our \sQuote{real} mpfr numbers, but defined to work
      correctly when used in \R code that also allows complex number input.}

    \item{all.equal}{\code{signature(target = "mpfr", current = "mpfr")},}
    \item{all.equal}{\code{signature(target = "mpfr", current = "ANY")}, and}
    \item{all.equal}{\code{signature(target = "ANY", current = "mpfr")}:
      methods for numerical (approximate) equality,
      \code{\link[base]{all.equal}} of multiple precision numbers.  Note
      that the default \code{tolerance} (argument) is taken to correspond
      to the (smaller of the two) precisions when both main arguments are
      of class \code{"mpfr"}, and hence can be considerably less than
      double precision machine epsilon \code{\link{.Machine}$double.eps}.}
    \item{coerce}{\code{signature(from = "numeric", to = "mpfr")}:
      \code{\link{as}(., "mpfr")} coercion methods are available for
      \code{\link{character}} strings, \code{\link{numeric}}, \code{\link{integer}},
      \code{\link{logical}}, and even \code{\link{raw}}.  Note however,
      that \code{\link{mpfr}(., precBits, base)} is more flexible.}
    \item{coerce}{\code{signature(from = "mpfr", to = "bigz")}: coerces
      to biginteger, see \code{\link[gmp]{bigz}} in package \pkg{gmp}.}
    \item{coerce}{\code{signature(from = "mpfr", to = "numeric")}: ... }
    \item{coerce}{\code{signature(from = "mpfr", to = "character")}: ... }

    \item{unique}{\code{signature(x = "mpfr")}, and corresponding S3 method
      (such that \code{unique(<mpfr>)} works inside \pkg{base} functions),
      see \code{\link{unique}}.

      Note that \code{\link{duplicated}()} works for \code{"mpfr"} objects
      without the need for a specific method.}

    \item{t}{\code{signature(x = "mpfr")}: makes \code{x} into an \eqn{n
	\times 1}{n x 1} \code{\linkS4class{mpfrMatrix}}.}

    \item{which.min}{\code{signature(x = "mpfr")}: gives the index of
	the first minimum, see \code{\link{which.min}}.}
    \item{which.max}{\code{signature(x = "mpfr")}: gives the index of
	the first maximum, see \code{\link{which.max}}.}
  }
}
\author{Martin Maechler}

\note{Many more methods (\dQuote{functions}) automagically work for
  \code{"mpfr"} number vectors (and matrices, see the
  \code{\linkS4class{mpfrMatrix}} class doc),
  notably
  \code{\link{sort}}, \code{\link{order}}, \code{\link{quantile}},
  \code{\link{rank}}.
}

\seealso{
  The \code{"\linkS4class{mpfrMatrix}"} class, which extends the
  \code{"mpfr"} one.

  \code{\link{roundMpfr}} to \emph{change} precision of an \code{"mpfr"}
  object which is typically desirable \emph{instead} of or in addition
  to \code{signif()} or \code{round()};
  \code{\link{is.whole}()} etc.

  Special mathematical functions such as some Bessel ones, e.g., \code{\link{jn}};
  further, \code{\link{zeta}(.)} \eqn{(= \zeta(.))}, \code{Ei()} etc.
  \code{\link{Bernoulli}} numbers and the Pochhammer function
  \code{\link{pochMpfr}}.
}
\examples{
## 30 digit precision
str(x <- mpfr(c(2:3, pi), prec = 30 * log2(10)))
x^2
x[1] / x[2] # 0.66666... ~ 30 digits

## indexing - as with numeric vectors
stopifnot(identical(x[2], x[[2]]),
	  ## indexing "outside" gives NA (well: "mpfr-NaN" for now):
	  is.na(x[5]),
	  ## whereas "[[" cannot index outside:
	  is(try(x[[5]]), "try-error"),
	  ## and only select *one* element:
	  is(try(x[[2:3]]), "try-error"))

## factorial() & lfactorial would work automagically via [l]gamma(),
## but factorial() additionally has an "mpfr" method which rounds
f200 <- factorial(mpfr(200, prec = 1500)) # need high prec.!
f200
as.numeric(log2(f200))# 1245.38 -- need precBits >~ 1246 for full precision

##--> see  factorialMpfr() for more such computations.

##--- "Underflow" **much** later -- exponents have 30(+1) bits themselves:

mpfr.min.exp2 <- - (2^30 + 1)
two <- mpfr(2, 55)
stopifnot(two ^ mpfr.min.exp2 == 0)
## whereas
two ^ (mpfr.min.exp2 * (1 - 1e-15))
## 2.38256490488795107e-323228497   ["typically"]

##--- "Assert" that {sort}, {order}, {quantile}, {rank}, all work :

p <- mpfr(rpois(32, lambda=500), precBits=128)^10
np <- as.numeric(log(p))
(sp <- summary(p))# using the print.summaryMpfr() method
stopifnot(all(diff(sort(p)) >= 0),
   identical(order(p), order(np)),
   identical(rank (p), rank (np)),
   all.equal(sapply(1:9, function(Typ) quantile(np, type=Typ, names=FALSE)),
      sapply(lapply(1:9, function(Typ) quantile( p, type=Typ, names=FALSE)),
	     function(x) as.numeric(log(x))),
      tol = 1e-3),# quantiles: interpolated in orig. <--> log scale
 TRUE)

m0 <- mpfr(numeric(), 99)
xy <- expand.grid(x = -2:2, y = -2:2) ; x <- xy[,"x"] ; y <- xy[,"y"]
a2. <- atan2(y,x)

stopifnot(identical(which.min(m0), integer(0)),
	  identical(which.max(m0), integer(0)),
          all.equal(a2., atan2(as(y,"mpfr"), x)),
	  max(m0) == mpfr(-Inf, 53), # (53 is not a feature, but ok)
	  min(m0) == mpfr(+Inf, 53),
	  sum(m0) == 0, prod(m0) == 1)

## unique(), now even base::factor()  "works" on <mpfr> :
set.seed(17)
p <- rlnorm(20) * mpfr(10, 100)^-999
pp <- sample(p, 50, replace=TRUE)
str(unique(pp)) # length 18 .. (from originally 20)
## Class 'mpfr' [package "Rmpfr"] of length 18 and precision 100
## 5.56520587824e-999 4.41636588227e-1000 ..
facp <- factor(pp)
str(facp) # the factor *levels* are a bit verbose :
# Factor w/ 18 levels "new(\"mpfr1\", ...........)" ...
# At least *some* factor methods work :
stopifnot(exprs = {
  is.factor(facp)
  identical(unname(table(facp)),
            unname(table(asNumeric(pp * mpfr(10,100)^1000))))
})

## ((unfortunately, the expressions are wrong; should integer "L"))
#
## More useful: levels with which to *invert* factor() :
## -- this is not quite ok:
## simplified from 'utils' :
deparse1 <- function(x, ...) paste(deparse(x, 500L, ...), collapse = " ")
if(FALSE) {
 str(pp.levs <- vapply(unclass(sort(unique(pp))), deparse1, ""))
 facp2 <- factor(pp, levels = pp.levs)
}
}
\keyword{classes}
