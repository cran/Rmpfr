\name{seqMpfr}
%% This is to contain the "true" seq() methods --- once we can ! ---
\Rdversion{1.1}
\alias{seqMpfr}
\title{"mpfr" Sequence Generation}
\description{
  Generate \sQuote{regular}, i.e., arithmetic sequences.  This is in
  lieu of methods for \code{\link{seq}} (dispatching on all three of
  \code{from}, \code{to}, and \code{by}.
}
\usage{
seqMpfr(from = 1, to = 1, by = ((to - from)/(length.out - 1)),
        length.out = NULL, along.with = NULL, \dots)
}
\arguments{
  \item{from, to}{the starting and (maximal) end value (numeric or
    \code{"\linkS4class{mpfr}"}) of the sequence.}
  \item{by}{number (numeric or \code{"\linkS4class{mpfr}"}): increment
    of the sequence.}
  \item{length.out}{desired length of the sequence.  A non-negative
    number, which will be rounded up if fractional.}
  \item{along.with}{take the length from the length of this argument.}
  \item{\dots}{arguments passed to or from methods.}
}
\details{
  see \code{\link[base]{seq}} (default method in package \pkg{base}),
  whose semantic we want to replicate (almost).% not  seq(9:6) => seq(4)
}
\value{
  a \sQuote{vector} of class \code{"\linkS4class{mpfr}"}, when one of
  the first three arguments was.
}
\author{Martin Maechler}
\seealso{
  The documentation of the \pkg{base} function \code{\link[base]{seq}};
  \code{\link{mpfr}}
}
\examples{
seqMpfr(0, 1, by = mpfr(0.25, prec=88))

seqMpfr(7, 3) # -> default prec.
\dontshow{
eq.test <- function(...) {
  args <- list(...)
  r <- do.call(seqMpfr, args)
  is(r, "mpfr") && all(r == do.call(seq, lapply(args, as.numeric)))
}
m <- mpfr(8, 60); m2 <- as(2, "mpfr")
stopifnot(seqMpfr(7,3) == 7:3, seqMpfr(to=1) == 1,
          eq.test(pi),
          eq.test(1,12,3),
          eq.test(1, by= 2, length= 5),
          eq.test(1, by=m2, length= 5),
          eq.test(1, length.out = 8),
          eq.test(1, length.out = m),
          eq.test(1, 8, length.out=8),
          eq.test(1, 8, length.out=m),
          eq.test(1, m, length.out=m),
          eq.test(1, m, length.out=8),
          eq.test(to=17, by= 2, length=4),
          eq.test(to=17, by=m2, length=4),
          TRUE)
}% {dontshow} -- but test
}
\keyword{manip}
