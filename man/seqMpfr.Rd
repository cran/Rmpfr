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
  \item{length.out}{desired length of the sequence.  A
    non-negative number, which for \code{seq} and \code{seq.int} will be
    rounded up if fractional.}
  \item{along.with}{take the length from the length of this argument.}
  \item{\dots}{arguments passed to or from methods.}
}
\details{
  see \code{\link[base]{seq}} (default method in package \pkg{base}),
  whose semantic we want to replicate.
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
stopifnot(seqMpfr(7,3) == 7:3)
}
\keyword{manip}
