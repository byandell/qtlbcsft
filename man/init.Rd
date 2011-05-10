\name{bcsft.init}
\alias{bcsft.init}
\alias{bcsft.step}
\alias{bcsft.nrec}
\alias{bcsft.rf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Examine low-level C routines for QTL calculations.
}
\description{
BC1F0 is a backcross, and BC0F2 is an intercross. The BCsFt calculations
are compared for the low-level C routines used in R/qtl.
}
\usage{
bcsft.init(s = 1, t = 0)
bcsft.step(rf = 0.5, s = 0, t = 2)
bcsft.nrec(rf = 0.5, s = 0, t = 2)
bcsft.rf(rf = 0.5, s = 0, t = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rf}{
Recombination frequence (between 0 and 0.5).
}
  \item{s}{
Number of backcross generations.
}
  \item{t}{
Number of intercross generations. If \code{s} = 0, then \code{t} must be
at least 2, as F1 is zero-th generation for all experimental crosses.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Brian S. Yandell <byandell@wisc.edu> and Laura Shannon.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
bcsft.init(s = 1, t = 0)
bcsft.step(rf = 0.5, s = 0, t = 2)
bcsft.nrec(rf = 0.5, s = 0, t = 2)
bcsft.rf(rf = 0.5, s = 0, t = 2)
}
\keyword{utilities}
