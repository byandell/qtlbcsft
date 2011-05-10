\name{cmp.rf}
\alias{cmp.rf}
\title{
Compare values for two objects created by est.rf.
}
\description{
Detailed comparison of two objects created by est.rf.
}
\usage{
cmp.rf(rf1, rf2, tol = 1e-06)
}
\arguments{
  \item{rf1}{
Recombination frequency (between 0 and 0.5).
}
  \item{rf2}{
Recombination frequency (between 0 and 0.5).
}
  \item{tol}{
Tolerance level.
}
}
\details{
Comparison of summaries and pattern of missing and zero values for
matrix that contains number of meioses in diagonal, \code{LOD} above diagonal,
and \code{rf} below diagonal.
}
\value{
List with values
  \item{lod.diff}{Matrix showing row, column and LOD values that differ}
  \item{rf.diff}{Matrix showing row, column and rf values that differ}
}
\author{
Brian S. Yandell <byandell@wisc.edu> and Laura Shannon.
}
\seealso{
\code{\link[qtl]{est.rf}}, \code{\link[qtl]{plot.rf}}
}
\examples{
data(fake.f2)
fake.ft <- fake.f2
class(fake.ft)[1] <- "bcsft"
attr(fake.ft, "scheme") <- c(0,2)
rf.f2 <- est.rf(fake.f2)
rf.ft <- est.rf(fake.ft)
cmp.rf(rf.f2,rf.ft)
}
\keyword{utilities}
