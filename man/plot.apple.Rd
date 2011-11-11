\name{plot.apple}
\alias{plot.apple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot apple path
}
\description{
Plot the solution path generated from apple.
}
\usage{
\method{plot}{apple}(x, col = "black", add = FALSE, 
main = "apple", type = "l", lty = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the object used to plot the path. object derived by \code{"apple"} and \code{"cv.apple"} are both available here. for \code{"cv.apple"} object, ebic and cv selected solutions are also lined in the graph. 
}
  \item{col}{
color used to plot the paths, for details please see the usage of \code{col} in \code{plot}.
}
  \item{add}{
whether or not to add this plot to an existing one.
}
  \item{main}{
title of the plot. default is \code{main="apple"}
}
  \item{type}{
what type of plot should be drawn, for details please see the usage of \code{plot}.
}
  \item{lty}{
the line type, for details please see the usage of \code{par}
}
  \item{...}{
see \code{matplot}.
}
}
\details{
if the object is a cv result, then both of ebic and cv selected result will be marked on the graph by solid and dotted vertical lines, respectively.
}
\references{
Yi Yu and Yang Feng, APPLE: Approximate Path for Penalized Likelihood Estimator, manuscript.
}
\author{
Yi Yu and Yang Feng
}

\seealso{
\code{\link{apple}}, \code{\link{cv.apple}} and \code{\link{predict.apple}}
}
\examples{
%% see \code{\link{apple}}, \code{\link{cv.apple}} 
}