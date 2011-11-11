\name{apple}
\alias{apple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximate Path for Penalized Likelihood Estimator
}
\description{
Fit a generalized linear model via penalized maximum likelihood. The regularization path is computed for the LASSO or MCP penalty at a grid of values for the regularization parameter lambda. Can deal all shapes of data, including very large sparse data matrices. Fits binomial-logistic and poisson regression models. 
}
\usage{
apple(X, y, family="binomial", penalty = "LASSO", gamma, cha.poi = 1, 
eps = 1e-15, lam.list, lambda.min.ratio, max.iter = 100, max.num, 
n.lambda = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
input matrix, of dimension nobs x nvars; each row is an observation vector.
}
  \item{y}{
response variable, of dimension nobs x 1. non-negative counts for 

\code{family="poisson"}, binary for \code{family="binomial"}.
}
  \item{family}{
response type.
}
  \item{penalty}{
\code{LASSO} and \code{MCP} are provided.
}
  \item{gamma}{
the MCP concavity parameter.
}
  \item{cha.poi}{
the value used to change from Newton Raphson correction to Coordinate Descent correction, which is the \eqn{\alpha} in the following inequality, \eqn{k> \alpha\sqrt{n}}, where \eqn{k} is the size of current active set. when this inequality holds, the correction method changes from Newton Raphson to Coordinate Descent.
}
  \item{eps}{
the precision used to test the convergence.
}
  \item{lam.list}{
a user supplied \eqn{\lambda} sequence. typical usage is to have the program compute its own \code{lambda} sequence based on \code{lambda.min.ratio} and \code{n.lambda}. supplying a value of \eqn{\lambda} overrides this.
}
  \item{lambda.min.ratio}{
optional input. smallest value for \code{lambda}, as a fraction of \code{max.lam}, the (data derived) entry value. the default depends on the sample size \code{n} relative to the number of variables \code{p}. if \code{n} > \code{p}, the default is 0.0001. otherwise, the default is 0.01.
}
  \item{max.iter}{
maximum number of iteration in the computation.
}
  \item{max.num}{
optional input. maximum number of nonzero coefficients.
}
  \item{n.lambda}{
the number of \code{lambda} values.
}
}
\value{
\item{a0}{intercept vector of length(\code{lambda}).}
\item{beta}{nvar x length(\code{lambda}) matrix of coefficients.}
\item{lambda}{the list of lambda derived the solution path.}
\item{ebic}{the list of EBIC values.}
\item{ebic.loc}{the location of the EBIC selected solution in the path.}
\item{family}{the family of the supplied dataset.}
}
\references{
Yi Yu and Yang Feng, APPLE: Approximate Path for Penalized Likelihood Estimator, manuscript.
}
\author{
Yi Yu and Yang Feng
}
\seealso{
\code{\link{plot.apple}}, \code{\link{cv.apple}} and \code{\link{predict.apple}}
}
\examples{
%%%generate data
p=10
n=200
d=5
coefs=c(3,1.5,0,0,2)
intercept=0
beta=rep(0,p)
beta[1:d]=coefs
X=matrix(rnorm(p*n), nrow=n)
mu=1/(1+exp(-X \%*\% beta-intercept))
y=rbinom(n,1,mu)
	
%%%fit the model
fit.apple=apple(X, y, family='binomial')

%%%plot the path
plot(fit.apple)
}

