\name{cv.apple}
\alias{cv.apple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cross validation for apple path
}
\description{
Does K-fold cross validation for apple.
}
\usage{
cv.apple(X, y, family="binomial", penalty="LASSO", gamma,  K = 10, 
alpha=0, seed = 1, cha.poi = 1, eps = 1e-15, lambda.min.ratio, 
max.iter = 100, max.num, n.lambda = 100)
}

\arguments{
  \item{X}{
input matrix, of dimension nobs x nvars; each row is an observation vector.}
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
  \item{K}{
number of folds used in cross validation. The default it 10.
}
  \item{alpha}{
weight used in the cross validation cost function, with \eqn{Q (\lambda) =\alpha Dev(\lambda)+ (1-\alpha) s(\lambda) \log n/n.}

}
  \item{seed}{
random seed used to sample the training sets and test sets
}
  \item{cha.poi}{
the value used to change from Newton Raphson correction to Coordinate Descent correction, which is the \eqn{\alpha} in the following inequality, \eqn{k> \alpha\sqrt{n}}, where \eqn{k} is the size of current active set. when this inequality holds, the correction method changes from Newton Raphson to Coordinate Descent.
}
  \item{eps}{
the precision used to test the convergence.
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
\item{cv}{list of cross validation loss}
\item{lambda}{list of lambda}
\item{a0}{the list of intercept}
\item{beta}{the list of coefficients}
\item{cv.loc}{location of cv selected solution in the path}
\item{ebic.loc}{the location of the EBIC selected solution in the path}
\item{cv.beta}{cross validation selected beta}
\item{ebic.beta}{ebic selected beta}
\item{cv.a0}{cv selected intercept}
\item{ebic.a0}{ebic selected intercept}
}
\references{
Yi Yu and Yang Feng, APPLE: Approximate Path for Penalized Likelihood Estimator, manuscript.
}
\author{
Yi Yu and Yang Feng
}
\seealso{
\code{\link{plot.apple}}, \code{\link{apple}} and \code{\link{predict.apple}}
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

%%%to do cross validation to apple 
fit.cv=cv.apple(X, y, family="binomial", alpha=0.25, K=5)

%%% plot paths
plot(fit.cv)
}