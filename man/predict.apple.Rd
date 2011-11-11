\name{predict.apple}
\alias{predict.apple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model prediction based on a fitted apple/cv.apple object.
}
\description{
Similar to other predict methods, this function returns predictions from a fitted \code{"apple"} or \code{"cv.apple"} object.
}
\usage{
\method{predict}{apple}(object, X, which = 1:length(object$lambda), 
type = c("link", "response", "class"),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
fitted \code{"apple"} or \code{"cv.apple"} model object.
}
  \item{X}{
matrix of values at which predictions are to be made.
}
  \item{which}{
indices of the penalty parameter \code{lambda} at which predictions are required. by default, all indices are returned.
}
  \item{type}{
type of prediction: \code{"link"} returns the linear predictors; \code{"response"} gives the fitted values; \code{"class"} returns the binomial outcome with the highest probability.
}
  \item{...}{
see \code{matplot}.
}
}
\value{
The object returned depends on type.
}
\references{
Yi Yu and Yang Feng, APPLE: Approximate Path for Penalized Likelihood Estimator, manuscript.
}
\author{
Yi Yu and Yang Feng
}

\seealso{
\code{\link{apple}}, \code{\link{cv.apple}} and \code{\link{plot.apple}}
}
\examples{
%%%generate training data
p=10
n=200
d=5
coefs=c(3,1.5,0,0,2)
intercept=0
beta=rep(0,p)
beta[1:d]=coefs
set.seed(2)
X=matrix(rnorm(p*n), nrow=n)
mu=1/(1+exp(-X \%*\% beta-intercept))
y=rbinom(n,1,mu)
	
%%%fit the model
fit.apple=apple(X, y, family="binomial")

%%% generate test data
set.seed(3)
testX=matrix(rnorm(p*n), nrow=n)

predict(fit.apple,testX,type="link")
predict(fit.apple,testX,type="response")
predict(fit.apple,testX,type="class")

%%% predictor for cv.apple object
fit=cv.apple(X, y, family="binomial", alpha=0)
predict(fit.apple,testX,type="link", which = fit$cv.loc)
predict(fit.apple,testX,type="response", which = fit$cv.loc)
predict(fit.apple,testX,type="class", which = fit$cv.loc)

}