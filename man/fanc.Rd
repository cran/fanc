\name{fanc}
\alias{fanc}
\alias{print.fanc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fanc (penalized maximum likelihood Factor Analysis via Non-Concave penalties)
}
\description{
This package computes the solution path of penalized maximum likelihood estimates via MC+ penalties.  
}
\usage{
fanc(x, factors, covmat, n.obs, length.rho=20, length.gamma=5, 
max.gamma=20, min.gamma=1.1, eta=0.0, max.count.em=2000, max.count.cd=2000, 
Delta=0.001, min.uniquevar=0.005, tol1=1e-5, 
tol2=1e-5, scale.x = TRUE, min.rho.zero = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix.
}
  \item{factors}{
The number of factors.
}
  \item{covmat}{
A covariance matrix, which is needed if the data matrix \code{"x"} is not available.
}
  \item{n.obs}{
The number of observations
}
  \item{length.rho}{
Candidates of tuning parameters which is used for grid search of reparametrization of MC+.
}
  \item{length.gamma}{
A length of tuning parameter which controls sparsenesses.  For each \code{rho}, \code{gamma=Inf} yields soft threshold operator (i.e., lasso penalty) and \code{gamma=+1} produces hard threshold operator.
}
  \item{max.gamma}{
A maximum value of gamma (excludes \code{Inf}.).
}
  \item{min.gamma}{
A minimum value of gamma.
}
  \item{eta}{
A tuning parameter used for preventing the occurrence of improper solutions.  \code{eta} must be non-negative.
}
  \item{max.count.em}{
A maximum number of iterations for EM algortihm. 
}
  \item{max.count.cd}{
A maximum number of iterations for coordinate descent algortihm. 
}
%  \item{max.count.initial}{
%A maximum number of iterations for choosing the initial values. 
%}
  \item{Delta}{
A proportion of maximum value of rho to minimum value of rho, i.e., \code{rho.min}
\code{=Delta*rho.max}. 
}
  \item{min.uniquevar}{
A minimum value of unique variances.
}
  \item{tol1}{
A positive scalar giving the tolerance at which the factor loadings in EM is considered close enough to zero to terminate the algorithm. 
}
  \item{tol2}{
A positive scalar giving the tolerance at which the factor loadings in coordinate descent is considered close enough to zero to terminate the algorithm. 
}
  \item{scale.x}{
If \code{"TRUE"}, each variable is standardized, otherwise it is left alone.
}
   \item{min.rho.zero}{
	If \code{"TRUE"}, the minimum value of \code{"rho"} is zero.
}
%  \item{initial.iter}{
%The number of candidates of initial values of factor loadings. 
%}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%Friedman, J. (2008).  Fast sparse regression and classification. \eqn{Technical report}, Standford University.\cr
%Hirose, K., Tateishi, S. and Konishi, S.. (2011).  Efficient algorithm to select tuning parameters in sparse regression modeling with regularization. arXiv:1109.2411 (arXiv).
%}
\author{
Kei Hirose\cr
\email{mail@keihirose.com}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{out} and \code{plot.fanc} objects.}
\examples{
#generate data
loadings0 <- matrix(c(rep(0.8,5),rep(0,5),rep(0,5),rep(0.8,5)),10,2)
common.factors0 <- matrix(rnorm(50*2),50,2)
unique.factors0 <- matrix(rnorm(50*10,sd=sqrt(0.36)),50,10)
x <- common.factors0 \%*\% t(loadings0) + unique.factors0

#fit data
fit <- fanc(x,2)

#output for fixed tuning parameters
out(fit, rho=0.1, gamma=Inf)

#plot solution path
plot(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
