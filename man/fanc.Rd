\name{fanc}
\alias{fanc}
\alias{print.fanc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fanc (penalized maximum likelihood factor analysis via nonconvex penalties)
}
\description{
This package computes the solution path of penalized maximum likelihood estimates via MC+ penalties.  
}
\usage{
fanc(x, factors, covmat, n.obs, cor.factor=FALSE, length.rho=20, rho.max, length.gamma=5, 
max.gamma=20, min.gamma=1.1, eta=0.0, initial.iter=10,  max.count.initial=500,  
max.count.em=10000, max.count.cd=500, max.count.bfgs=500, all.random.start=FALSE, 
Delta=0.001, min.uniquevar=0.005, tol.em=1e-5, tol.cd=1e-5, tol.bfgs=1e-5, scale.x = TRUE, 
min.rho.zero=FALSE, zita=0.0, trace.em=FALSE)}
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
  \item{cor.factor}{
An indicator of the factor correlation.  If \code{"TRUE"}, the factor correlation is considered.
}
  \item{length.rho}{
Candidates of tuning parameters which is used for grid search of reparametrization of MC+.
}
  \item{rho.max}{
Maximum value of rho.
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
  \item{max.count.bfgs}{
A maximum number of iterations for BFGS algorithm used in the update of factor correlation. 
}
  \item{max.count.initial}{
A maximum number of iterations for choosing the initial values. 
}
  \item{all.random.start}{
	If \code{"TRUE"}, the initial value of factor loadings is randomly chosen for each tuning parameter, which can be very slow.
}
  \item{Delta}{
A proportion of maximum value of rho to minimum value of rho, i.e., \code{rho.min}
\code{=Delta*rho.max}. 
}
  \item{min.uniquevar}{
A minimum value of unique variances.
}
  \item{tol.em}{
A positive scalar giving the tolerance at which the parameter in EM is considered close enough to zero to terminate the algorithm. 
}
  \item{tol.cd}{
A positive scalar giving the tolerance at which the factor loadings in coordinate descent is considered close enough to zero to terminate the algorithm. 
}
  \item{tol.bfgs}{
A positive scalar giving the tolerance at which the factor correlation in BFGS algorithm is considered close enough to zero to terminate the algorithm. 
}
  \item{scale.x}{
If \code{"TRUE"}, each variable is standardized, otherwise it is left alone.
}
   \item{min.rho.zero}{
	If \code{"TRUE"}, the minimum value of \code{"rho"} is zero.
}
  \item{initial.iter}{
The number of candidates of initial values of factor loadings. 
}
   \item{zita}{
	A value of hyper-parameter of factor correlation.
}
   \item{trace.em}{
	If \code{"TRUE"}, tracing information on the progress for each tuning parameter is displayed.
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{loadings}{factor loadings}
\item{uniquenesses}{unique variances}
\item{Phi}{factor correlation}
\item{rho}{rho}
\item{AIC}{AIC}
\item{BIC}{BIC}
\item{CAIC}{CAIC}
\item{df}{degrees of freedom (number of non-zero parameters for the lasso estimation)}
\item{criteria}{values of AIC, BIC and CAIC}
\item{goodness.of.fit}{values of GFI and AGFI}
\item{gamma}{a value of gamma}
\item{Npflag}{If the number of observation is larger than the number of variables, 1, otherwise 0.}
\item{factors}{the number of factors}
\item{cor.factor}{An indicator of the factor correlation}
\item{x}{data matrix}
\item{convergence}{indicator of convergence of EM algorithm, coordinate descent and BFGS.  If all of these variables are 0, the algorithm has been converged}
}

\references{
%% ~put references to the literature/web site here ~
%Friedman, J. (2008).  Fast sparse regression and classification. \eqn{Technical report}, Standford University.\cr
Hirose, K. and Yamamoto, M. (2012). 
\emph{Sparse estimation via nonconcave penalized likelihood in a factor analysis model},\cr
  \url{http://arxiv.org/abs/1205.5868}\cr
  \emph{arXiv:1205.5868 (arXiv).}
}
\author{
Kei Hirose\cr
\email{mail@keihirose.com}
}

%\note{
%  ~~further notes~~
%}

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{out} and \code{plot.fanc} objects.}
\examples{
#generate data
set.seed(0)
loadings0 <- matrix(c(rep(0.8,5),rep(0,5),rep(0,5),rep(0.8,5)),10,2)
common.factors0 <- matrix(rnorm(50*2),50,2)
unique.factors0 <- matrix(rnorm(50*10,sd=sqrt(0.36)),50,10)
x <- common.factors0 \%*\% t(loadings0) + unique.factors0

#fit data
fit <- fanc(x,2)
fit2 <- fanc(x,2,cor.factor=TRUE) #factor correlation is estimated

#print candidates of gamma and rho
print(fit)

#output for fixed tuning parameters
out(fit, rho=0.1, gamma=Inf)

#plot solution path
#plot(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
