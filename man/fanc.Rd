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
fanc(x, factors, n.obs, cor.factor=FALSE, normalize = TRUE, rho.max, covmat, 
control=list())}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix.
}
  \item{factors}{
The number of factors.
}
  \item{cor.factor}{
An indicator of the factor correlation.  If \code{"TRUE"}, the factor correlation is considered.  Default is \code{"FALSE"}.
}
  \item{normalize}{
If \code{"TRUE"}, each variable is normalized, otherwise it is left alone.
}
  \item{rho.max}{
Maximum value of rho.
}
  \item{covmat}{
A covariance matrix, which is needed if the data matrix \code{"x"} is not available.
}
  \item{n.obs}{
The number of observations, which is needed to calculate the model selection criteria and goodness-of-fit indices when the data matrix \code{"x"} is not available.
}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
}

\details{
The \code{control} argument is a list that can supply any of the following components:
\describe{
	 \item{\code{length.rho}}{
Candidates of tuning parameters which is used for grid search of reparametrization of MC+.
}
  \item{\code{length.gamma}}{
A length of tuning parameter which controls sparsenesses.  For each \code{rho}, \code{gamma=Inf} yields soft threshold operator (i.e., lasso penalty) and \code{gamma=+1} produces hard threshold operator.
}
  \item{\code{max.gamma}}{
A maximum value of gamma (excludes \code{Inf}.).
}
  \item{\code{min.gamma}}{
A minimum value of gamma.
}
  \item{\code{eta}}{
A tuning parameter used for preventing the occurrence of improper solutions.  \code{eta} must be non-negative.
}
  \item{\code{ncand.initial}}{
The number of candidates of initial values of factor loadings. 
}
  \item{\code{maxit.em}}{
A maximum number of iterations for EM algortihm. 
}
  \item{\code{maxit.cd}}{
A maximum number of iterations for coordinate descent algortihm. 
}
  \item{\code{maxit.bfgs}}{
A maximum number of iterations for BFGS algorithm used in the update of factor correlation. 
}
  \item{\code{maxit.initial}}{
A maximum number of iterations for choosing the initial values. 
}
  \item{\code{start}}{
	Type of start.  If \code{"cold"}, the initial value of factor loadings is randomly chosen for each tuning parameter, which can be slow.
}
  \item{\code{Delta}}{
A proportion of maximum value of rho to minimum value of rho, i.e., \code{rho.min}
\code{=Delta*rho.max}. 
}
  \item{\code{min.uniquevar}}{
A minimum value of unique variances.
}
  \item{\code{tol.em}}{
A positive scalar giving the tolerance at which the parameter in EM is considered close enough to zero to terminate the algorithm. 
}
  \item{\code{tol.cd}}{
A positive scalar giving the tolerance at which the factor loadings in coordinate descent is considered close enough to zero to terminate the algorithm. 
}
  \item{\code{tol.bfgs}}{
A positive scalar giving the tolerance at which the factor correlation in BFGS algorithm is considered close enough to zero to terminate the algorithm. 
}
   \item{\code{min.rhozero}}{
	If \code{"TRUE"}, the minimum value of \code{"rho"} is zero.
}
   \item{\code{zita}}{
	A value of hyper-parameter of factor correlation.
}
   \item{\code{progress}}{
	If \code{"TRUE"}, the progress for each tuning parameter is displayed.
}
   \item{\code{openmp}}{
	If \code{"TRUE"}, the parallel computation via OpenMP is excecuted.
}
   \item{\code{num.threads}}{
	The number of threads of the openmp. Only used when \code{openmp} is \code{"TRUE"}, 
}
   \item{\code{gamma.ebic}}{
	The value of gamma used in the extended BIC 
}

}
}
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
Hirose, K. and Yamamoto, M. (2014). 
\emph{Sparse estimation via nonconcave penalized likelihood in a factor analysis model},\cr
  \emph{Statistics and Computing, in press}
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

#select a model via model selection criterion
select(fit, criterion="BIC", gamma=Inf)

#plot solution path
#plot(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
