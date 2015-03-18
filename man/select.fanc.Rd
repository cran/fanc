\name{select}
\alias{select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
select from a "fanc" object for fixed value of gamma.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This functions give us the loadings from a "fanc" object for fixed value of gamma.  
}
\usage{
select(x, criterion=c("BIC","AIC","CAIC","EBIC"), 
	gamma,  scores=FALSE, df.method="active")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{object} here~~
Fitted \code{"fanc"} model object.
}
%%  \item{type}{
%%     ~~Describe \code{tuning} here~~
%%The type of \code{"solution"} depicts the solution path along with the tuning parameter selected by model selection criterion of \code{criteion}. The type of \code{"df"} shows the degrees of freedom.
%%}
  \item{criterion}{
%%     ~~Describe \code{tuning} here~~	
The criterion by which to select the tuning parameter rho. One of "AIC", "BIC", "CAIC", or "EBIC". Default is "BIC".
}
  \item{gamma}{
%%     ~~Describe \code{tuning} here~~
The value of gamma.
}
  \item{scores}{
%%     ~~Describe \code{tuning} here~~
Logical flag for outputting the factor scores. Defalut is FALSE.
}
  \item{df.method}{
Two types of degrees of freedom are supported.  If \code{"active"}, the degrees of freedom of are the number of nonzero parameters.  If \code{"reparametrization"}, the degrees of freedom of the MC+ are reparametrized based on the degrees of freedom of the lasso.  
}
%  \item{\dots}{Other graphical parameters to plot}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{loadings}{factor loadings}
\item{uniquenesses}{unique variances}
\item{Phi}{factor correlation}
\item{scores}{factor scores}
\item{df}{degrees of freedom (number of non-zero parameters for the lasso estimation)}
\item{criteria}{values of AIC, BIC and CAIC}
\item{goodness.of.fit}{values of GFI and AGFI}
\item{rho}{a value of rho}
\item{gamma}{a value of gamma}
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
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{fanc} and \code{plot.fanc} objects.}

%\examples{
%#generate data
%loadings0 <- matrix(c(rep(0.8,5),rep(0,5),rep(0,5),rep(0.8,5)),10,2)
%common.factors0 <- matrix(rnorm(50*2),50,2)
%unique.factors0 <- matrix(rnorm(50*10,sd=sqrt(0.36)),50,10)
%x <- common.factors0 \%*\% t(loadings0) + unique.factors0

%#fit data
%fit <- fanc(x,2)

%#output for fixed tuning parameters
%out(fit, rho=0.1, gamma=Inf)
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
