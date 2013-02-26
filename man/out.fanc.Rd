\name{out}
\alias{out}
\alias{out.aic}
\alias{out.bic}
\alias{out.caic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
output from a "fanc" object for fixed value of gamma.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This functions give us the loadings from a "fanc" object for fixed value of gamma.  
}
\usage{
out(x, rho, gamma, scores=FALSE)
out.aic(x, gamma, scores=FALSE)
out.bic(x, gamma, scores=FALSE)
out.caic(x, gamma, scores=FALSE)
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
  \item{gamma}{
%%     ~~Describe \code{tuning} here~~
The value of gamma.
}
  \item{rho}{
%%     ~~Describe \code{tuning} here~~
The value of rho.
}
  \item{scores}{
%%     ~~Describe \code{tuning} here~~
Logical flag for outputting the factor scores. Defalut is FALSE.
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
