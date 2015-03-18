\name{plot.fanc}
\alias{plot.fanc}
\alias{MakeInterface}
\alias{cbExposeCanvas}
\alias{cbExposeLabel}
\alias{cbValueChanged}
\alias{cbValueChangedLabel}
\alias{cbExposeLabelGamma}
\alias{cbExposeLabelLambda}
\alias{cbValueChangedGamma}
\alias{cbValueChangedLabelGamma}
\alias{cbValueChangedLabelLambda}
\alias{cbValueChangedLambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot the solution path from a "fanc" object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This functions plots the solution paths from a "fanc" object for fixed value of gamma.  
}
\usage{
\method{plot}{fanc}(x, Window.Height=500, type=NULL, df.method="reparametrization", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{object} here~~
Fitted \code{"fanc"} model object.
}
  \item{Window.Height}{
%%     ~~Describe \code{object} here~~
A window height. The default is 500.
}
  \item{type}{
%%     ~~Describe \code{tuning} here~~
Two plot types are supported.  If \code{"path"}, the path diagram is depicted.  If \code{"heatmap"}, the heatmap is depicted.
}
  \item{df.method}{
Two types of degrees of freedom are supported.  If \code{"reparametrization"}, the degrees of freedom of the MC+ are reparametrized based on the degrees of freedom of the lasso.  If \code{"active"}, the degrees of freedom of are the number of nonzero parameters.
}
%  \item{gamma}{
%%     ~~Describe \code{tuning} here~~
%The value of gamma.
%}
  \item{\dots}{Other graphical parameters to plot}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
NULL
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

\seealso{\code{fanc} and \code{out} objects.}

%\examples{
%#generate data
%loadings0 <- matrix(c(rep(0.8,5),rep(0,5),rep(0,5),rep(0.8,5)),10,2)
%common.factors0 <- matrix(rnorm(50*2),50,2)
%unique.factors0 <- matrix(rnorm(50*10,sd=sqrt(0.36)),50,10)
%x <- common.factors0 \%*\% t(loadings0) + unique.factors0

%#fit data
%fit <- fanc(x,2)

%#plot solution path
%plot(fit)
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
