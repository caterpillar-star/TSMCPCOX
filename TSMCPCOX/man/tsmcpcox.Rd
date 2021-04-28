\name{tsmcpcox}
\alias{tsmcpcox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two stage multiple change point detection in the Cox model.
}
\description{
This function provides a two stage procedure for simultaneously detecting multiple change points in the Cox model. In the splitting stage, the change point problem is converted into a model selection problem so that a modern model selection method can be applied. In the refining stage, the change points obtained in the splitting stage are finalized via a refining method. The tuning parameter lambda
is chosen by BIC.
}




\usage{
tsmcpcox(Y, Delta, X1, X2 = NULL, Z, method = c("MCP", "SCAD"), c = 25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
the right-censored failure time.
}
  \item{Delta}{
%%     ~~Describe \code{Delta} here~~
}
  \item{X1}{
the n-by-p design matrix which has the threshold effect.
}
  \item{X2}{
the n-by-p2 design matrix which doesn't has the threshold effect.
}
  \item{Z}{
the covariate where the change points are.
}
  \item{method}{
the method to be used by  mcp or scad.
}
  \item{c}{
The number of candidates for the segment length.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
tsmcpcox returns an object of class "tsmcpcox".
An object of class "tsmcpcox" is a list containing the following components:
\item{the number of change.points}{estimators of change points.}
\item{change.points}{estimators of change points.}
\item{coefficients}{estimators of coefficients.}
\item{m.opt}{the selected segment length.}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Bing Wang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
survival ncvreg SIS NMOF
}
\examples{
#true change point location:
## 0.3,0.7
library(TSMCPCOX)
 library(survival)
 library(NMOF)
library(ncvreg)
library(SIS)
set.seed(10123)
n=1000
u<-runif(n,0,1)
 X<-cbind(rbinom(n,1,0.5),rnorm(n,0,1))
 Z=sort(runif(n,0,1))
 C<-rexp(n,0.5)  #censoring time
 tau<-c(0.3,0.7) ## true change point location
beta<-c(1,1) #cofficients in first segment
 b1<-c(-2,0)  #cofficients in second segment
 b2<-c(3,-3)  #cofficients in third segment
 r<-drop(X\%*\%beta)+drop(X\%*\%b1)*(Z>tau[1])+drop(X\%*\%b2)*(Z>tau[2])
 T<--log(u)/exp(r)# set lambda0(t)=1
 Y<-ifelse(T<=C,T,C)
 Delta<-ifelse(T<=C,1,0)

##estimate change points
#   tsmcpcox(Y=Y, Delta=Delta,X1=X,Z=Z, method = "MCP")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
