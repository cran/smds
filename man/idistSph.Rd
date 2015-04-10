\name{idistSph}
\alias{idistSph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute the interval distance of the hypersphere model.
%%  ~~function to do ... ~~
}
\description{
	Compute the interval distance for given hypersphere objects.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
idistSph(X,r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  	The center coordinate matrix
}
  \item{r}{
  	The radius vector.
}
}
\value{
\code{idistSph} returns the interval-valued dissimilarity matrix \code{IDM} (an object of class "array": \code{IDM[1,,]}: the lower dissmilarity matrix; \code{IDM[2,,]}: the upper dissmilarity matrix). 
}
\author{
	Yoshikazu Terada
%%  ~~who you are~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
