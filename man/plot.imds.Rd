\name{plot.imds}
\alias{plot.imds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	plotting the estimated spheres or boxes.
%%  ~~function to do ... ~~
}
\description{
	plotting the estimated hyperspheres or hyperboxes.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{plot}{imds}(x, xylim="auto",clab=1:nrow(X),lab.cex=1,lab.col="black",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  	Object class "imds"
}
  \item{xylim}{
  	The x limits (\code{xylim[1,1]}, \code{xylim[1,2]}) of the plot. 
  	The y limits (\code{xylim[2,1]}, \code{xylim[2,2]}) of the plot. 
  	Note that \code{xylim[i,1]} > \code{xylim[i,2]} (i=1, 2) is allowed and leads to a 'reversed axis'.
}
  \item{clab}{
  	A character vector for objects.
}
  \item{lab.cex}{
  	numeric character expansion factor.
}
  \item{lab.col}{
  	label color.
}
  \item{...}{
  	further graphical parameters (from par).
}
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
