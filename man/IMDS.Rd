\name{IMDS}
\alias{IMDS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Multidimensional scaling of interval-valued dissimilarities.
}
\description{
	Performs MDS for given interval-valued dissimilarities.
}
\usage{
IMDS(IDM, p=2,eps= 1e-5 ,maxit =1000,model=c("sphere","box"),
opt.method=c("MM", "BFGS"), ini = "auto",report=100,grad.num=FALSE,
rel=0, dil=1)
}
\arguments{
		\item{IDM}{
		The interval-valued dissimilarity matrix (an object of class "array": \code{IDM[1,,]}: the lower dissmilarity matrix; \code{IDM[2,,]}: the upper dissmilarity matrix). 
	}
		\item{p}{
		Number of dimensions.
	}
		\item{eps}{
		Convergence criterion for the majorization minimization algorithm.
	}
		\item{maxit}{
		Maximum number of iteretions.
	}
		\item{model}{
		If "sphere", then the hypersphere model is used.
		If "box", then the hyperbox model is used.
	}
		\item{opt.method}{
		If "BFGS", then the BFGS method is used for optimizing the stress function.
		If "MM", then the majorization minimization algortihm is used.
	}
		\item{ini}{
		List which consists of an initial center coordinate matrix \code{ini[[1]]} and an initial radius vector \code{ini[[2]]} (optional). If "auto", then a configuration of the classical MDS for \code{(IDM[1,,]+IDM[2,,])/2} is used as the center coordinate matrix and radii are drawn from the uniform distribution \eqn{U(0,1)}.
	}
		\item{report}{
		The frequency of reports. Defaults to every 100 iterations.
	}
		\item{grad.num}{
		If \code{FALSE}, then exact gradient function is used in the BFGS method.
		If \code{TRUE}, then a numerical gradient is used in the BFGS method.
	}
		\item{rel}{
		If \code{acc=1}, accelerate by the relaxed update.
		If \code{acc=0}, the relaxed update is not used.
		Not used for the BFGS method. Defaults to \code{0}.
	}
		\item{dil}{
		If \code{acc=1}, accelerate by the optimal dilation of the configration.
		If \code{acc=0}, the optimal dilation is not used.
		Not used for the BFGS method. Defaults to \code{1}.
	}
}

\details{
  The default optimization method is a majorization-minimization algorithm with the optimal dilation. 
  
  Method \code{"MM"} is a majorization-minimization (MM) algortihm for the specified model. 
  If \code{model="box"}, method \code{"MM"} is a MM algorithm, called I-Scal, which is proposed by Groenen et al. (2006).
  If \code{model="sphere"}, method \code{"MM"} is a MM algorithm which can be considered as I-Scal for the hypersphere model.

  Method \code{"BFGS"} is a quasi-Newton method (also known as a variable
  metric algorithm), specifically that published simultaneously in 1970
  by Broyden, Fletcher, Goldfarb and Shanno. For more details, see Chapter 15 of Nash (1990).
}

\value{
\code{IMDS} returns a list with components:
\item{X }{ The best corrdinate matrix with p columns whose rows give the coordinates of the vertexes.}
\describe{
	\item{If \code{model="sphere"},}{\code{r } The best radius vector.}
	\item{If \code{model="box"},}{\code{R } The best radius matrix with p columns whose rows give the radii of objects.}
}
\item{str}{The value of the stress function of IMDS corresponding to \code{X} is returned.}
\item{str.vec}{If "MM", then the vector of values on each iteration is returned.}
\item{EIDM}{If "MM", then the interval-valued dissimilarity matrix correspondint to the estimated parameters.}
}
\references{
	 Groenen, P. J. F., Winsberg, S., Rodriguez, O., and Diday, E. (2006). I- scal: Multidimensional scaling of interval dissimilarities. \emph{Computational Statistics & Data Analysis}, \bold{51}, 360--378.
	 
	 Nash, J. C. (1990) \emph{Compact Numerical Methods for
    Computers. Linear Algebra and Function Minimisation.} Adam Hilger.
}
\author{
Yoshikazu Terada
}
\examples{
###################
#Fats and Oil data#
###################
###############################
data(oil.idiss)
#Apply the hypersphere model via the BFGS method
set.seed(1)
res.bfgs <- IMDS(IDM=oil.idiss, p=2,model="sphere",opt.method="BFGS", ini = "auto")
plot(res.bfgs,main="Sph_bfgs")
#Apply the hypersphere model via the MM algorithm
set.seed(1)
res.mm <- IMDS(IDM=oil.idiss, p=2,model="sphere",opt.method="MM", ini = "auto")
plot(res.mm,main="Sph_MM")
#Apply the hyperbox model via the BFGS method
set.seed(1)
res.bfgs <- IMDS(IDM=oil.idiss, p=2,model="box",opt.method="BFGS", ini = "auto")
plot(res.bfgs,main="Box_bfgs")
#Apply the hyperbox model via the MM algorithm
set.seed(1)
res.mm <- IMDS(IDM=oil.idiss, p=2,model="box",opt.method="MM", ini = "auto")
plot(res.mm,main="Box_MM")
###############################
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

