\name{pers.rate.obj}

\alias{pers.rate.obj}

\title{Calculating the global disease persistence rate in a metapopulation.}

\description{Calculating the global disease persistence rate in a metapopulation, by using the Kaplan–Meier survival curve of the metapopulation.}

\usage{pers.rate.obj<-function(object)}

\arguments{
  \item{object}{a seir object with the condition, this object has to have the value of the parameter 'persistence'. Because, based on this parameter, we can estimated the survival rate.} 
}

\details{For this function, we have a set of the survival time of all subpopulation in the metapopulation, due to the value of the parameter 'persistence'. Based on this value, we use the function 'survreg' of the package 'survival' to estimate the global persistence rate for a metapopulation. The function 'survival' is used for the Parametric Survival Model. Therefore, we get the estimated global persistence rate.
}

\value{
	A matrix (or vector) with columns giving lower and upper confidence limits for the estimated value. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%). 
}

\author{TRAN Thi Cam Giang}

\seealso{
	'pers.rate.obj' function in the 'dizzys' package.
}

\examples{
	#STO, STO
	sto<-seir(N=10e5,type="stoch",duration=30*365,beta1=0.1,nbVilles=10)
	objper<- persistence(sto)
	pers.obj<-pers.rate.obj(objper)
		pause()
}

\keyword{ seir model }
\keyword{ R package }% __ONLY ONE__ keyword per line
