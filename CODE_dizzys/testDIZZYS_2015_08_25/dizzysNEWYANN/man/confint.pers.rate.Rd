\name{confint.pers.rate}

\alias{confint.pers.rate}

\title{Confidence interval for the estimated global disease persistence rate in a metapopulation.}

\description{This function permits us to find the confidence interval for  the estimated global disease persistence rate in a metapopulation.}

\usage{confint.pers.rate<-function(object, level)}

\arguments{
  \item{object}{a seir object with the condition, this object has to have the value of the parameter 'persistence'. Because, based on this parameter, we can estimated the survival rate.} 
  \item{level}{the confidence level required.}
}

\details{For this function, in the step 1, we also use the function 'survreg' of the package 'survival' to estimate the global persistence rate for a metapopulation. Then, we use the generic function 'confint' to find the confidence interval for the estimated global disease persistence rate in the metapopulation.
}

\value{
	Result returned is a number.
}

\author{TRAN Thi Cam Giang}

\seealso{
	'confint.pers.rate' function in the 'dizzys' package.
}

\examples{
	#STO, STO
	sto<-seir(N=10e5,type="stoch",duration=30*365,beta1=0.1,nbVilles=10)
	objper<- persistence(sto)
	pers.obj<-pers.rate.obj(objper)
	confint.pers.rate(pers.obj)
	pause()
}

\keyword{ seir model }
\keyword{ R package }% __ONLY ONE__ keyword per line
