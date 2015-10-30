\name{plot.comp.surv.estim}

\alias{plot.comp.surv.estim}


\title{Plotting Survival Probability upon Time}

\description{Plotting Kaplan–Meier survival curve based on persistence value of seir object and curve of survival function estimated by calculating cumulative hazard function.}

\usage{plot.comp.surv.estim(this, ...)}

\arguments{
  \item{this}{a seir object with the value of the slot "persistence".}
  \item{\dots}{other arguments.}
}

\details{To create Kaplan-Meier survival curve based on persistence value of seir object in the following:
First, we get the persistence vlue of the seir object. Then, we use the estimation method of KM to find the survival values S(t) depending on the time t, we obtain the survival curve. Second, based on the values S(t) by time t, we calculate the value for cumulative hazard function h(t) with h(t) = -log(s(t)), so we obtain the cumulative hazard curve h(t). Third, from the values of cumulative hazard h(t) depending in time t, we can find the formula h(t) with h(t) = h0*exp(beta*t) by using the COX model or COX regression. So now, we have the formula h(t) with the value of h0 and beta. Finally, according to the formula h(t), we can find the formula for s(t), with s(t) = exp(-h(t)). We call s(t) the estimated survival function. To perform this fucntion, we use the 'survfit' function of the 'survival' package, and 'coxph' , 'Surv' functions of the 'KMsurv' package.
}

\value{Result returned is  a plot of the Kaplan–Meier survival curve and the curve of estimated survival function.}

\references{David Diez, Survival Analysis in R.
	     Alex Cook (9/2008), Cox proportional hazards model, Introduction to Survival Analysis.
}

\author{TRAN Thi Cam Giang}

\seealso{'survival' package and 'KMsurv' package}

\examples{
	p<-persistence(seir(type="sto",nbVilles=15,N=1e5))	
	plot.comp.surv.estim(p)
}

\keyword{ persistence }
\keyword{ Cox regression }% __ONLY ONE__ keyword per line
