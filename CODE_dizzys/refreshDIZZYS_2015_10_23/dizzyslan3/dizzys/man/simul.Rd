\name{simul}

\alias{simul}

\title{Redoing or Continuing a Simulation.}

\description{Redoing or continuing a simulation by using values of parameters given.}

\usage{
simul(object,type="deter",continue=F,duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/8,gamma=1/5,
T=365,phi=0,nbVilles=1,epsilon=0.0,rho=0.0,seed=as.numeric(Sys.time()),rng="good",append=TRUE,t0=NULL,S=NULL,E=NULL,I=NULL,R=NULL,N=1e7)
}

\arguments{
  \item{object}{a seir object. If missing this object, we create a new seir object by using the values of state variables and of parameters given.}
  \item{type}{type of simulation, 'deterministic' or 'stochastic'.}
  \item{continue}{logical (T/F). If continue=F, it means to redo silulation. Else, if continue=T, it means to continue to do silulation.}
  \item{duration}{time values over which to perform the numerical integration.}
  \item{method}{simulation algorithm,'direct' or 'adaptivetau'. If we chose method="direct", it means that we implement the direct algorithm of Gillespie in 1977. In contrast, if we chose method="adaptivetau", it means that we implement the adaptive tau-leaping approximation for simulating the trajectory of a continuous-time Markov process.}
  \item{unitTIME}{unit of time.}
  \item{mu}{per capita birth and death rates (per day).}
  \item{beta0}{mean value of contact rate (per day).}
  \item{beta1}{amplitude of contact rate.}
  \item{sigma}{transition rate from exposed (E) to infected (I), per day (inverse of the average duration of the latency period).}
  \item{gamma}{recovery rate (per day).}
  \item{T}{period of the contact rate (in days).}
  \item{phi}{phase of the contact rate (in radians).}
  \item{nbVilles}{number of subpopulations in the metapopulation.}
  \item{epsilon}{proportion of infections by contacts with infected from other subpopulations. The value of epsilon between subpopulation i and j is equal to the value of epsilon between subpopulation j and i. epsilon is always greater than or equal to zero and less than or equal to one,specially, epsilon, of the subpopulation i and i, is equal to zero.}
  \item{rho}{coupling rate between subpopulations i and j. rho between subpopulation i and j is equal to rho between subpopulation j and i. rho is always greater than or equal to zero and less than or equal to one, specially, rho, of the subpopulation i and i, is equal to one.}
  \item{seed}{random seed for random number generator.}
  \item{rng}{random number generator, 'good' or 'fast'. If we chose rng="good", it means that we use the random number generator in C++. In the other side, we choose rng="fast", it means that we use the random number generator of the professor Yann Chevaleyre at university of Paris 13. With the random generator's Yann Chevaleyre, we find that it is faster than the random generator in C++, in contrast, its accuracy is less than that of the random generator's Yann Chevaleyre.}
  \item{S}{initial value of the state variable 'S', this is the number of susceptible individuals at the time t0.}
  \item{E}{initial value of the state variable 'E', this is the number of exposed individuals at the time t0.}
  \item{I}{initial value of the state variable 'I', this is the number of infected individuals at the time t0.}
  \item{R}{initial value of the state variable 'R', this is the number of recovered individuals in immunity at the time t0.}
  \item{N}{number of population.}
  \item{append}{logical(T/F). It is available when continue=T.
	If continue=T and append=T, it means that we want our new object contains the old data + the new.
        If continue=T and append=F, we only get the new.}
  \item{t0}{numeric. It means that if we want the start time after what has already been simulated ("t0 = NULL") or if we prefer that time restarts any value (for example "t0 = 1").}
}

\details{For this fucntion, to get the values of arguments, first, for the arguments missed, we get its values from the seir object, in contrast, we get its given values. Moreover, in the set of arguments, there are some arguments we can give each subpopulation each value's.
  \bold{S, E, I, R, N}: the initial values of state variables for 'nbVilles' subpopulations. Multiple values can be specified so that each subpopulation can be given its own state variables. If there are fewer values than subpopulations they are recycled in the standard fashion. Subpopulations will all be simulated in the first value specified.
  \bold{mu, beta0, beta1, sigma, gamma, phi}: the parameters of simulation for 'nbVilles' subpopulations. Multiple values can be specified so that each subpopulation can be given its own parameters. If there are fewer values than subpopulations they are recycled in the standard fashion. Subpopulations will all be simulated in the first value specified.
}

\value{
	Result returned is an object of the 'seir' class.
}

\author{TRAN Thi Cam Giang}

\seealso{
	'simul' function in the 'dizzys' package.
}

\examples{
	#STO, STO
	sto<-seir(N=10e5,type="stoch",duration=5*365,beta1=0.1,nbVilles=2)
	plot(simul(sto,type="stoch",continue=TRUE,duration=5*365,beta1=0,phi=c(pi/2,0)))
		pause()
	#DET, DET
	det<-seir(N=10e5,type="deter",duration=50*365)
	plot(simul(det,type="deter",continue=TRUE,duration=50*365,beta1=0,phi=pi/2))
		pause()
}

\keyword{ seir model }
\keyword{ R package }% __ONLY ONE__ keyword per line
