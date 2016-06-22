\name{seir}

\alias{seir}

\title{Creat a seir Object}

\description{This is important function, it creates a seir object for 'seir' class.}

\usage{
seir(type="stoch",duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/8,gamma=1/5,
	T=365,phi=0,nbVilles=1,epsilon=0.0,rho=0.0,seed=as.numeric(Sys.time()),rng="good",S=NULL,E=NULL,I=NULL,R=NULL,N=1e7)
}

\arguments{
  \item{type}{type of simulation, 'deterministic' or 'stochastic'.}
  \item{duration}{time values over which to perform the numerical integration.}
  \item{S}{initial value of the state variable 'S', this is the number of susceptible individuals at the time t0.}
  \item{E}{initial value of the state variable 'E', this is the number of exposed individuals at the time t0.}
  \item{I}{initial value of the state variable 'I', this is the number of infected individuals at the time t0.}
  \item{R}{initial value of the state variable 'R', this is the number of recovered individuals in immunity at the time t0.}
  \item{N}{number of population.}
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
  \item{seed}{random seed.}
  \item{rng}{random number generator, 'good' or 'fast'. If we chose rng="good", it means that we use the random number generator in C++. In the other side, we chose rng="fast", it means that we use the random number generator of the professor Yann Chevaleyre at university of Paris 13. With the random generator's Yann Chevaleyre, we find that it is faster than the random generator in C++, in contrast, its accuracy is less than that of the random generator's Yann Chevaleyre.}
  \item{method}{simulation algorithm,'direct' or 'adaptivetau'. If we chose method="direct", it means that we implement the direct algorithm of Gillespie in 1977. In contrast, if we chose method="adaptivetau", it means that we implement the adaptive tau-leaping approximation for simulating the trajectory of a continuous-time Markov process.}
  \item{unitTIME}{unit of time for simulation.}  
  \item{\dots}{further parameters.}
}

\details{To implement this function, we have three step. The first step, in R,  is to get the initial values of state variables and the values of parameter, then, sends those values to the C++ functions. The second step is to do simulation in the C++ implementation. The third step is to send the result from the C++ function to R. We will get the simulation resultat in R.

Moreover, in the set of arguments, there are some arguments we can give each subpopulation each value's.

  \bold{S, E, I, R, N}: the initial values of state variables for 'nbVilles' subpopulations. Multiple values can be specified so that each subpopulation can be given its own state variables. If there are fewer values than subpopulations they are recycled in the standard fashion. Subpopulations will all be simulated in the first value specified.

  \bold{mu, beta0, beta1, sigma, gamma, phi}: the parameters of simulation for 'nbVilles' subpopulations. Multiple values can be specified so that each subpopulation can be given its own parameters. If there are fewer values than subpopulations they are recycled in the standard fashion. Subpopulations will all be simulated in the first value specified.
}

\value{Result returned is a seir object of the 'seir' class with all values for all slots.}

\references{Norman Matloff (2009), The Art of R Programming.}

\author{TRAN Thi Cam Giang}

\seealso{
	'seir', 'seir.det', 'seir.stoch' fuctions in 'dizzys' package.
}

\examples{
	obj<-seir()
	plot(obj)
	obj<-seir(type="sto",nbVilles=3,N=c(1e7,1e6))
	plot(obj)
	obj<-seir(type="deter",nbVilles=3,N=c(1e7,1e6))
	plot(obj)	
}

\keyword{ deterministic model }
\keyword{ stochastic model }% __ONLY ONE__ keyword per line
