\name{equilibrium}

\alias{equilibrium}

\title{Finding Endemic Equilibrium of a seasonally-forced SEIR/SIR model}

\description{'equilibrium' returns the values of the state variables of an SEIR/SIR model at endemic equilibrium.}

\usage{equilibrium(duration, unitTIME, N, mu, beta0, beta1, sigma, gamma, phi, T)}

\arguments{
  \item{duration}{time values over which to perform the numerical integration.}
  \item{unitTIME}{time unit of simulation.}
  \item{N}{number of population.}
  \item{mu}{per capita birth and death rates (per day).}
  \item{beta0}{mean value of contact rate (per day).}
  \item{beta1}{amplitude of contact rate.}
  \item{sigma}{transition rate from exposed (E) to infected (I), per day (inverse of the average duration of the latency period.}
  \item{gamma}{recovery rate (per day).}
  \item{phi}{phase of the contact rate (in radians).}
  \item{T}{period of the contact rate (in days).}
}
\details{The host population is supposed to be at demographic equilibrium with the births balancing the deaths (births and death occur with the same rate 'mu') and a total population size constant and equal to 1. The contact rate is forced by a sinusoid of period 'T' and phase 'phi'. The SEIR epidemiological model is defined by the following set of differential equations:
	
	dS/dt =  mu - beta0 * (1 + beta1 * cos(2 * pi * t / T + phi)) * S * I - mu * S

	dE/dt = beta0 * (1 + beta1 * cos(2 * pi * t / T + phi)) * S * I - mu * E - sigma * E

	dI/dt = sigma *  E - mu * I - gamma * I

	dR/dt = gamma * I - mu * R
	

'beta0' is the mean value of the contact rate and 'beta1' is the amplitude (in percentage of the mean) around the mean. The endemic equilibrium value of this sinusoidally-forced epidemiological system is found in two steps. First the endemic equilibium point (S*,E*,I*,R*) of the unforced system ('beta1 = 0') is found from the following analytical equations:

	S* = (gamma + mu)*(sigma + mu)/(beta0 * sigma)

	E* = mu*((1/(sigma + mu))-((gamma + mu)/(beta0 * sigma)))

	I* = mu*((beta0*sigma - (gamma + mu)*(sigma + mu))/(beta0*(gamma + mu)*(sigma + mu)))

	R* = 1 - S* - E* - I*

This endemic equilibium point is used as a starting point for numerical integration of the system of differential equations. This is done by the 'seir.det' function that uses the 'ode' function of the 'deSolve' package. The last values of the of the state variable are returned. These correspond to a given point on the equilibrium limit cycle (assuming that this equilibrium limit cycle is reached) and the value of this point will depend on the phase 'phi'.}

\value{'equilibrium' returns a named vector of 3 numerical proportions containing the values of the state variables S, E and I on one point of the equilibrium limit cycle (provided that this limit cycle is reached). In case this limit cycle is reached, the value of the S, E and I state variables will depend on the input phase value 'phi'.}

\author{TRAN Thi Cam Giang}

\references{Anderson RM & May RM (1991) Infectious Diseases of Humans - Dynamics and Control. Oxford University Press.}

\seealso{'seir' in package 'dizzyz' and 'ode' in package 'deSolve'.}

\examples{
## The point on the limit cycle depends on the input phase value 'phi':
	equilibrium(mu=.000456, beta0=1000/365, beta1=.1, sigma=1/7, gamma=1/7, phi=0, T=365,duration=100,unitTIME=1)
	equilibrium(mu=.000456, beta0=1000/365, beta1=.1, sigma=1/7, gamma=1/7, phi=pi/2, T=365, duration=100,unitTIME=7)
}

\keyword{SEIR/SIR model}
\keyword{limit cycle equilibrium}% __ONLY ONE__ keyword per line
