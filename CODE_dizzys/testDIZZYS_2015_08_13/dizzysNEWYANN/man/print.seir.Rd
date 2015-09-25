\name{print.seir}

\alias{print.seir}

\title{Printing  seir Object}

\description{Generic function allowing us to print a seir object.}

\usage{print.seir(object, ...)}

\arguments{
  \item{object}{a seir object.}
  \item{\dots}{further arguments.}
}

\details{The function 'print' allows us to print all slots of a seir object of the 'seir' class. It contains  in the following:

	'pop' is a list of the values of the state variales according to time.

	'duration' is time to do simulation.

	'S, E, I, R'  are the initial values of the state variables, number of susceptible, exposed, infected, recovered individuals, respectively.

	'N' is the initial number of population.

	'T' is period of the contact rate (in days).

	'mu' is per capita birth and death rates (per day).

	'beta0' is mean value of contact rate (per day).

	'beta1' is amplitude of contact rate.

	'sigma' is transition rate from exposed (E) to infected (I), per day (inverse of the average duration of the latency period.

	'gamma' is recovery rate (per day).

	'unitTIME' is unit of time.

	'phi' is phase of the contact rate (in radians).

	'nbVilles' is number of subpopulation in metapopulation.

	'type' is type of simulation, 'deterministic' or 'stochastic'.

	'epsilon' is proportion of infections by contacts with infected from other subpopulations.

	'rho' is coupling rate between subpopulations i and j.

	'seed' is seed number of random number generator.

	'rng' is random number generator, 'good' or 'fast'.

	'method' is simulation algorithm,'direct' or 'adaptivetau'.

	'persistence' is data frame about the information of the persistence of the metapopulation.
}

\author{TRAN Thi Cam Giang}

\seealso{
	'print.seir' function in 'dizzys' package.
}

\examples{
	obj<-seir(nbVilles=3, N=1e7)
	print(obj)
}
