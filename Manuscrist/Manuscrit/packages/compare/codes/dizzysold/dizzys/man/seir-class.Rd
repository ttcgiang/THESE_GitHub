\name{seir-class}

\Rdversion{1.1}

\docType{class}

\alias{seir-class}

\title{seir Class \code{"seir"}}

\description{This package permit us to simulate infectious diseases by using the SIR/SEIR models.To do that, we implement the direct algorithm of Gillespie in 1977 as described in the book 'Modeling Infectious Diseases IN HUMANS AND ANIMALS' by Keeling and Pejman Rohani (2008), and the adaptive tau leaping to approximate the trajectory of a continuous-time stochastic process as described by Cao et al. (2007) The Journal of Chemical Physics.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("seir", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{pop}:}{list. Object of class}
    \item{\code{duration}:}{numeric. time values over which to perform the numerical integration.}
    \item{\code{S}:}{numeric. Initial value of the state variable 'S', this is the number of susceptible individuals at the time t0.}
    \item{\code{E}:}{numeric. Initial value of the state variable 'E', this is the number of exposed individuals at the time t0.}
    \item{\code{I}:}{numeric. Initial value of the state variable 'I', this is the number of infected individuals at the time t0.}
    \item{\code{R}:}{numeric. Initial value of the state variable 'R', this is the number of recovered individuals at the time t0.}
    \item{\code{N}:}{numeric. Initial value of the state variable 'N', this is the number of population at the time t0.}
    \item{\code{T}:}{numeric. Period of the contact rate (in days).}
    \item{\code{mu}:}{numeric. Per capita birth and death rates (per day).}
    \item{\code{beta0}:}{numeric. Mean value of contact rate (per day).}
    \item{\code{beta1}:}{numeric. Amplitude of contact rate.}
    \item{\code{sigma}:}{numeric. Transition rate from exposed (E) to infected (I), per day (inverse of the average duration of the latency period).}
    \item{\code{gamma}:}{numeric. Recovery rate (per day).}
    \item{\code{unitTIME}:}{numeric. Unit of time.}
    \item{\code{phi}:}{numeric. Phase of the contact rate (in radians).}
    \item{\code{nbVilles}:}{numeric. Number of subpopulations in metapopulation.}
    \item{\code{type}:}{character. Type of simulation, 'deterministic' or 'stochastic'.}
    \item{\code{epsilon}:}{numeric. Proportion of infections by contacts with infected from other subpopulations. The value of epsilon between subpopulation i and j is equal to the value of epsilon between subpopulation j and i. epsilon is always greater than or equal to zero and less than or equal to one,specially, epsilon, of the subpopulation i and i, is equal to zero.}
    \item{\code{rho}:}{numeric. Coupling rate between subpopulations i and j. rho between subpopulation i and j is equal to rho between subpopulation j and i. rho is always greater than or equal to zero and less than or equal to one, specially, rho, of the subpopulation i and i, is equal to one.}
    \item{\code{seed}:}{numeric. Random seed for random number generator.}
    \item{\code{rng}:}{character. Random number generator, 'good' or 'fast'. If we chose rng="good", it means that we use the random number generator in C++. In the other side, we chose rng="fast", it means that we use the random number generator of the professor Yann Chevaleyre at university of Paris 13. With the random generator's Yann Chevaleyre, we find that it is faster than the random generator in C++, in contrast, its accuracy is less than that of the random generator's Yann Chevaleyre.}
    \item{\code{method}:}{character. Simulation algorithm,'direct' or 'adaptivetau'. If we chose method="direct", it means that we implement the direct algorithm of Gillespie in 1977. In contrast, if we chose method="adaptivetau", it means that we implement the adaptive tau-leaping approximation for simulating the trajectory of a continuous-time Markov process.}
    \item{\code{persistence}:}{data.frame. Data frame about the information of the persistence of the metapopulation.}
  }
}
\section{Methods}{
No methods defined with class "seir" in the signature.
}

\references{Matt J. Keeling and Pejman Rohani (2008) Modeling Infectious Diseases IN HUMANS AND ANIMALS. Princeton University Press.
	Cao Y, Gillespie DT, Petzold LR,The Journal of Chemical Physics, 2007.
	Norman Matloff (2009), The Art of R Programming.}

\author{TRAN Thi Cam Giang}

\seealso{
	'seir' class in 'dizzys' package.
}

\examples{
showClass("seir")
}
\keyword{classes}
