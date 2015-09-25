#
# Doing stochastic simulation for SEIR and SIR model
# by using the "adaptive tau-leaping" algorithm described by Cao Y, Gillespie DT, Petzold LR, the Journal of Chemical Physics (2007). This programme interchanges C++ and R.
# Here, we call the function 'ssesAdaptiveTau' written from the C++ program.
# seir/sir sses.adaptivetau
#
sses.adaptivetau <- function(duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),
			beta0=1000/365,beta1=.1,sigma=1/8,gamma=1/5,
			T=365,phi=0,nbVilles=1,epsilon=0.0,rho=0.0,seed=as.numeric(Sys.time()),
			rng="good",S=NULL,E=NULL,I=NULL,R=NULL,N=NULL,...) {

	# initial values of variables
	init.values <- c(S=S,E=E,I=I,R=R)	
	#number of variables	
	relratechange <- rep(1, length(init.values))	
	#calling the function 'ssesAdaptiveTau' from C++
	allPops <- .Call('ssesAdaptiveTau',
			#variables
      		      init.values, NULL, NULL, 
	                #parameters
			nbVilles,beta0,beta1,mu,sigma,gamma,phi,rho,epsilon,T,
		        duration, FALSE, relratechange, NULL, NULL)
	# The result returned is a 2D matrix for all subpopulations
	# now, we divide the global result into 'nbVilles' small data.frame for each subpopulation
	populations <- vector("list",nbthVilles)
	names(populations)<-paste("pop",c(1:nbVilles))
	for(vil in 1:nbVilles){
		populations[[vil]] <- data.frame(allPops[,1],allPops[,0*nbVilles+vil+1],
					allPops[,1*nbVilles+vil+1],allPops[,2*nbVilles+vil+1],allPops[,3*nbVilles+vil+1])
		names(populations[[vil]]) <- c("time","S","E","I","R")
	}	
	#returning the result that is a list containing 'nbVilles' data.frame for all subpopulations	
	return(populations)	
}
#the end
