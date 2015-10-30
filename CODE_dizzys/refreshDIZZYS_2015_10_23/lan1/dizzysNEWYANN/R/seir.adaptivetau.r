#
# Doing stochastic simulation for SEIR and SIR model
# by using the "adaptive tau-leaping" algorithm described by Cao Y, Gillespie DT, Petzold LR, the Journal of Chemical Physics (2007). This programme interchanges C++ and R.
# Here, we call the function 'ssesAdaptiveTau' written from the C++ program.
# seir/sir sses.adaptivetau
sses.adaptivetau <- function(sigma=1/8,gamma=1/5,mu=1/(70*365),seed=as.numeric(Sys.time()),S=NULL,E=NULL,I=NULL,R=NULL,N=1e5,
		nbCONTACT0=300, nbCONTACT1=0.1, nbMulCONTACT=1, phiPHASE=c(0),
		probVISITER=0.01, probINFECTER=0.01, duration=1000, nbVilles=2, unitTIME=1,periDISE=365,
		typeRNG="good",typeSIMU="stoch",method="direct",statSTATE=FALSE,...) {

	# initial values of variables
	# initial value of variables for all populations	
	initVar<-initVarSEIRN(nbVilles=nbVilles,S=S,E=E,I=I,R=R,N=N,mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,sigma=sigma,gamma=gamma,periDISE=periDISE,phiPHASE=phiPHASE)
	initS <- initVar$S
	initE <- initVar$E
	initI <- initVar$I
	initR <- initVar$R
	initN <- initVar$N

	init.values <- c(S=initS,E=initE,I=initI,R=initR)
	print(init.values)	
	#number of variables	
	relratechange <- rep(1, length(init.values))	
	#calling the function 'ssesAdaptiveTau' from C++
	allPops <- .Call('ssesAdaptiveTauRETRY',
			#variables
      		      init.values, NULL, NULL, 
	                #parameters
		   sigma,gamma,mu,seed,nbCONTACT0,nbCONTACT1,nbMulCONTACT,
                         phiPHASE,probVISITER,probINFECTER,nbVilles,unitTIME,typeRNG,periDISE,
                         duration, FALSE, relratechange, NULL, NULL)
	print(allPops)
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
