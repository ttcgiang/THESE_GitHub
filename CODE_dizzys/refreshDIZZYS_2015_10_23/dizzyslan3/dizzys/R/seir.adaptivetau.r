#seir/sir sses.adaptivetau
sses.adaptivetau <- function(duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/8,gamma=1/5,
			T=365,phi=0,nbVilles=1,epsilon=0.0,rho=0.0,seed=as.numeric(Sys.time()),rng="good",
			S=NULL,E=NULL,I=NULL,R=NULL,N=NULL,...) {

	#tmax<-duration*365
	#params		
	init.values <- c(S=S,E=E,I=I,R=R)		
	relratechange <- rep(1, length(init.values))	
	allPops <- .Call('ssesAdaptiveTau',
      		      init.values, NULL, NULL, 
		      #params,
			nbVilles,beta0,beta1,mu,sigma,gamma,phi,rho,epsilon,T,
		      #
		      duration, FALSE, relratechange, NULL, NULL)

	#matplot(allPops[,"time"], allPops[,c("I")], type='l',xlab='Time', ylab='infecteds')		
	populations <- vector("list",nbVilles)
	names(populations)<-paste("pop",c(1:nbVilles))
	for(vil in 1:nbVilles){
		populations[[vil]] <- data.frame(allPops[,1],allPops[,0*nbVilles+vil+1],
					allPops[,1*nbVilles+vil+1],allPops[,2*nbVilles+vil+1],allPops[,3*nbVilles+vil+1])
		names(populations[[vil]]) <- c("time","S","E","I","R")
	}		
	return(populations)	
}
#the end
