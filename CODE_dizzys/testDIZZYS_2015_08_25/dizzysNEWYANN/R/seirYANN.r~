#
# Defining the seir class in R
# this defines necessary parameters and variables 
#
setClass("seirYANN",
	# names of parameters and variables 
   representation(sigma="numeric", # average latent period
		gamma="numeric", # average infectious period
		mu="numeric", # birth and death rate	
		seed="numeric", # number 'seed'		
		S="numeric",E="numeric",I="numeric",R="numeric",N="numeric",	# value of variables 
  	        nbCONTACT0="numeric", 
		nbCONTACT1="numeric", 
		nbMulCONTACT="numeric",
		phiMIN="numeric", 
		phiMAX="numeric",
		probVISITER="numeric",
		probINFECTER="numeric", 
		duration="numeric",# simulation time
		nbVilles="numeric", #number of population in a metapopulation
		unitTIME="numeric", #unit of time
		periDISE="numeric", # disease period 
		typeSIMU="character", # type of simulation, 'deterministic' or 'stochastic'
		method="character",  # method of simulation, 'direct' or 'adaptivetau'
		typeRNG="character", # type of the random number generator 0: generator of C++ (good), 1: generator of Yann (fast)
		disSTAT="numeric", #stationary distribution
		lsLocalEXTPOINT="list",
		lsRecolEXTPOINT="list",
		pop="list" # containing (seir) data.frames of all populations over time				
		),

	# initial values of parameters and variables
  prototype(sigma=1/7,gamma=1/7,mu=1/(70*365),seed=23456,S=NULL,E=NULL,I=NULL,R=NULL,N=1e7,
		nbCONTACT0=300, nbCONTACT1=0.1, nbMulCONTACT=1, phiMIN=0.0, phiMAX=0.0,
		probVISITER=0.01, probINFECTER=0.01, duration=1000, nbVilles=2, unitTIME=1,periDISE=365,
		typeRNG="good",typeSIMU="stoch",method="direct",disSTAT=c(0),pop=NULL)
)

###################################################
# additional fucntion for other main fucntions
#
####
# checking a number is integer or no
#
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}
####
# sequence of time due to simulation time
#
timeseq <- function(duration=10*365, unitTIME=1){
	return(seq(0,duration,le=duration/unitTIME))
}
#####
###################################################
#
# Initializing the values of al variables for one METAPOPULATION
# according to the input or the equilibrium values after 100 years we do a deterministic simulation
#
initVarYANN <- function(nbVilles=1,S=NULL,E=NULL,I=NULL,R=NULL,N=NULL,
				mu=1/(70*365),sigma=1/8,gamma=1/5,periDISE=365,phi=0,
				nbCONTACT0=300, nbCONTACT1=0.1, probINFECTER=0.01){
	#resize vectors of variables	
	S<-resizeVector(V=S,n=nbVilles)
	E<-resizeVector(V=E,n=nbVilles)
	I<-resizeVector(V=I,n=nbVilles)
	R<-resizeVector(V=R,n=nbVilles)
	N<-resizeVector(V=N,n=nbVilles)
	#parameter
	mu<-resizeVector(V=mu,n=nbVilles)
	nbCONTACT0<-resizeVector(V=nbCONTACT0,n=nbVilles)
	nbCONTACT1<-resizeVector(V=nbCONTACT1,n=nbVilles)
	sigma<-resizeVector(V=sigma,n=nbVilles)
	gamma<-resizeVector(V=gamma,n=nbVilles)
	phi<-resizeVector(V=phi,n=nbVilles)
	
	
	vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
	nbVar <- length(vecVar)/nbVilles
	
	# checking values of variables
	if(nbVar==0) stop("There are not any initial variables!")
	else if(nbVar==5){# case 1: there are enough values of variables
		validN = TRUE
		for(i in 1:nbVilles){
			if(N[i]!=(S[i]+E[i]+I[i]+R[i])) validN <- FALSE
		}
		if(validN) return(vecVar)
		else
			stop("Values of initial variables should verify N = S + E + I + R!")

	}
	else if(nbVar==4){# case 2: there are 4 values of variables
		if(is.integer0(which(vecVar==S))){ # missing  S
			#print("missing  S")
			S<-N-E-I-R	
			vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
			return(vecVar)		
		}
		else if(is.integer0(which(vecVar==E))){# missing  E
			#print("missing  E")
			E<-N-S-I-R	
			vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
			return(vecVar)		
		}
		else if(is.integer0(which(vecVar==I))){# missing  I
			#print("missing  I")
			I<-N-S-E-R	
			vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
			return(vecVar)		
		}
		else if(is.integer0(which(vecVar==R))){# missing  R
			#print("missing  R")
			R<-N-S-I-E	
			vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
			return(vecVar)		
		}
		else if(is.integer0(which(vecVar==N))){# missing  N
			#print("missing  N")
			N<-S+E+I+R
			vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
			return(vecVar)		
		}
	}
	else if(nbVar<4){# case 3: there are less than 4 values of variables
		if((nbVar==1)&(!is.null(N))){ # there is N (population size for population)
		# here, we do a deterministic simulation in 100 years,
		# then get initial valeus for variables at the last time 
			S<-c(); E<-c(); I<-c(); R<-c()
			for(i in 1:nbVilles){				
				valEqual <- equilibriumYANN(duration=100*365,N=N[i],mu=mu[i],nbCONTACT0=nbCONTACT0[i],probINFECTER=probINFECTER,sigma=sigma[i],gamma=gamma[i],phi=phi[i],periDISE=periDISE)

				Sdet <- valEqual[1]; S = c(S,round(Sdet));
				Edet <- valEqual[2]; E = c(E,round(Edet));
				Idet <- valEqual[3]; I = c(I,round(Idet));
				R <- c(R,N[i] - round(Sdet) - round(Edet) - round(Idet)); 

			}
			vecVar<-c(S=S,E=E,I=I,R=R,N=N)	
			#print(vecVar)
			return(vecVar)			
		}
		else # there is no N, there is error
			stop("Number of initial variables should be more than 4!")
	}

}
###################
# Basic function
# function does simulation 'determionistic' or 'stochastic'
# according to the parameter 'type'
# 
seirYANN<-function(type="stoch",duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),nbCONTACT0=1000/365,nbCONTACT1=0.0,probINFECTER=0.01,sigma=1/8,gamma=1/5,
			periDISE=365,phiMIN=0,phiMAX=0,nbVilles=1,seed=as.numeric(Sys.time()),typeRNG="good",
			S=NULL,E=NULL,I=NULL,R=NULL,N=1e5,...){
	# stochastic simulation
	if(!is.integer0(grep(type,"stochastic")))
		seir.stochYANN(method=method,duration=duration,unitTIME=unitTIME,mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1, probINFECTER=probINFECTER,
			sigma=sigma,gamma=gamma,periDISE=periDISE,phiMIN=phiMIN,phiMAX=phiMAX,nbVilles=nbVilles,
			seed=seed,typeRNG=typeRNG,S=S,E=E,I=I,R=R,N=N)
	# deterministic simulation
	else if(!is.integer0(grep(type,"deterministic")))
		seir.detYANN(nbVilles=nbVilles,duration=duration,unitTIME=unitTIME,S=S,E=E,I=I,R=R,N=N,periDISE=periDISE,
			mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,probINFECTER=probINFECTER,sigma=sigma,gamma=gamma,phiMIN=phiMAX)
	# error
	else
		stop("invalid value of type!")
}
#######
#
# small function does deterministic simulation, that is called in the main function 'seir' above
# seir.det() produit un object de seir/sir en déterminist
#sigma=1/8,gamma=1/5,mu=1/(70*365),seed=as.numeric(Sys.time()),S=NULL,E=NULL,I=NULL,R=NULL,N=1e7,
#		nbCONTACT0=300, nbCONTACT1=0.1, nbMulCONTACT=1, phiMIN=0.0, phiMAX=0.0,
#		probVISITER=0.01, probINFECTER=0.01, duration=1000, nbVilles=2, unitTIME=1,periDISE=365,
seir.detYANN<-function(nbVilles=1,duration=10*365,unitTIME=1,S=NULL,E=NULL,I=NULL,R=NULL,N=1e7,periDISE=365,
			mu=1/(70*365),nbCONTACT0=300,nbCONTACT1=.1,probINFECTER=0.01,sigma=1/8,gamma=1/5,phi=0,...){
	# checking the parametres
	sigma<-resizeVector(sigma,nbVilles); #sigma

	for(i in 1:length(sigma))
		if(sigma[i]<0) stop("Error, invalid value of sigma, sigma should be more than zero.")
	gamma<-resizeVector(gamma,nbVilles); #gamma
	for(i in 1:length(gamma))
		if(gamma[i]<0) stop("Error, invalid value of gamma, gamma should be more than zero.")
	mu<-resizeVector(mu,nbVilles); #mu
	for(i in 1:length(mu))
		if(mu[i]<0) stop("Error, invalid value of mu, mu should be more than zero.")
	nbCONTACT0<-resizeVector(nbCONTACT0,nbVilles); #nbCONTACT0
	for(i in 1:length(nbCONTACT0))
		if(nbCONTACT0[i]<0) stop("Error, invalid value of nbCONTACT0, nbCONTACT0 should be more than zero.")
	nbCONTACT1<-resizeVector(nbCONTACT1,nbVilles); #nbCONTACT1
	for(i in 1:length(nbCONTACT1))
		if(nbCONTACT1[i]<0) stop("Error, invalid value of nbCONTACT1, nbCONTACT1 should be more than zero.")
	# phases	
	phi<-resizeVector(phi,nbVilles); #phi
	# initial values for all variables
	initVar<-initVarYANN(nbVilles=nbVilles,S=S,E=E,I=I,R=R,N=N,mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,probINFECTER=probINFECTER,
				sigma=sigma,gamma=gamma,periDISE=periDISE,phi=phi)
	#S
	S<-initVar[1:(1*nbVilles)]
	if(!is.integer(as.integer(S))||(S<0)) stop("Error, invalid value of S, S should be an positive integer.")
	#E
	E<-initVar[(1*nbVilles+1):(2*nbVilles)]
	if(!is.integer(as.integer(E))||(E<0)) stop("Error, invalid value of E, E should be an positive integer.")
	#I
	I<-initVar[(2*nbVilles+1):(3*nbVilles)]
	if(!is.integer(as.integer(I))||(I<0)) stop("Error, invalid value of I, I should be an positive integer.")
	#R
	R<-initVar[(3*nbVilles+1):(4*nbVilles)]
	if(!is.integer(as.integer(R))||(R<0)) stop("Error, invalid value of R, R should be an positive integer.")
	#N
	N<-initVar[(4*nbVilles+1):(5*nbVilles)]
	if(!is.integer(as.integer(N))||(N<0)) stop("Error, invalid value of N, N should be an positive integer.")

	init.varib <- c(S=S,E=E,I=I,R=R)

	times <- timeseq(duration,unitTIME)
	
	# Solving the system of differential equations 
	# (with deSolve package and code in integration between C++ and R)
	populations<-vector("list",nbVilles)
	for(i in 1:nbVilles){
		#parameters
		pars<-c(nbCONTACT0=nbCONTACT0[i],nbCONTACT1=nbCONTACT1[i],probINFECTER=probINFECTER,periDISE=periDISE,phi=phi[i],mu=mu[i],sigma=sigma[i],gamma=gamma[i],N=N[i])

#		print(pars)
		#variables
		yini<-c(t=times[1],S=S[i],E=E[i],I=I[i])
		# calling the function 'ode' in the package 'deSolve' 
		require("deSolve")
		# and at the same time, the function 'derivs' in the code file C++
		pops.det<-as.data.frame(ode(func = "derivs", y = yini, parms = pars, times = times, dllname = "dizzysNEWYANN", initfunc = "initmod"))
		valueN<-rep(N[i],nrow(pops.det))
		# Gathering the result
		valueR<-valueN-pops.det[,3]-pops.det[,4]-pops.det[,5]
		pops.det<-data.frame(pops.det[,1], pops.det[,3],pops.det[,4],pops.det[,5],valueR,valueN)
		# adding the names of variables into the result
		names(pops.det)<-c("time","S","E","P","R","N")		
		names(populations)<-paste("pop",i)
		populations[[i]]<-pops.det
	}
	#return a seir object
	return(new("seirYANN",nbVilles=nbVilles,duration=duration,unitTIME=unitTIME,
				N=N,S=S,E=E,I=I,R=R,periDISE=periDISE,mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,probINFECTER=probINFECTER,
					sigma=sigma,gamma=gamma,phiMIN=phi,phiMAX=phi,pop=populations,typeSIMU="deterministic"))
}
##############
#R on the server can't work with C++11 (this is a new library in C++ from 2011)
#This function uses the uniform distribution of R to generate random numbers. (trainsition matrix)
#After, we find the stationary distribution
#

#transition matrix
transMATRIX<-function(nbVilles,lower, upper)
{
	lower<-0.0
	upper<-1/nbVilles
	for (i in 1:nbVilles) {
		vecrow<-seq(nbVilles,0.0)
		n<-(nbVilles-1)
		vecGENEtp<-runif(n,min=lower,max=upper)
		valVIL=1.0-sum(vecGENEtp)
		if(i==1) vecGENE<-c(valVIL,vecGENEtp)
		else 
		if(i<nbVilles)
			vecGENE<-c(vecGENEtp[1:(i-1)],valVIL,vecGENEtp[i:n])
		else
		if(i==nbVilles)
			vecGENE<-c(vecGENEtp,valVIL)
	
		if(i==1){
#			vecGENE<-as.vector(vecGENE)
			resMATRIX<-as.data.frame(t(vecGENE))
		}
		else{
			resMATRIX<-rbind(resMATRIX,vecGENE)
		}
    }

    return (as.matrix(resMATRIX));
}
# stationary distribution
# multifly 100times
statDIS<-function(transMATRIX=matrix()){
	x<-transMATRIX
	y<-x
	for(i in 1:100)
		y<-y%*%x
	return(y)
}
###########
#
# function does one stochastic simulation,
# result is a stochastic seir object
# parameter: 'typeRNG' has two names, "good" or "fast", 'good' is the random number generator of C++, 'fast' is this of Yann
# parameter: 'method' also has two names, "direct" or "adaptivetau"
# all parameters left are vectors,
# sigma may be 'infinity' if we want SIR simulation
seir.stochYANN<-function(sigma=1/8,gamma=1/5,mu=1/(70*365),seed=as.numeric(Sys.time()),S=NULL,E=NULL,I=NULL,R=NULL,N=1e5,
		nbCONTACT0=300, nbCONTACT1=0.1, nbMulCONTACT=1, phiMIN=0.0, phiMAX=0.0,
		probVISITER=0.01, probINFECTER=0.01, duration=1000, nbVilles=2, unitTIME=1,periDISE=365,
		typeRNG="good",typeSIMU="stoch",method="direct",disSTAT=NULL,...){
	# initial value of variables for all populations	
	initVar<-initVarYANN(nbVilles=nbVilles,S=S,E=E,I=I,R=R,N=N,mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,sigma=sigma,gamma=gamma,periDISE=periDISE,phi=phiMIN)
	#S
	S<-initVar[1:(1*nbVilles)]
	if(!is.integer(as.integer(S))||(S<0)) stop("Error, invalid value of S, S should be an positive integer.")
	#E
	E<-initVar[(1*nbVilles+1):(2*nbVilles)]
	if(!is.integer(as.integer(E))||(E<0)) stop("Error, invalid value of E, E should be an positive integer.")
	#I
	I<-initVar[(2*nbVilles+1):(3*nbVilles)]
	if(!is.integer(as.integer(I))||(I<0)) stop("Error, invalid value of I, I should be an positive integer.")
	#R
	R<-initVar[(3*nbVilles+1):(4*nbVilles)]
	if(!is.integer(as.integer(R))||(R<0)) stop("Error, invalid value of R, R should be an positive integer.")
	#N
	N<-initVar[(4*nbVilles+1):(5*nbVilles)]
	if(!is.integer(as.integer(N))||(N<0)) stop("Error, invalid value of N, N should be an positive integer.")

	init.varib <- c(S=S,E=E,I=I,R=R)
	#print(init.varib)
	S0<-S[1]
	E0<-E[1]
	I0<-I[1]
	R0<-R[1]
	N0<-S0+E0+I0+R0

	# stochastic simulation of n populations
	# for the DIRECT algorithm
	if(!is.integer0(grep(method,"direct"))){

		if(!is.integer0(grep(typeRNG,"good"))) rtypeRNG=0
		else
		  if(!is.integer0(grep(typeRNG,"fast"))) rtypeRNG=1
		else
			stop("invalid value of typeRNG, it should be 'good' or 'fast'")

		#finding the stationary distribution
		if (is.null(disSTAT)){
			matrix01<-transMATRIX(nbVilles,0,1)
			statMatriw<-statDIS(matrix01)
			statline<-statMatriw[1,]			
		}
		else{
			statline=resizeVector(disSTAT,nbVilles)			
		}
			
		#gathering the result
		# here, we call the function 'getValeurPOPS' from C++
		resMetaPOP <- .Call("getValeurPOPSYANN",sigma,gamma,mu,seed,S0,E0, I0,R0,
				nbCONTACT0,nbCONTACT1,nbMulCONTACT,phiMIN,phiMAX,probVISITER,probINFECTER,
				duration,nbVilles,unitTIME, rtypeRNG, periDISE,statline)

		populations <- vector("list",nbVilles)
		names(populations)<-paste("pop",c(1:nbVilles))
		nbElement<-nbVilles+2

		for(i in 1:nbElement){
			if(i==(nbVilles+1)) {				
				resLocalExt<- resMetaPOP[[i]]				
			}
			else if(i==(nbVilles+2)) {
				resRecolExt<-resMetaPOP[[i]]
			}
			else {
				populations[[i]] <- data.frame(resMetaPOP[[i]])
				names(populations[[i]]) <- c("time","S","E","P","R","N","I")
			}
		}	
	}
	else
	stop("Object method should only be 'direct'")
	#lsLocalEXTPOINT="list",
	#	lsRecolEXTPOINT="list",
	#	lsLocalEXTDIS="list",
	#	lsRecolEXTDIS="list",
	#calculating the difference between the moment of the beginning of the local extinction and of the end of the local extinction$
	lsLocalEXTPOINT <- vector("list",nbVilles)
	lsRecolEXTPOINT <- vector("list",nbVilles)
	for( ivil in 1:nbVilles){
		nbRecolEXT<-length(resRecolExt[[ivil]])
		if(nbRecolEXT==0){
			lsLocalEXTPOINT[[ivil]]<-NULL
			lsRecolEXTPOINT[[ivil]]<-NULL
		}
		else{
			lsLocalEXTPOINT[[ivil]]<-resLocalExt[[ivil]]
			lsRecolEXTPOINT[[ivil]]<-resRecolExt[[ivil]]
			#			
		}
		
		
	}
	
	#return an object_stoch
	return(new("seirYANN",sigma=sigma,gamma=gamma,mu=mu,seed=seed,S=S0,E=E0,I=I0,R=R0,N=N0,
		nbCONTACT0=nbCONTACT0, nbCONTACT1=nbCONTACT1, nbMulCONTACT=nbMulCONTACT, phiMIN=phiMIN, phiMAX=phiMAX,
		probVISITER=probVISITER, probINFECTER=probINFECTER, duration=duration, nbVilles=nbVilles, unitTIME=unitTIME,periDISE=periDISE,
		typeRNG=typeRNG,typeSIMU="stochastic",method="direct",disSTAT=c(as.vector(statline)),lsLocalEXTPOINT=lsLocalEXTPOINT,lsRecolEXTPOINT=lsRecolEXTPOINT,
		pop=populations))
#the end
}

###############


