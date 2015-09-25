#
# Defining the seir class in R
# this defines necessary parameters and variables 
#
setClass("seir",
	# names of parameters and variables 
   representation(pop="list", # containing (seir) data.frames of all populations over time
		# variables S: susceptible, E: exposed, I: infected, R: recovered, N: population size
		S="numeric",E="numeric", I="numeric",R="numeric",N="numeric",		
		T="numeric", # disease period 
		duration="numeric", # simulation time
		mu="numeric", # birth and death rate
		beta0="numeric",beta1="numeric", # mean value and amplitude of contact rate \beta
		sigma="numeric", # average latent period
		gamma="numeric", # average infectious period
		unitTIME="numeric", #unit of time
		phi="numeric", # phase of forcing
		nbVilles="numeric", #number of population in a metapopulation
		type="character", # type of simulation, 'deterministic' or 'stochastic'
		epsilon="numeric",# infection rate from outside
		rho="numeric", # rate of coupling
		seed="numeric", # number 'seed'
		typeRNG="character", # type of the random number generator 0: generator of C++, 1: generator of Yann
		method="character", # algorithm of simulation, 'direct' or 'adaptivetau'
		persistence="data.frame" # global persistence in a metapopulation
		),

	# ial values of parameters and variables
  prototype(duration=1,S=NULL,E=NULL,I=NULL,R=NULL,N=1e7,T=365,mu=1/(70*365),
		beta0=1000/365,beta1=.1,sigma=1/7,gamma=1/7,pop=NULL,
		unitTIME=1,phi=0,nbVilles=1,epsilon=0.0,rho=0.0,
		seed=as.numeric(Sys.time()),typeRNG="good",method="direct")
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
# sequence of the phase of forcing for each population
# due to the interval [phiMIN,phiMAX]
#
phiseq <- function(phiMIN=0, phiMAX=0, nbVilles=1)
	 return(seq(from=phiMIN,to=phiMAX,length.out=nbVilles))
#####
# function 'pause'
#
pause <- function() readline("press Enter to continue!")
#####
# resize a vector with length given
# exemple:
# V=c(1,10),n=6, so Output:  1  1  1 10 10 10
# V=c(1,10),n=5, so Output:  1 10  1 10  1
# V=c(1,10),n=3, so Output:  1 10  1 
resizeVector <- function(V=c(1,10),n=6)
{
	if(is.null(V)) return(V)
	else{
		legV <- length(V)
		size <- ceiling(n/legV)
		if(n > legV){
			res<-c()
			if((n %% 2)==0){			
				for(i in 1:legV) res<-cbind(res,rep(V[i],size))
				res <- res[1:n]		
			}
			else{
				res<-rep(V,n)[1:n]		
			}
			return(res)
		}
		else
			if(n<legV) return(V[1:n])
		else return(V)
	}
}
###################################################
#
# Initializing the values of al variables for all populations
# according to the input or the equilibrium values after 100 years we do a deterministic simulation
#
initVariables <- function(nbVilles=1,S=NULL,E=NULL,I=NULL,R=NULL,N=NULL,
				mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/8,gamma=1/5,T=365,phi=0){
	#resize vectors of variables	
	S<-resizeVector(V=S,n=nbVilles)
	E<-resizeVector(V=E,n=nbVilles)
	I<-resizeVector(V=I,n=nbVilles)
	R<-resizeVector(V=R,n=nbVilles)
	N<-resizeVector(V=N,n=nbVilles)
	#parameter
	mu<-resizeVector(V=mu,n=nbVilles)
	beta0<-resizeVector(V=beta0,n=nbVilles)
	beta1<-resizeVector(V=beta1,n=nbVilles)
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
				valEqual <- equilibrium(duration=100*365,N=N[i],mu=mu[i],beta0=beta0[i],beta1=beta1[i],sigma=sigma[i],gamma=gamma[i],phi=phi[i],T=T)
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
seir<-function(type="stoch",duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),beta0=1000/365,beta1=0.0,sigma=1/8,gamma=1/5,
			T=365,phi=0,nbVilles=1,epsilon=0.0,rho=0.0,seed=as.numeric(Sys.time()),typeRNG="good",
			S=NULL,E=NULL,I=NULL,R=NULL,N=1e5,...){
	# stochastic simulation
	if(!is.integer0(grep(type,"stochastic")))
		seir.stoch(method=method,duration=duration,unitTIME=unitTIME,mu=mu,beta0=beta0,beta1=beta1,
			sigma=sigma,gamma=gamma,T=T,phi=phi,nbVilles=nbVilles,epsilon=epsilon,rho=rho,
			seed=seed,typeRNG=typeRNG,S=S,E=E,I=I,R=R,N=N)
	# deterministic simulation
	else if(!is.integer0(grep(type,"deterministic")))
		seir.det(nbVilles=nbVilles,duration=duration,unitTIME=unitTIME,S=S,E=E,I=I,R=R,N=N,T=T,
			mu=mu,beta0=beta0,beta1=beta1,sigma=sigma,gamma=gamma,phi=phi)
	# error
	else
		stop("invalid value of type!")
}
#######
#
# small function does deterministic simulation, that is called in the main function 'seir' above
# seir.det() produit un object de seir/sir en déterminist
seir.det<-function(nbVilles=1,duration=10*365,unitTIME=1,S=NULL,E=NULL,I=NULL,R=NULL,N=1e7,T=365,
			mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/8,gamma=1/5,phi=0,...){
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
	beta0<-resizeVector(beta0,nbVilles); #beta0
	for(i in 1:length(beta0))
		if(beta0[i]<0) stop("Error, invalid value of beta0, beta0 should be more than zero.")
	beta1<-resizeVector(beta1,nbVilles); #beta1
	for(i in 1:length(beta1))
		if(beta1[i]<0) stop("Error, invalid value of beta1, beta1 should be more than zero.")
	phi<-resizeVector(phi,nbVilles); #phi

	# initial values for all variables
	initVar<-initVariables(nbVilles=nbVilles,S=S,E=E,I=I,R=R,N=N,mu=mu,beta0=beta0,beta1=beta1,sigma=sigma,gamma=gamma,T=T,phi=phi)
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


	times <- timeseq(duration,unitTIME)
	
	# Solving the system of differential equations 
	# (with deSolve package and code in integration between C++ and R)
	populations<-vector("list",nbVilles)
	for(i in 1:nbVilles){
		#parameters
		pars<-c(beta0=beta0[i],beta1=beta1[i],T=T,phi=phi[i],mu=mu[i],sigma=sigma[i],gamma=gamma[i],N=N[i])
#		print(pars)
		#variables
		yini<-c(t=times[1],S=S[i],E=E[i],I=I[i])
		# calling the function 'ode' in the package 'deSolve' 
		# and at the same time, the function 'derivs' in the code file C++
		pops.det<-as.data.frame(ode(func = "derivsGENE", y = yini, parms = pars, times = times, dllname = "dizzysNEWYANN", initfunc = "initparmsGENE"))
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
	return(new("seir",nbVilles=nbVilles,duration=duration,unitTIME=unitTIME,
				N=N,S=S,E=E,I=I,R=R,T=T,mu=mu,beta0=beta0,beta1=beta1,
					sigma=sigma,gamma=gamma,phi=phi,pop=populations,type="deterministic"))
}
###########
#
# function does one stochastic simulation,
# result is a stochastic seir object
# parameter: 'typeRNG' has two names, "good" or "fast", 'good' is the random number generator of C++, 'fast' is this of Yann
# parameter: 'method' also has two names, "direct" or "adaptivetau"
# all parameters left are vectors,
# sigma may be 'infinity' if we want SIR simulation
seir.stoch<-function(method="direct",nbVilles=1,sigma=1/8,gamma=1/5,mu=1/(70*365),beta0=1000/365,beta1=0.0,
			phi=0,rho=0.0,seed=as.numeric(Sys.time()),unitTIME=1,
			duration=5*365,typeRNG="good",T=365,S=NULL,E=NULL,I=NULL,R=NULL,N=1e5,...){
	# initial value of variables for all populations	
	initVar<-initVariables(nbVilles=nbVilles,S=S,E=E,I=I,R=R,N=N,mu=mu,beta0=beta0,beta1=beta1,sigma=sigma,gamma=gamma,T=T,phi=phi)
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
	print(init.varib)

	#phi
	phi=resizeVector(phi,nbVilles)

	# stochastic simulation of n populations
	# for the DIRECT algorithm
	if(!is.integer0(grep(method,"direct"))){

		if(!is.integer0(grep(typeRNG,"good"))) rtypeRNG=0
		else
		  if(!is.integer0(grep(typeRNG,"fast"))) rtypeRNG=1
		else
			stop("invalid value of typeRNG, it should be 'good' or 'fast'")

		#gathering the result
		# here, we call the function 'getValeurPOPS' from C++
		pops <- .Call("getValeurPOPSGENE",nbVilles,sigma,gamma,mu,beta0,beta1,phi,rho,
						seed,unitTIME,duration,rtypeRNG,T,init.varib)
		populations <- vector("list",nbVilles)
		names(populations)<-paste("pop",c(1:nbVilles))
		for(i in 1:nbVilles){
			populations[[i]] <- data.frame(pops[i])
			names(populations[[i]]) <- c("time","S","E","P","R","N","I")
		}	
	}
	else
	# for the ADAPTIVETAU algorithm
	if(!is.integer0(grep(method,"adaptivetau"))){	
		flagSEIR = TRUE;
		if(is.infinite(sigma))	flagSEIR=FALSE;		
			
		if(flagSEIR){# stochastic simulation for SEIR model
			init.values <- c(S=S,E=E,I=I,R=R)	
			relratechange <- rep(1,length(init.values))	
			#calling the function 'ssesAdaptiveTau' from C++
			allPops <- .Call('ssesAdaptiveTauYANN',init.values,NULL,NULL, 
			      #params,
			      nbVilles,beta0,beta1,mu,sigma,gamma,phi,rho,T,
			      #
			      tmax,FALSE,relratechange,NULL,NULL)
			populations <- vector("list",nbVilles)
			names(populations)<-paste("pop",c(1:nbVilles))
			for(vil in 1:nbVilles){
				time<-allPops[,1]
				S<-allPops[,0*nbVilles+vil+1]
				E<-allPops[,1*nbVilles+vil+1]
				I<-allPops[,2*nbVilles+vil+1]
				R<-allPops[,3*nbVilles+vil+1]
				N<-S+E+I+R
			populations[[vil]] <- data.frame(time,S,E,I,R,N,I)
			names(populations[[vil]]) <- c("time","S","E","P","R","N","I")
			}	
		}
		else{	# stochastic simulation for SIR model
			init.values <- c(S=S,I=I,R=R)
			# sigma == 'infinity'	
			#sigma<-resizeVector(Inf,nbVilles);
			relratechange <- rep(1,length(init.values))	
			allPops <- .Call('ssesAdaptiveTau',init.values,NULL,NULL, 
			      #params,
			      nbVilles,beta0,beta1,mu,sigma,gamma,phi,rho,T,
			      #
			      tmax,FALSE,relratechange,NULL,NULL)
			populations <- vector("list",nbVilles)
			names(populations)<-paste("pop",c(1:nbVilles))
			for(vil in 1:nbVilles){
				time<-allPops[,1]
				S<-allPops[,0*nbVilles+vil+1]
				E<-rep(E[vil],length(allPops[,1]))
				I<-allPops[,1*nbVilles+vil+1]
				R<-allPops[,2*nbVilles+vil+1]
				N<-S+E+I+R
				populations[[vil]] <- data.frame(time,S,E,I,R,N,I)
				names(populations[[vil]]) <-c("time","S","E","P","R","N","I")
			}	
		}
	}
	else
	stop("Object method should be 'direct' or 'adaptivetau'!")
	
	#return an object_stoch
	return(new("seir",duration=duration,unitTIME=unitTIME,T=T,S=S,E=E,I=I,R=R,N=N,
				mu=resizeVector(mu,nbVilles),beta0=resizeVector(beta0,nbVilles),beta1=resizeVector(beta1,nbVilles),
				sigma=resizeVector(sigma,nbVilles),gamma=resizeVector(gamma,nbVilles),
				phi=phi,nbVilles=nbVilles,
				rho=rho,seed=seed,
				typeRNG=typeRNG,method=method,pop=populations,type="stochastic"))
#the end
}
################################################
#
# Here, we are interested in the global disease persistence in a metapopulation
# We will arrange persistence time of each population in an increase.
# We save the name of population and its disease persistence time
# 
# Input: a seir object
#
setGeneric("persistence", function(object,...) standardGeneric("persistence"))
setMethod("persistence", "seir",
function(object,...){
	# number of population
	nbVilles <- object@nbVilles
	# population in the metapopulation
	populations <- object@pop
	#time of simulation
	dura<-object@duration
	# calculating the persistence time of each population
	persis <- sort(sapply(populations,function(x) tail(subset(x,x$I>0),1)[1,1]))
	nbresVil<-nbVilles
	vecresVil<-c()
	nbVildie<-c()
	reVecTime<-c()
	for(i in 1:nbVilles)
	{
		 if(persis[i]<dura){
			reVecTime<-c(reVecTime,persis[i])
			vecresVil<-c(vecresVil,nbresVil)
			nbVildie<-c(nbVildie,1)
			nbresVil<-nbresVil-1
		}
	}
	if(is.vector(reVecTime)){
		pers <- cbind(time=reVecTime)
		allVilles <- dimnames(pers)[[1]]
		allVilles <- as.numeric(lapply(strsplit(allVilles," "), function(x) x[2]))
		pers <- cbind(time=reVecTime,resVil=vecresVil,ndie=nbVildie,ville=allVilles)
		object@persistence <-as.data.frame(pers)
		
	}
	else
		object@persistence <-data.frame(time=c(),resVil=c(),ndie=c(),ville=c())
	# returning the result
	return(object)
}
)
##################
#
# plotting trajectory of global persistence time in a metapopulation
# There are two types we can plot
# type I: Kaplan–Meier curve of the disease persistence time
# type II: finding the entier extinction position of each population 
#
setGeneric("plot.pers", function(object,...) standardGeneric("plot.pers"))
setMethod("plot.pers", "seir",
function(object,x="time",y="P",type="s",col="red",xlim=c(),ylim=c(),curvetype="KM",vilabline=c(),add=F,
		xlab="time",ylab="#populations non extinct",unitTIME=1,...){
	# gathering the value of persistence in a metapopulation
	pers <- object@persistence
	#time according to unit of time
	pers$time<-pers$time/unitTIME
	leng=nrow(pers)
	if(leng>0){# if there is extinction
		#number of population
		nbVilles <- object@nbVilles
		#time of simulation
		duration <- object@duration
		# TYPE I: Kaplan–Meier curve
		if(!is.integer0(grep(curvetype,"KM"))){
			if(add==F) # creating a new plot
			plot(pers$time,pers$resVil,type=type,col=col,xlab=xlab,ylab=ylab,...)
			else	# adding a new Kaplan–Meier curve on a plot given
				lines(pers$time,pers$resVil,type=type,col=col,...)
		}
		else # TYPE I: population curve
		if(!is.integer0(grep(curvetype,"population"))){
			# plotting all fluctuations of the number of infected for all populations
			plot.seir(object,y=y,type=type,col=col,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
			pop<-c(1:nbVilles)
			allVillesEX <- pers$ville
			nbVilExt<-length(allVillesEX)
			nbvilabline <- length(vilabline)	
			if(is.null(pop)){# 'pop', population that we want plot
				for(i in 1:nbvilabline)
					for(j in 1:nbVilExt) 
						if(vilabline[i]==allVillesEX[j])
							abline(v=pers$time[j])	
			}
			else{
				for(i in 1:nbvilabline){# finding position of extinction
					flagEX_i <- FALSE
					if(!is.integer0(which(pop==vilabline[i]))){
						for(j in 1:nbVilExt)
							if(vilabline[i]==allVillesEX[j]){
								abline(v=pers$time[j])	
								flagEX_i=TRUE
							}							
					}
					if(!flagEX_i) print(paste("city ",vilabline[i],"has no local extinction"))
				}
			}
		}
		else
			stop("Error, we have two curve types, curvetype should be 'KM' or 'population' !")
	}
	else  #if there is no extinction
		stop("Error, there is no extinction!")
}
)
###################################################################################
# Estimating the survival parameters.
# calculing Survival Probability
# using the package 'survival'
#
####
# Kaplan-Meier estimate with 95% confidence bounds
setGeneric("surv.prob", function(object,...) standardGeneric("surv.prob"))
setMethod("surv.prob", "seir",
   function(object,...){	
	pers<-object@persistence
 	my.fit <- survfit(Surv(pers[,1],pers[,3])~1)
	return(my.fit)
  }
)
# plotting  estimated Kaplan-Meier curve with 95% confidence bounds
setGeneric("plot.surv.prob", function(object,...) standardGeneric("plot.surv.prob"))
setMethod("plot.surv.prob", "seir",
   function(object,...){
	my.fit <- surv.prob(object)
	survival:::plot.survfit(my.fit, main="Kaplan-Meier estimate with 95% confidence bounds",xlab="time", ylab="survival probability")
}
)
# estimating the global persistence rate
# by basing on Parametric Regression Models
# 
# INPUT: a seir object that has persistence
# OUTPUT: estimated persistence rate
pers.rate.obj<-function(object){
	perobj<-persistence(object)@persistence
	leng<-nrow(perobj)
	# time and ndie
	time <- perobj$time	
	ndie<- perobj$ndie
	para.obj<-survreg(Surv(time,ndie)~1, dist="exp")
	return(para.obj)
}
#####
# Confidence interval of estimated persistence rate
# INPUT:
#	pers.rate.obj : is a object of the estimated persistence rate
#	level: level of confidence

# OUTPUT:
#	confidence interval of the estimated persistence rate		
confint.pers.rate<-function(pers.rate.obj, level=0.95){
	if(is.null(pers.rate.obj)){
		print(pers.rate.obj)
		stop("See the persistence, no local extinction in the metapopulation!")		
	}
	else{
		return(confint(pers.rate.obj,level=level))

	}
}

#############################################################################################

