# "simul" remakes simulations 
# This function is a global function, it can integrate a lot of simulation types, respectively.
#
#
globSEIRSimulNewInfec<-function(object,type="deter",continue=F,
			duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),
			nbCONTACT0=1000/365,nbCONTACT1=.1,probINFECTER=0.01,probVISITER=0.1,
			sigma=1/8,gamma=1/5,
			T=365,phiPHASE=c(0),nbVilles=1,epsilon=0.0,rho=0.0,seed=as.numeric(Sys.time()),
			rng="good", append=TRUE, t0=NULL,
			S=NULL,E=NULL,I=NULL,R=NULL,N=1e7){
	
	if (missing(object)){#if no seir object, we create a new seir object by simulating
		return(globSEIRNewInfec(type=type,method=method,duration=duration,unitTIME=unitTIME,mu=mu,
			nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,
			probVISITER=probVISITER, probINFECTER=probINFECTER,
			sigma=sigma,gamma=gamma,T=T,phiPHASE=phiPHASE,nbVilles=nbVilles,epsilon=epsilon,rho=rho,
			seed=seed,rng=rng,S=S,E=E,I=I,R=R,N=N))
	}
	else{	# if there is already a seir object

		# creating a template seir object 
		objtmp<-new("seir")
		# number of population
		nbVilObj<-object@nbVilles
		#simulation time of the object given
		oldrep<- length(object@duration)

		# ckecking the  validation of parameters in the object given
		if(missing(nbVilles))	nbVilles<-object@nbVilles
		else		nbVilles<-nbVilles
		if(is.null(nbVilles)) stop("invalid value of nbVilles")
		#type of simulation
		if(missing(type))	type<-object@type
		else		type<-type
		if(is.null(type)) stop("invalid value of type")
		#duration : time of simulation	
		if(missing(duration)) duration<-sum(object@duration)
		else	duration<-duration
		objtmp@duration<-c(object@duration,duration)
		if(is.null(duration)) stop("invalid value of duration")
		#unitTIME : unit of time
		if(missing(unitTIME))	unitTIME<-object@unitTIME	
		else 		unitTIME<-unitTIME
		if(is.null(unitTIME)) stop("invalid value of unitTIME")
		#mu : birth and death rate
		if(missing(mu)){
			if(oldrep==1) mu<-object@mu
			else mu<-object@mu[((oldrep-1)*nbVilObj+1) : (oldrep*nbVilObj)]
		}	
		else mu<-mu

		if(is.null(mu)) stop("invalid value of mu")
		else{
			newmu<-c()
			for(i in 1:oldrep){
				if((i==1)||(nbVilObj==1)) oldmu<-object@mu[i:(i*nbVilObj)]
				else oldmu<-object@mu[(i+1):(i*nbVilObj)]
				newmu<-c(newmu,resizeVector(oldmu,nbVilles))
			}
			objtmp@mu<-c(newmu,resizeVector(mu,nbVilles))
		}
		#nbCONTACT0 : mean value of the number of contact per unit time \nbCONTACT0
		if(missing(nbCONTACT0)){
			if(oldrep==1) nbCONTACT0<-object@nbCONTACT0
			else nbCONTACT0<-object@nbCONTACT0[((oldrep-1)*nbVilObj+1) :(oldrep*nbVilObj)]
		}	
		else nbCONTACT0<-nbCONTACT0

		if(is.null(nbCONTACT0)) stop("invalid value of nbCONTACT0")
		else{
			newnbCONTACT0<-c()
			for(i in 1:oldrep){
				if((i==1)||(nbVilObj==1)) oldnbCONTACT0<-object@nbCONTACT0[i:(i*nbVilObj)]
				else oldnbCONTACT0<-object@nbCONTACT0[(i+1):(i*nbVilObj)]
				newnbCONTACT0<-c(newnbCONTACT0,resizeVector(oldnbCONTACT0,nbVilles))
			}
			objtmp@nbCONTACT0<-c(newnbCONTACT0,resizeVector(nbCONTACT0,nbVilles))
		}
		#nbCONTACT1 : amplitude of the number of contact per unit time \nbCONTACT1
		if(missing(nbCONTACT1)){
			if(oldrep==1) nbCONTACT1<-object@nbCONTACT1
			else nbCONTACT1<-object@nbCONTACT1[((oldrep-1)*nbVilObj+1) : (oldrep*nbVilObj)]
		}	
		else nbCONTACT1<-nbCONTACT1
		if(is.null(nbCONTACT1)) stop("invalid value of nbCONTACT1")
		else{
			newnbCONTACT1<-c()
			for(i in 1:oldrep){
				if((i==1)||(nbVilObj==1)) oldnbCONTACT1<-object@nbCONTACT1[i:(i*nbVilObj)]
				else oldnbCONTACT1<-object@nbCONTACT1[(i+1):(i*nbVilObj)]
				newnbCONTACT1<-c(newnbCONTACT1,resizeVector(oldnbCONTACT1,nbVilles))
			}
			objtmp@nbCONTACT1<-c(newnbCONTACT1,resizeVector(nbCONTACT1,nbVilles))
		}
		# probINFECTER[i,j]: is the probability that a susceptible individual native from i being in contact with another infected individual native from k gets infected. 
		if(missing(probINFECTER)){
			if(oldrep==1) probINFECTER<-object@probINFECTER
			else probINFECTER<-object@probINFECTER[((oldrep-1)*nbVilObj+1) : (oldrep*nbVilObj)]
		}	
		else probINFECTER<-probINFECTER
		if(is.null(probINFECTER)) stop("invalid value of probINFECTER")
		else{
			newprobINFECTER<-c()
			for(i in 1:oldrep){
				if((i==1)||(nbVilObj==1)) oldprobINFECTER<-object@probINFECTER[i:(i*nbVilObj)]
				else oldprobINFECTER<-object@probINFECTER[(i+1):(i*nbVilObj)]
				newprobINFECTER<-c(newprobINFECTER,resizeVector(oldprobINFECTER,nbVilles))
			}
			objtmp@probINFECTER<-c(newprobINFECTER,resizeVector(probINFECTER,nbVilles))
		}
		##
		#probVISITER[i,j] : the probability that an individual from subpopulation i visits subpopulation j
		if(missing(probVISITER)){
			if(oldrep==1) probVISITER<-object@probVISITER
			else probVISITER<-object@probVISITER[((oldrep-1)*nbVilObj+1) : (oldrep*nbVilObj)]
		}	
		else probVISITER<-probVISITER
		if(is.null(probVISITER)) stop("invalid value of probVISITER")
		else{
			newprobVISITER<-c()
			for(i in 1:oldrep){
				if((i==1)||(nbVilObj==1)) oldprobVISITER<-object@probVISITER[i:(i*nbVilObj)]
				else oldprobVISITER<-object@probVISITER[(i+1):(i*nbVilObj)]
				newprobVISITER<-c(newprobVISITER,resizeVector(oldprobVISITER,nbVilles))
			}
			objtmp@probVISITER<-c(newprobVISITER,resizeVector(probVISITER,nbVilles))
		}
		#sigma : average latent period
		if(missing(sigma)){
			if(oldrep==1) sigma<-object@sigma
			else sigma<-object@sigma[((oldrep-1)*nbVilObj+1) : (oldrep*nbVilObj)]
		}	
		else sigma<-sigma

		if(is.null(sigma)) stop("invalid value of sigma")
		else{
			newsigma<-c()
			for(i in 1:oldrep){
				if((i==1)||(nbVilObj==1)) oldsigma<-object@sigma[i:(i*nbVilObj)]
				else oldsigma<-object@sigma[(i+1):(i*nbVilObj)]
				newsigma<-c(newsigma,resizeVector(oldsigma,nbVilles))
			}
			objtmp@sigma<-c(newsigma,resizeVector(sigma,nbVilles))
		}	
		#gamma : average infectious period
		if(missing(gamma)){
			if(oldrep==1) gamma<-object@gamma
			else gamma<-object@gamma[((oldrep-1)*nbVilObj+1) : (oldrep*nbVilObj)]
		}	
		else gamma<-gamma

		if(is.null(gamma)) stop("invalid value of gamma")
		else{
			newgamma<-c()
			for(i in 1:oldrep){
				if((i==1)||(nbVilObj==1)) oldgamma<-object@gamma[i:(i*nbVilObj)]
				else oldgamma<-object@gamma[(i+1):(i*nbVilObj)]
				newgamma<-c(newgamma,resizeVector(oldgamma,nbVilles))
			}
			objtmp@gamma<-c(newgamma,resizeVector(gamma,nbVilles))
		}
		#T : period of year
		if(missing(T))	T<-object@T
		else		T<-T
		if(is.null(T)) stop("invalid value of T")
		#phi : phase of forcing
		if(missing(phiPHASE)){
			if(oldrep==1) phiPHASE<-object@phiPHASE
			else phiPHASE<-object@phiPHASE[((oldrep-1)*nbVilObj+1) : (oldrep*nbVilObj)]
		}	
		else phiPHASE<-phiPHASE

		if(is.null(phiPHASE)) stop("invalid value of phiPHASE")
		else{
			newphiPHASE<-c()
			for(i in 1:oldrep){
				if((i==1)||(nbVilObj==1)) oldphiPHASE<-object@phiPHASE[i:(i*nbVilObj)]
				else oldphiPHASE<-object@phiPHASE[(i+1):(i*nbVilObj)]
				newphiPHASE<-c(newphiPHASE,resizeVector(oldphiPHASE,nbVilles))
			}
			objtmp@phiPHASE<-c(newphiPHASE,resizeVector(phiPHASE,nbVilles))
		}
		#seed : number 'seed'
		if(missing(seed)) 	seed<-object@seed
		else 			seed<-seed
		if(is.null(seed)) stop("invalid value of seed")	
		#epsilon : infection rate from outside
		if(missing(epsilon)) 	epsilon<-object@epsilon
		else epsilon<-epsilon
		if(is.null(epsilon)) stop("invalid value of epsilon")	
		#rho : rate of coupling
		if(missing(rho)) 	rho<-object@rho
		else rho<-rho
		if(is.null(rho)) stop("invalid value of rho")	

		#rng : random number generator
		if(missing(rng)) rng<-object@rng
		else	rng = rng			
		if(length(rng) == 0L)  print("invalid value of rng, review the value rng!")

	}
	
	#CASE: there is a seir object given 
	# If continue=T: we continue the simulation of the object
	# If continue=F: we redo the simulation of the object
	if(!continue){#continue=F		
		if(missing(S)&missing(E)&missing(I)&missing(R)&missing(N)){#missing S,E,I,R,N
			S<-object@S
			E<-object@E
			I<-object@I
			R<-object@R
			N<-object@N
		}
		else{# There are S, E, I, R, N
			N<-N; 	S<-S; E<-E; I<-I; R<-R
		}
		# result
		return(globSEIRNewInfec(type=type,method=method,duration=duration,unitTIME=unitTIME,mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,
			probINFECTER=probINFECTER,probVISITER=probVISITER,
			sigma=sigma,gamma=gamma,T=T,phiPHASE=phiPHASE,nbVilles=nbVilles,epsilon=epsilon,rho=rho,
			seed=seed,rng=rng,S=S,E=E,I=I,R=R,N=N))
	}
	else{#continue=T, we will add new simulation in the tail of the simulation given at the time t0 given
		#at t0
		oldPops <- object@pop
		nbVilObj<-object@nbVilles
		tpSEIR <- data.frame()
		# checking t0
		if(is.null(t0)){# if t0 is null
			for(i in 1:nbVilObj){
			oldpopi <- oldPops[[i]]				
			tpSEIR <- rbind(tpSEIR,tail(oldpopi,1))					
			}
		}
		else if((t0>=0)&&(t0<= sum(object@duration)*365)){#if t0 is out of the simulation time of the object given
			for(i in 1:nbVilObj){
			oldpopi <- oldPops[[i]]				
			tpSEIR <- rbind(tpSEIR,tail(subset(oldpopi,oldpopi$time<=t0),1))
			}
		}
		else # error
			stop("Error, t0 should be in [0, duration] time")
		#missing S,E,I,R,N
		if(missing(S)&missing(E)&missing(I)&missing(R)&missing(N)){
			S <- tpSEIR[,2]		
			E <- tpSEIR[,3]				
			I <- tpSEIR[,4]	
			R <- tpSEIR[,5]	
			N<-S+E+I+R
		}
		else{#
			N<-N; 	S<-S; E<-E; I<-I; R<-R
		}
		
		#appending the new simulation in the tail of the simulation given at time t0
			# checking the type of simulation
				#if TYPE of the object given is 'stochastic', but TYPE of the new obj is 'deterministic'
		if(!is.integer0(grep(object@type,"stochastic")) && !is.integer0(grep(type,"deterministic"))){
			# nex object
			newObj <- globSEIRNewInfec(type=type,duration=duration,unitTIME=unitTIME,
					S=S,E=E,I=I,R=R,N=N,
					mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,
					probINFECTER=probINFECTER,probVISITER=probVISITER,
					sigma=sigma,gamma=gamma,T=T,phiPHASE=phiPHASE,	
					nbVilles=nbVilles,seed=seed)

			newPops<-newObj@pop
			for(i in 1:nbVilles){
				newpopi<-newPops[[i]]
				newpopi$time<-newpopi$time + tpSEIR[1,1]
				newPops[[i]]<-newpopi	
			}			
			if(!append){ # don't append
				newObj@pop<-newPops								
				return(newObj)
			}
			else{	# append
					
				seqvil<-resizeVector(seq(1:nbVilObj),nbVilles)
				for(i in 1:nbVilles){	
					oldpopi <- oldPops[[seqvil[i]]]
					if(!is.null(t0)) oldpopi<-subset(oldpopi,oldpopi$time<=t0)	
					oldpopi<-data.frame(time=oldpopi[,1],S=oldpopi[,2],E=oldpopi[,3],
								P=oldpopi[,4],R=oldpopi[,5],N=oldpopi[6])			
					newPops[[i]]<-rbind(oldpopi,newPops[[i]])	
				}		
				newObj@pop<-newPops
			}		
		}
		else	#if TYPE of the object given is 'deterministic', but TYPE of the new obj is 'stochastic'
		if(!is.integer0(grep(object@type,"deterministic")) && !is.integer0(grep(type,"stochastic"))){
			#nex object
			newObj <- globSEIRNewInfec(type=type,duration=duration,unitTIME=unitTIME,S=S,E=E,I=I,R=R,N=N,
				mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,
				probINFECTER=probINFECTER,probVISITER=probVISITER,
				sigma=sigma,gamma=gamma,T=T,phiPHASE=phiPHAS,
				nbVilles=nbVilles,epsilon=epsilon,rho=rho,seed=seed,rng=rng)
			newPops <- newObj@pop	
			for(i in 1:nbVilles){	
				newpopi <- newPops[[i]]
				newpopi$time <- newpopi$time + tpSEIR[1,1]
				newPops[[i]]<-newpopi
			}	

			if(!append){#don't append
				newObj@pop<-newPops
				 return(newObj)
			}
			else{#append
				#oldpop
				oldpop1<-oldPops[[1]]
				if(!is.null(t0)) oldpop1<-subset(oldpop1,oldpop1$time<=t0)				
				oldpop1<- cbind(time=oldpop1[,1],S=oldpop1[,2],E=oldpop1[,3],
										P=oldpop1[,4],R=oldpop1[,5],N=oldpop1[,6],I=oldpop1[,4])
								
				for(i in 1:nbVilles){						
					newPops[[i]] <- rbind(oldpop1,newPops[[i]])				
				}
				newObj@pop <- newPops
			}
		}
		else	#if TYPE of the object given is 'deterministic', but TYPE of the new obj is 'deterministic'
		if(!is.integer0(grep(object@type,"deterministic")) && !is.integer0(grep(type,"deterministic"))){			
			# new object
			newObj <- globSEIRNewInfec(type=type,duration=duration,unitTIME=unitTIME,
			S=S,E=E,I=I,R=R,N=N,
			mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,
			probINFECTER=probINFECTER,probVISITER=probVISITER,
			sigma=sigma,gamma=gamma,T=T,phiPHASE=phiPHASE[i],	
			nbVilles=nbVilles,seed=seed)
			newPops<-newObj@pop	
			for(i in 1:nbVilles){
				newpopi<-newPops[[i]]
				newpopi$time<-newpopi$time + tpSEIR[1,1]
				newPops[[i]]<-newpopi
			}
			if(!append){	#don't append
				newObj@pop<-newPops
				return(newObj)
			}
			else{	#append		
				seqvil<-resizeVector(seq(1:nbVilObj),nbVilles)
				for(i in 1:nbVilles){
					oldpopi <- oldPops[[seqvil[i]]]
					if(!is.null(t0)) oldpopi<-subset(oldpopi,oldpopi$time<=t0)				
					newPops[[i]]<-rbind(oldpopi,newPops[[i]])	
				}	
				newObj@pop<-newPops					
			}
		}
		else
		#if TYPE of the object given is 'stochastic', but TYPE of the new obj is 'stochastic'
		if(!is.integer0(grep(object@type,"stochastic")) && !is.integer0(grep(type,"stochastic"))){
			# new object
			newObj <- globSEIRNewInfec(type=type,duration=duration,unitTIME=unitTIME,S=S,E=E,I=I,R=R,N=N,
				mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,
				probINFECTER=probINFECTER,probVISITER=probVISITER,
				sigma=sigma,gamma=gamma,T=T,phiPHASE=phiPHASE,
				nbVilles=nbVilles,epsilon=epsilon,rho=rho,seed=seed,rng=rng)
			newPops <- newObj@pop
			for(i in 1:nbVilles){
				newpopi <- newPops[[i]]
				newpopi$time <- newpopi$time + tpSEIR[1,1]
				newPops[[i]]<-newpopi			
			}
			if(!append){#don't append
				newObj@pop<-newPops
				return(newObj)
			}
			else{#append
				
				seqvil<-resizeVector(seq(1:nbVilObj),nbVilles)
				for(i in 1:nbVilles){
					oldpopi <- oldPops[[seqvil[i]]]
					if(!is.null(t0)) oldpopi <- subset(oldpopi,oldpopi$time<=t0)
					newPops[[i]] <- rbind(oldpopi,newPops[[i]])	
				}				
						
					newObj@pop <- newPops
			}				
		}
		else	#error
			stop("Error, type should be 'stochastic' or 'deterministic'.")
		#append
		newObj@duration<-objtmp@duration
		newObj@mu<-objtmp@mu	
		newObj@nbCONTACT0<-objtmp@nbCONTACT0	
		newObj@nbCONTACT1<-objtmp@nbCONTACT1
		newObj@probINFECTER <-objtmp@probINFECTER
		newObj@probVISITER  <- newObj@probVISITER
		newObj@sigma<-objtmp@sigma
		newObj@gamma<-objtmp@gamma
		newObj@phiPHASE<-objtmp@phiPHASE	
		#result			
		return(newObj)	
	}
	
}
######################################################################################
