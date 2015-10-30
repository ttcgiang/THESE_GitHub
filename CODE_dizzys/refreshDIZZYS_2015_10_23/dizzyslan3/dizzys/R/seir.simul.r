# "simul" remakes simulations 
#revalidparm<-function(changed=0,oldparm=0.1,newparm=0.2){
#	return(c(changed=1,oldparm=oldparm,newparm=newparm))
#}
simul<-function(object,type="deter",continue=F,
			duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/8,gamma=1/5,
			T=365,phi=0,nbVilles=1,epsilon=0.0,rho=0.0,seed=as.numeric(Sys.time()),rng="good",
			append=TRUE, t0=NULL,
			S=NULL,E=NULL,I=NULL,R=NULL,N=1e7){
	#step 1: get all values of parameters
	if (missing(object)){
		return(seir(type=type,method=method,duration=duration,unitTIME=unitTIME,mu=mu,beta0=beta0,beta1=beta1,
			sigma=sigma,gamma=gamma,T=T,phi=phi,nbVilles=nbVilles,epsilon=epsilon,rho=rho,
			seed=seed,rng=rng,S=S,E=E,I=I,R=R,N=N))
	}
	else{
		objtmp<-new("seir")
		#nbVilles
		nbVilObj<-object@nbVilles
		oldrep<- length(object@duration)

		if(missing(nbVilles))	nbVilles<-object@nbVilles
		else		nbVilles<-nbVilles
		if(is.null(nbVilles)) stop("invalid value of nbVilles")
		#type
		if(missing(type))	type<-object@type
		else		type<-type
		if(is.null(type)) stop("invalid value of type")
		#duration		
		if(missing(duration)) duration<-sum(object@duration)
		else	duration<-duration
		objtmp@duration<-c(object@duration,duration)
		if(is.null(duration)) stop("invalid value of duration")
		#unitTIME
		if(missing(unitTIME))	unitTIME<-object@unitTIME	
		else 		unitTIME<-unitTIME
		if(is.null(unitTIME)) stop("invalid value of unitTIME")
		#mu
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
		#beta0
		if(missing(beta0)){
			if(oldrep==1) beta0<-object@beta0
			else beta0<-object@beta0[((oldrep-1)*nbVilObj+1) :(oldrep*nbVilObj)]
		}	
		else beta0<-beta0

		if(is.null(beta0)) stop("invalid value of beta0")
		else{
			newbeta0<-c()
			for(i in 1:oldrep){
				if((i==1)||(nbVilObj==1)) oldbeta0<-object@beta0[i:(i*nbVilObj)]
				else oldbeta0<-object@beta0[(i+1):(i*nbVilObj)]
				newbeta0<-c(newbeta0,resizeVector(oldbeta0,nbVilles))
			}
			objtmp@beta0<-c(newbeta0,resizeVector(beta0,nbVilles))
		}
		#beta1
		if(missing(beta1)){
			if(oldrep==1) beta1<-object@beta1
			else beta1<-object@beta1[((oldrep-1)*nbVilObj+1) : (oldrep*nbVilObj)]
		}	
		else beta1<-beta1
		if(is.null(beta1)) stop("invalid value of beta1")
		else{
			newbeta1<-c()
			for(i in 1:oldrep){
				if((i==1)||(nbVilObj==1)) oldbeta1<-object@beta1[i:(i*nbVilObj)]
				else oldbeta1<-object@beta1[(i+1):(i*nbVilObj)]
				newbeta1<-c(newbeta1,resizeVector(oldbeta1,nbVilles))
			}
			objtmp@beta1<-c(newbeta1,resizeVector(beta1,nbVilles))
		}
		#sigma
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
		#gamma
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
		#T
		if(missing(T))	T<-object@T
		else		T<-T
		if(is.null(T)) stop("invalid value of T")
		#phi
		if(missing(phi)){
			if(oldrep==1) phi<-object@phi
			else phi<-object@phi[((oldrep-1)*nbVilObj+1) : (oldrep*nbVilObj)]
		}	
		else phi<-phi

		if(is.null(phi)) stop("invalid value of phi")
		else{
			newphi<-c()
			for(i in 1:oldrep){
				if((i==1)||(nbVilObj==1)) oldphi<-object@phi[i:(i*nbVilObj)]
				else oldphi<-object@phi[(i+1):(i*nbVilObj)]
				newphi<-c(newphi,resizeVector(oldphi,nbVilles))
			}
			objtmp@phi<-c(newphi,resizeVector(phi,nbVilles))
		}
		#seed
		if(missing(seed)) 	seed<-object@seed
		else 			seed<-seed
		if(is.null(seed)) stop("invalid value of seed")	
		#epsilon
		if(missing(epsilon)) 	epsilon<-object@epsilon
		else epsilon<-epsilon
		if(is.null(epsilon)) stop("invalid value of epsilon")	
		#rho
		if(missing(rho)) 	rho<-object@rho
		else rho<-rho
		if(is.null(rho)) stop("invalid value of rho")	

		#rng
		if(missing(rng)) rng<-object@rng
		else	rng = rng			
		if(length(rng) == 0L)  print("invalid value of rng, review the value rng!")

	}
	
	#Si l'object fourni
	#Si continue=T: on continue la simulation de l'object
	#Si continue=F: on refait la simulation de l'object
	if(!continue){#continue=F
		#missing S,E,I,R,N
		if(missing(S)&missing(E)&missing(I)&missing(R)&missing(N)){
			S<-object@S
			E<-object@E
			I<-object@I
			R<-object@R
			N<-object@N
		}
		else{
			N<-N; 	S<-S; E<-E; I<-I; R<-R
		}
		return(seir(type=type,method=method,duration=duration,unitTIME=unitTIME,mu=mu,beta0=beta0,beta1=beta1,
			sigma=sigma,gamma=gamma,T=T,phi=phi,nbVilles=nbVilles,epsilon=epsilon,rho=rho,
			seed=seed,rng=rng,S=S,E=E,I=I,R=R,N=N))
	}
	else{#continue=T
		#at t0
		oldPops <- object@pop
		nbVilObj<-object@nbVilles
		tpSEIR <- data.frame()

		if(is.null(t0)){
			for(i in 1:nbVilObj){
			oldpopi <- oldPops[[i]]				
			tpSEIR <- rbind(tpSEIR,tail(oldpopi,1))					
			}
		}
		else if((t0>=0)&&(t0<= sum(object@duration)*365)){
			for(i in 1:nbVilObj){
			oldpopi <- oldPops[[i]]				
			tpSEIR <- rbind(tpSEIR,tail(subset(oldpopi,oldpopi$time<=t0),1))
			}
		}
		else
			stop("Error, t0 should be in [0, duration] time")
		#missing S,E,I,R,N
		if(missing(S)&missing(E)&missing(I)&missing(R)&missing(N)){
			S <- tpSEIR[,2]		
			E <- tpSEIR[,3]				
			I <- tpSEIR[,4]	
			R <- tpSEIR[,5]	
			N<-S+E+I+R
		}
		else{
			N<-N; 	S<-S; E<-E; I<-I; R<-R
		}
		
		#append
		if(!is.integer0(grep(object@type,"stochastic")) && !is.integer0(grep(type,"deterministic"))){
			newObj <- seir(type=type,duration=duration,unitTIME=unitTIME,
					S=S,E=E,I=I,R=R,N=N,
					mu=mu,beta0=beta0,beta1=beta1,sigma=sigma,gamma=gamma,T=T,phi=phi,	
					nbVilles=nbVilles,seed=seed)

			newPops<-newObj@pop
			for(i in 1:nbVilles){
				newpopi<-newPops[[i]]
				newpopi$time<-newpopi$time + tpSEIR[1,1]
				newPops[[i]]<-newpopi	
			}			
			if(!append){ 
				newObj@pop<-newPops								
				return(newObj)
			}
			else{	
					
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
		else
		if(!is.integer0(grep(object@type,"deterministic")) && !is.integer0(grep(type,"stochastic"))){
			newObj <- seir(type=type,duration=duration,unitTIME=unitTIME,S=S,E=E,I=I,R=R,N=N,
				mu=mu,beta0=beta0,beta1=beta1,sigma=sigma,gamma=gamma,T=T,phi=phi,
				nbVilles=nbVilles,epsilon=epsilon,rho=rho,seed=seed,rng=rng)
			newPops <- newObj@pop	
			for(i in 1:nbVilles){	
				newpopi <- newPops[[i]]
				newpopi$time <- newpopi$time + tpSEIR[1,1]
				newPops[[i]]<-newpopi
			}	

			if(!append){
				newObj@pop<-newPops
				 return(newObj)
			}
			else{
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
		else
		if(!is.integer0(grep(object@type,"deterministic")) && !is.integer0(grep(type,"deterministic"))){			
			newObj <- seir(type=type,duration=duration,unitTIME=unitTIME,
						S=S,E=E,I=I,R=R,N=N,
						mu=mu,beta0=beta0,beta1=beta1,sigma=sigma,gamma=gamma,T=T,phi=phi[i],	
						nbVilles=nbVilles,seed=seed)
			newPops<-newObj@pop	
			for(i in 1:nbVilles){
				newpopi<-newPops[[i]]
				newpopi$time<-newpopi$time + tpSEIR[1,1]
				newPops[[i]]<-newpopi
			}
			if(!append){	
				newObj@pop<-newPops
				return(newObj)
			}
			else{				
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
		if(!is.integer0(grep(object@type,"stochastic")) && !is.integer0(grep(type,"stochastic"))){
			newObj <- seir(type=type,duration=duration,unitTIME=unitTIME,S=S,E=E,I=I,R=R,N=N,
				mu=mu,beta0=beta0,beta1=beta1,sigma=sigma,gamma=gamma,T=T,phi=phi,
				nbVilles=nbVilles,epsilon=epsilon,rho=rho,seed=seed,rng=rng)
			newPops <- newObj@pop
			for(i in 1:nbVilles){
				newpopi <- newPops[[i]]
				newpopi$time <- newpopi$time + tpSEIR[1,1]
				newPops[[i]]<-newpopi			
			}
			if(!append){
				newObj@pop<-newPops
				return(newObj)
			}
			else{
				
				seqvil<-resizeVector(seq(1:nbVilObj),nbVilles)
				for(i in 1:nbVilles){
					oldpopi <- oldPops[[seqvil[i]]]
					if(!is.null(t0)) oldpopi <- subset(oldpopi,oldpopi$time<=t0)
					newPops[[i]] <- rbind(oldpopi,newPops[[i]])	
				}				
						
					newObj@pop <- newPops
			}				
		}
		else
			stop("Error, type should be 'stochastic' or 'deterministic'.")
		newObj@duration<-objtmp@duration
		newObj@mu<-objtmp@mu	
		newObj@beta0<-objtmp@beta0	
		newObj@beta1<-objtmp@beta1
		newObj@sigma<-objtmp@sigma
		newObj@gamma<-objtmp@gamma
		newObj@phi<-objtmp@phi				
		return(newObj)	
	}
	
}

