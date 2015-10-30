 setClass("seir",
   representation(pop="list",duration="numeric",S="numeric",E="numeric",I="numeric",R="numeric",N="numeric",T="numeric",
		mu="numeric",beta0="numeric",beta1="numeric",sigma="numeric",gamma="numeric",
		unitTIME="numeric",phi="numeric",nbVilles="numeric",type="character",
		epsilon="numeric",rho="numeric",seed="numeric",rng="character",method="character",persistence="data.frame"),

   prototype(duration=1,S=NULL,E=NULL,I=NULL,R=NULL,N=1e7,T=365,mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/7,gamma=1/7,pop=NULL,
		unitTIME=1,phi=0,nbVilles=1,epsilon=0.0,rho=0.0,
		seed=as.numeric(Sys.time()),rng="good",method="direct")
)

###################################################
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}
####
timeseq <- function(duration=10*365, unitTIME=1){
	#dayDURA <- duration*365
	return(seq(0,duration,le=duration/unitTIME))
}
#####
phiseq <- function(phiMIN=0, phiMAX=0, nbVilles=1)
	 return(seq(from=phiMIN,to=phiMAX,length.out=nbVilles))
##pause()
pause <- function() readline("press Enter to continue!")
#####
# resize un vector with length given
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
####
initVariables <- function(nbVilles=1,S=NULL,E=NULL,I=NULL,R=NULL,N=NULL,
				mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/8,gamma=1/5,T=365,phi=0){
	S<-resizeVector(V=S,n=nbVilles)
	E<-resizeVector(V=E,n=nbVilles)
	I<-resizeVector(V=I,n=nbVilles)
	R<-resizeVector(V=R,n=nbVilles)
	N<-resizeVector(V=N,n=nbVilles)
	
	vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
	nbVar <- length(vecVar)/nbVilles

	if(nbVar==0) stop("There are not any initial variables!")
	else if(nbVar==5){
		#print("les 5 arguments sont spécifiés")
		validN = TRUE
		for(i in 1:nbVilles){
			if(N[i]!=(S[i]+E[i]+I[i]+R[i])) validN <- FALSE
		}
		if(validN) return(vecVar)
		else
			stop("Values of initial variables should verify N = S + E + I + R!")

	}
	else if(nbVar==4){
		#print("les 4 arguments sont spécifiés")
		if(is.integer0(which(vecVar==S))){
			#print("missing  S")
			S<-N-E-I-R	
			vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
			return(vecVar)		
		}
		else if(is.integer0(which(vecVar==E))){
			#print("missing  E")
			E<-N-S-I-R	
			vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
			return(vecVar)		
		}
		else if(is.integer0(which(vecVar==I))){
			#print("missing  I")
			I<-N-S-E-R	
			vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
			return(vecVar)		
		}
		else if(is.integer0(which(vecVar==R))){
			#print("missing  R")
			R<-N-S-I-E	
			vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
			return(vecVar)		
		}
		else if(is.integer0(which(vecVar==N))){
			#print("missing  N")
			N<-S+E+I+R
			vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
			return(vecVar)		
		}
	}
	else if(nbVar<4){
		if((nbVar==1)&(!is.null(N))){
			S<-c(); E<-c(); I<-c(); R<-c()
			for(i in 1:nbVilles){				
				valEqual <- equilibrium(duration=100*365,N=N[i],mu=mu[i],beta0=beta0[i],beta1=beta1[i],sigma=sigma[i],gamma=gamma[i],phi=phi[i],T=T)
				Sdet <- valEqual[1]; S = c(S,round(Sdet));
				Edet <- valEqual[2]; E = c(E,round(Edet));
				Idet <- valEqual[3]; I = c(I,round(Idet));
				R <- c(R,N[i] - round(Sdet) - round(Edet) - round(Idet)); 
			}
			vecVar<-c(S=S,E=E,I=I,R=R,N=N)	
			return(vecVar)			
		}
		else
			stop("Number of initial variables should be more than 4!")
	}

}
###################
#forme S, E, P, N, I
formeTSEPRNI<-function(nbVilles=1,pop=list(),unitTIME=1)
{
	nbPops<-length(pop)
	if(nbPops!=nbVilles) stop("Error, nbVilles shouble be the same as the number of populations.")
	else{
		flagTSEPRNI<-TRUE;
		for(i in 1:nbPops)
		{
			if(ncol(pop[[i]])!=7){
				flagTSEPRNI<-FALSE;
				break;
			}
		}	
		if(flagTSEPRNI) return(pop)
		else{
			for(i in 1:nbPops)
			{
				tpTIME<-0
				tmpopi<-data.frame()
				popi<-pop[[i]]
				nbrow<-nrow(popi)
				for(j in 1:nbrow)
				{
					
					if(j==1){
						tmpopi<-rbind(tmpopi,c(time=popi[1,1],S=popi[1,2],E=popi[1,3],P=popi[1,4],
								R=popi[1,5],N=popi[1,5],I=0))
						oldj<-1
					}
					else{
						if(popi[j,1]>=tpTIME)
							tptotal<-0
							subpopi<-popi[oldj:j,]
							for(k in 2:nrow(subpopi)) tptotal<-tptotal+abs(subpopi[k,4]-subpopi[(k-1),4])
							tmpopi<-rbind(tmpopi,c(time=tpTIME,S=popi[j,2],E=popi[j,3],P=popi[j,4],
								R=popi[j,5],N=popi[j,6],I=tptotal))
							oldj<-j
					}			
					tpTIME=tpTIME+unitTIME
				}
				pop[[i]]<-tmpopi	
				names(pop[[i]]) <- c("time","S","E","P","R","N","I")
			}
			return(pop)
		}
	}
}
#################
#Basic function
seir<-function(type="stoch",duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/8,gamma=1/5,
			T=365,phi=0,nbVilles=1,epsilon=0.0,rho=0.0,seed=as.numeric(Sys.time()),rng="good",
			S=NULL,E=NULL,I=NULL,R=NULL,N=1e7,...){

	if(!is.integer0(grep(type,"stochastic")))
		seir.stoch(method=method,duration=duration,unitTIME=unitTIME,mu=mu,beta0=beta0,beta1=beta1,
			sigma=sigma,gamma=gamma,T=T,phi=phi,nbVilles=nbVilles,epsilon=epsilon,rho=rho,
			seed=seed,rng=rng,S=S,E=E,I=I,R=R,N=N)
	else if(!is.integer0(grep(type,"deterministic")))
		seir.det(nbVilles=nbVilles,duration=duration,unitTIME=unitTIME,S=S,E=E,I=I,R=R,N=N,T=T,
			mu=mu,beta0=beta0,beta1=beta1,sigma=sigma,gamma=gamma,phi=phi)
	else
		stop("invalid value of type!")
}
#seir.det() produit un object de seir/sir en déterminist
seir.det<-function(nbVilles=1,duration=10*365,unitTIME=1,S=NULL,E=NULL,I=NULL,R=NULL,N=1e7,T=365,
			mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/8,gamma=1/5,phi=0,...){
	#parametres
	sigma<-resizeVector(sigma,nbVilles);
	for(i in 1:length(sigma))
		if(sigma[i]<0) stop("Error, invalid value of sigma, sigma should be more than zero.")
	gamma<-resizeVector(gamma,nbVilles);
	for(i in 1:length(gamma))
		if(gamma[i]<0) stop("Error, invalid value of gamma, gamma should be more than zero.")
	mu<-resizeVector(mu,nbVilles);
	for(i in 1:length(mu))
		if(mu[i]<0) stop("Error, invalid value of mu, mu should be more than zero.")
	beta0<-resizeVector(beta0,nbVilles);
	for(i in 1:length(beta0))
		if(beta0[i]<0) stop("Error, invalid value of beta0, beta0 should be more than zero.")
	beta1<-resizeVector(beta1,nbVilles);
	for(i in 1:length(beta1))
		if(beta1[i]<0) stop("Error, invalid value of beta1, beta1 should be more than zero.")
	phi<-resizeVector(phi,nbVilles);

	#init.varib	
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
	times <- timeseq(duration,unitTIME)
	
	# Solving the system of differential equations (with deSolve package):
	populations<-vector("list",nbVilles)
	for(i in 1:nbVilles){
		pars<-c(beta0=beta0[i],beta1=beta1[i],T=T,phi=phi[i],mu=mu[i],sigma=sigma[i],gamma=gamma[i],N=N[i])
		yini<-c(t=times[1],S=S[i],E=E[i],I=I[i])
		#print(yini)
		pops.det<-as.data.frame(ode(func = "derivs", y = yini, parms = pars, times = times, dllname = "dizzys", initfunc = "initparms"))
		valueN<-rep(N[i],nrow(pops.det))
		#print(head(pops.det))
		valueR<-valueN-pops.det[,3]-pops.det[,4]-pops.det[,5]
		pops.det<-data.frame(pops.det[,1], pops.det[,3],pops.det[,4],pops.det[,5],valueR,valueN)
		#print(head(pops.det))
		names(pops.det)<-c("time","S","E","P","R","N")		
		names(populations)<-paste("pop",i)
		populations[[i]]<-pops.det
	}
	#populations<-formeTSEPRNI(1,populations,unitTIME)
	#return a seir object
	return(new("seir",nbVilles=nbVilles,duration=duration,unitTIME=unitTIME,
				N=N,S=S,E=E,I=I,R=R,T=T,mu=mu,beta0=beta0,beta1=beta1,
					sigma=sigma,gamma=gamma,phi=phi,pop=populations,type="deterministic"))
}
######
#seir.stoch produit un object de seir en stochastic
#rng has two names, "good" or "fast"
#method also has two names, "direct" or "adaptivetau"
#phi=[phiMIN,phiMAX]
seir.stoch<-function(duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/8,gamma=1/5,
			T=365,phi=0,nbVilles=1,epsilon=0.0,rho=0.0,seed=as.numeric(Sys.time()),rng="good",
			S=NULL,E=NULL,I=NULL,R=NULL,N=1e7,...){
	if(is.null(phi)) stop("invalid value of phi!")
	else		phi<-resizeVector(phi,nbVilles)
	#parametres
	#sigma
	sigma<-resizeVector(sigma,nbVilles);
	for(i in 1:length(sigma))
		if(sigma[i]<0) stop("Error, invalid value of sigma, sigma should be more than zero.")
	#gamma
	gamma<-resizeVector(gamma,nbVilles);
	for(i in 1:length(gamma))
		if(gamma[i]<0) stop("Error, invalid value of gamma, gamma should be more than zero.")
	#mu
	mu<-resizeVector(mu,nbVilles);
	for(i in 1:length(mu))
		if(mu[i]<0) stop("Error, invalid value of mu, mu should be more than zero.")
	#beta0
	beta0<-resizeVector(beta0,nbVilles);
	for(i in 1:length(beta0))
		if(beta0[i]<0) stop("Error, invalid value of beta0, beta0 should be more than zero.")
	#beta1
	beta1<-resizeVector(beta1,nbVilles);
	for(i in 1:length(beta1))
		if(beta1[i]<0) stop("Error, invalid value of beta1, beta1 should be more than zero.")
	#
	phi<-resizeVector(phi,nbVilles);
	arr_rho<-matrix(rho[1])
	arr_epsilon<-matrix(epsilon[1])
	tmax<-duration
	
	#init.varib	
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
	#simuler n populations
	# for the direct algorithm
	if(!is.integer0(grep(method,"direct"))){

		if(!is.integer0(grep(rng,"good"))) typeRNG=0
		else
		  if(!is.integer0(grep(rng,"fast"))) typeRNG=1
		else
			stop("invalid value of rng, it should be 'good' or 'fast'")

		#récupérer les parametres
		pops <- .Call("getValeurPOPS",nbVilles,sigma,gamma,mu,beta0,beta1,phi,arr_rho,arr_epsilon,
						seed,unitTIME,duration,typeRNG,T,init.varib)
		populations <- vector("list",nbVilles)
		names(populations)<-paste("pop",c(1:nbVilles))
		for(i in 1:nbVilles){
			populations[[i]] <- data.frame(pops[i])
			names(populations[[i]]) <- c("time","S","E","P","R","N","I")
		}	
	}
	else
	if(!is.integer0(grep(method,"adaptivetau"))){	
		flagSEIR = TRUE;
		for(i in 1:nbVilles){
			if(is.infinite(sigma[i])){
				flagSEIR=FALSE;
				break;
			}
		}		
		if(flagSEIR){
			init.values <- c(S=S,E=E,I=I,R=R)	
			relratechange <- rep(1,length(init.values))	
			allPops <- .Call('ssesAdaptiveTau',init.values,NULL,NULL, 
			      #params,
			      nbVilles,beta0,beta1,mu,sigma,gamma,phi,arr_rho,arr_epsilon,T,
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
		else{
			init.values <- c(S=S,I=I,R=R)	
			sigma<-resizeVector(Inf,nbVilles);
			relratechange <- rep(1,length(init.values))	
			allPops <- .Call('ssesAdaptiveTau',init.values,NULL,NULL, 
			      #params,
			      nbVilles,beta0,beta1,mu,sigma,gamma,phi,arr_rho,arr_epsilon,T,
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
				mu=mu,beta0=beta0,beta1=beta1,sigma=sigma,gamma=gamma,
				phi=phi,nbVilles=nbVilles,
				epsilon=epsilon,rho=rho,seed=seed,rng=rng,method=method,pop=populations,type="stochastic"))
#the end
}
#######
#Persistence
#KMsurv
#survival
setGeneric("persistence", function(object,...) standardGeneric("persistence"))
setMethod("persistence", "seir",
	function(object,...){

	nbVilles <- object@nbVilles
	populations <- object@pop
	persis <- sort(sapply(populations,function(x)tail(which(x$I>0),1)))
	#print(persis)
	pers <- cbind(time=persis,nbVilles=nbVilles:1)
	#print(pers)
	allVilles <- dimnames(pers)[[1]]
	allVilles <- as.numeric(lapply(strsplit(allVilles," "), function(x) x[2]))
	pers <- cbind(time=persis,nbVilles=nbVilles:1,ndie=1,dieVille=allVilles)
	remain <- pers[,2]-pers[,3]
	pers <- cbind(pers,remain)
	object@persistence <-as.data.frame(pers)
	return(object)
}
)
####################################################################################"
#plot.pers
setGeneric("plot.pers", function(object,...) standardGeneric("plot.pers"))
setMethod("plot.pers", "seir",
	function(object,x="time",y="P",type="s",col="red",xlim=c(),ylim=c(),curvetype="KM",vilabline=c(),...){
	pers <- object@persistence
	nbVilles <- object@nbVilles
	duration <- object@duration
	if(!is.integer0(grep(curvetype,"KM"))){
		pers <- rbind(data.frame(time=c(0),nbVilles=nbVilles,ndie=0,dieVille=0,remain=nbVilles),pers)
		last <- tail(pers,1)
		if(last[1,1]<duration) pers <- rbind(pers,data.frame(time=duration,
			nbVilles=last[1,2]-1,ndie=0,dieVille=0,remain=0))
		plot(pers$time,pers$remain,type=type,col=col,xlab="time",ylab="#populations non extinct")
	}
	else
	if(!is.integer0(grep(curvetype,"population"))){
		plot.seir(object,y=y,type=type,col=col,xlim=xlim,ylim=ylim,...)
		pop<-c(1:nbVilles)
		allVilles <- pers$dieVille
		nbvilabline <- length(vilabline)	
		if(is.null(pop)){
			for(i in 1:nbvilabline)
				for(j in 1:nbVilles)
					if(vilabline[i]==allVilles[j])
						abline(v=pers$time[j])	
		}
		else{
			for(i in 1:nbvilabline){
				if(!is.integer0(which(pop==vilabline[i]))){
					#print(vilabline[i])
					for(j in 1:nbVilles)
						if(vilabline[i]==allVilles[j])
							abline(v=pers$time[j])		
				}
			}
		}
	}
	else
		stop("We have two curve types, curvetype should be 'KM' or 'population' !")
}
)
# Estimating the survival parameters.
# calculing Survival Probability
#
setGeneric("surv.prob", function(object,...) standardGeneric("surv.prob"))
setMethod("surv.prob", "seir",
   function(object,...){	
	pers<-object@persistence
 	my.fit <- survfit(Surv(pers[,1],pers[,3])~1)
	return(my.fit)
  }
)
#
setGeneric("plot.surv.prob", function(object,...) standardGeneric("plot.surv.prob"))
setMethod("plot.surv.prob", "seir",
   function(object,...){
	my.fit <- surv.prob(object)
	survival:::plot.survfit(my.fit, main="Kaplan-Meier estimate with 95% confidence bounds",xlab="time", ylab="survival probability")
}
)
#
setGeneric("surv.cox", function(object, method='efron',...) standardGeneric("surv.cox"))
setMethod("surv.cox", "seir",
   function(object, method='efron',...){
	ps<-object@persistence
	time <-ps[,1]
	ndies <-ps[,3]
	x <-ps[,3]
	phm <- coxph(Surv(time,ndies)~x,method=method)
	return(phm)
}
)
#
#R does not provide baseline hazard
# R provides average hazard
# i.e. for individual with all covariates equal to average x — over whole sample
setGeneric("detail.cox", function(object,...) standardGeneric("detail.cox"))
setMethod("detail.cox", "seir",
   function(object,...){
	phm <- surv.cox(object)
	d.phm <- coxph.detail(phm)
	return(d.phm)
}
)
#
setGeneric("average.hasard", function(object,...) standardGeneric("average.hasard"))
setMethod("average.hasard", "seir",
   function(object,...){
	d.phm <- detail.cox(object)
	meanx <-c(mean(d.phm$x, na.rm=T))
	x <-c(0) - meanx
	return(x)
 }
)
#
#The basic Cox PH model attempts to fit survival data with covariates z to a hazard function of the form
#h(t|z) = h0(t)exp{β'z}
#where β is an unknown vector
setGeneric("surv.beta", function(object,...) standardGeneric("surv.beta"))
setMethod("surv.beta", "seir",
   function(object,...){
	phm <- surv.cox(object)
	beta <-phm$coef
	return(beta)
 }
)
#
#h0 (t) is the baseline hazard, which is non-parametric
setGeneric("baseline.hazard", function(object,...) standardGeneric("baseline.hazard"))
setMethod("baseline.hazard", "seir",
   function(object,...){
	d.phm <- detail.cox(object)
	h0 <- c(0,d.phm$hazard)
	return(h0)
 }
)
#
setGeneric("baseline.surv", function(object,...) standardGeneric("baseline.surv"))
setMethod("baseline.surv", "seir",
   function(object,...){
	h0 <- baseline.hazard(object)
	S0 <-exp(-cumsum(h0))
	return(S0)
 }
)
#If the estimate of the baseline survival function, S0(t), is provided, then the estimate of the
# survival function for an individual with covariates zk may be obtained via
#	S(t) = S0^exp(t(beta) %*% x) 
setGeneric("plot.comp.surv.estim", function(object,...) standardGeneric("plot.comp.surv.estim"))
setMethod("plot.comp.surv.estim", "seir",
   function(object,...){
	d.phm <- detail.cox(object)
	times <-c(0,d.phm$time)	
	S0 <-baseline.surv(object)
	beta <-surv.beta(object)	
	x <-average.hasard(object)
	Sx <-S0^exp(t(beta) %*% x)
	plot.surv.prob(object)
	lines(times,Sx,type="s",col="red", pch=22)
	box()
	legend(100, 0.1, c("old","estimated"), cex=0.8,col=c("black","red"), lty=c(1,1))
 }
)
#############################################################################################

