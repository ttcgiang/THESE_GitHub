#####################################################################################
pop <- function (object,...) 
  UseMethod("pop")

pop.seir <- function(object,subset=c(),fct=NULL,...){
	pops.stoch <- object@pop
	nbVilles <- length(pops.stoch)
	if(nbVilles == 1) return(pops.stoch)
	else
	if(nbVilles >1){
		if(is.null(subset)){
			pops.sto.subset <- pops.stoch
		}
		else{
			pops.sto.subset <-  vector("list")
			namepop<-c()
			nbvalidvil <- 0; nsubset <-length(subset)
			#checking city at subset[i] existing or non in popsall list
			for(i in 1:nsubset){			
				ixvil <- which(c(1:nbVilles)==subset[i])
				#print(ixvil)			
				if(!is.integer0(ixvil)){
					namepop <- cbind(namepop,paste("pop",ixvil))
					nbvalidvil <- nbvalidvil+1
					pops.sto.subset[[nbvalidvil]] <- pops.stoch[[ixvil]]				
				}
				else
				print("Invalid city index in subset!")			
			}
			names(pops.sto.subset)<-namepop				
		}
		#working with fct
		if(is.null(fct)){#returning all 
			return(pops.sto.subset)
		}
		else{#returning a seul dataframe according to the parametre "fct"
			#rejecting the parametre "time" in list		
			tmp <- pops.sto.subset
			# le nombre de populations est
			nbpop <- length(tmp)
			# la longueur est
			l <- nrow(tmp[[1]])
			# le nombre de variable est
			nbvar <- ncol(tmp[[1]])
			# Pour transformer ta liste en array:
			tmp <- array(unlist(tmp),c(l,nbvar,nbpop))
			# Pour appliquer n'importe fonction "fct" à travers toutes les populations:
			#print(class(apply(tmp,c(1,2),fct)))
			return(apply(tmp,c(1,2),fct))				
		}
	}
	else
		stop("invalid number of cities!")	
}
# printing data of an object of seir.stoch
########################
'+.seir'<-function(a=object, b=object){
	cat("addition of two seir objects.\n")
}

'-.seir'<-function(a=object, b=object){
	cat("subtraction of two seir objects.\n")
}

'*.seir'<-function(a=object, b=object){
	cat("multiplication of two seir objects.\n")
}

'/.seir'<-function(a=object, b=object){
	cat("division of two seir objects.\n")
}

'^.seir'<-function(a=object, b=object){
	cat("exponentiation of two seir objects.\n")
}

'**.seir'<-function(a=object, b=object){
	cat("exponentiation of two seir objects.\n")
}

'%%.seir'<-function(a=object, b=object){
	cat("modulus of two seir objects.\n")
}

'%/%.seir'<-function(a=object, b=object){
	cat("integer division of two seir objects.\n")
}

'<.seir'<-function(a=object, b=object){
	cat("'less than' operator of two seir objects.\n")
}

'<=.seir'<-function(a=object, b=object){
	cat("'less than or equal to' operator of two seir objects.\n")
}

'>.seir'<-function(a=object, b=object){
	cat("'more than or equal to' operator of two seir objects.\n")
}

'>=.seir'<-function(a=object, b=object){
	cat("'more than or equal to' operator of two seir objects.\n")
}

'==.seir'<-function(a=object, b=object){
	cat("'exactly equal to' operator of two seir objects.\n")
}

'!=.seir'<-function(a=object, b=object){
	cat("'not equal to' operator of two seir objects.\n")
}

'!.seir'<-function(a=object){
	if(is.null(a)) return(TRUE)
	else return(FALSE)
}

'|.seir'<-function(a=object, b=object){
	if(is.null(a)&&is.null(b)) return(FALSE)
	else return(TRUE)
}

'&.seir'<-function(a=object, b=object){
	if(!is.null(a)&& !is.null(b)) return(TRUE)
	else return(FALSE)
}


########################
paramduration<-function(ville=1,nbVilles=2,param=c(0.2,0.2,0.1,0.1,0.3,0.3),duration=c(5,10,8)){
	nbfois<-length(param)/nbVilles
	pamvil<-c()
	for(i in 0:(nbfois-1)) pamvil<-c(pamvil,param[ville+i*nbVilles]) 
	
	newpamvil<-c()
	newdurvil<-c()
	if(length(pamvil)>1){
		for(i in 1:(length(pamvil)-1)){
			if(pamvil[i]!=pamvil[i+1]){
				newpamvil<-c(newpamvil,pamvil[i])
				newdurvil<-c(newdurvil,sum(duration[1:i]))
			}
			#else newdurvil<-c(newdurvil,sum(duration[i],duration[i+1]))

			if(i==(length(pamvil)-1)&&(pamvil[i]!=pamvil[i+1])){
				newpamvil<-c(newpamvil,pamvil[i+1])
				newdurvil<-c(newdurvil,sum(duration[1:(i+1)]))
			}
			if(i==(length(pamvil)-1)&&(pamvil[i]==pamvil[i+1])){
				newpamvil<-c(newpamvil,pamvil[i])
				 newdurvil<-c(newdurvil,sum(duration[1:(i+1)]))
			}	
		}
	
	}
	else{
		newpamvil<-pamvil
		newdurvil<-duration
	}
	return(list(newpamvil,newdurvil))
}
###########
summary.seir<-function(object,...)
{
	pops.stoch <- object@pop
	nbVilles <- object@nbVilles

	cat(paste("Type:",object@type),fill=60)	
	cat(paste("Method:",object@method),fill=60)
	cat(paste("RNG: built-in C++"),fill=60)
	cat(paste("Duration:",sum(object@duration),"days"),fill=60)
	unitTIME <- object@unitTIME
	if(unitTIME<=1) cat(paste("Time unit:",unitTIME,"day"),fill=60)
	else cat(paste("Time unit:",unitTIME,"days"),fill=60)
	cat(paste("Seed:",round(object@seed,digits=5)),fill=60)
	cat(paste("Number of populations:",object@nbVilles),fill=60)
	unityear<-365
	
	for(i in 1: nbVilles){
		cat(paste("$-----pop",i),fill=60)
		cat("Coefficients:\n")
		cat(paste(format(paste("Population size N of population ",i,":"),width=50),format(object@N[i],width=10)),fill=60)
		#beta0
		beta0dura<-paramduration(ville=i,nbVilles=nbVilles,object@beta0,object@duration)		
		vilbeta0<-beta0dura[[1]]
		vilduration<-round(beta0dura[[2]]/unityear,digits=2)
		nbsimulbeta0<-length(vilbeta0)
		if(nbsimulbeta0==1){
			cat(paste(format("Mean contact rate beta0:",width=50),format(round(vilbeta0,digits=5),width=10),
					"/day/ind",sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulbeta0){
			if(j==1) {
				    cat(paste(format("Mean contact rate beta0:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(paste(round(vilbeta0[j],digits=5),"/day/ind",sep=""),width=22)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(paste(round(vilbeta0[j],digits=5),"/day/ind",sep=""),width=22)),fill=60)
			}
			}
		}	
		#cat(paste("  amplitude of  contact rate beta1:",object@beta1),fill=60)
		beta1dura<-paramduration(ville=i,nbVilles=nbVilles,object@beta1,object@duration)		
		vilbeta1<-beta1dura[[1]]
		vilduration<-round(beta1dura[[2]]/unityear,digits=2)
		nbsimulbeta1<-length(vilbeta1)
		if(nbsimulbeta1==1){
			cat(paste(format("Amplitude of contact rate beta1:",width=50),format(round(vilbeta1,digits=5),width=10),
				sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulbeta1){
			if(j==1) {
				    cat(paste(format("Amplitude of contact rate beta1:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(round(vilbeta1[j],digits=5),width=7)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(round(vilbeta1[j],digits=5),width=7)),fill=60)
			}
			}
		}	
		#cat(paste("  latency rate sigma:",object@sigma,"/day (mean latency duration:",1/object@sigma,"days"),fill=60)
		sigmadura<-paramduration(ville=i,nbVilles=nbVilles,object@sigma,object@duration)		
		vilsigma<-sigmadura[[1]]
		vilduration<-round(sigmadura[[2]]/unityear,digits=2)
		nbsimulsigma<-length(vilsigma)
		if(nbsimulsigma==1){
			cat(paste(format("Latency rate sigma:",width=50),format(round(vilsigma,digits=5),width=10),
					"/day\n(mean duration of latency period: ",1/vilsigma," days)",sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulsigma){
			if(j==1) {
				    cat(paste(format("Latency rate sigma:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(paste(round(vilsigma[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(paste(round(vilsigma[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			}
		}	
		#cat(paste("  recovery rate gamma:",object@gamma,"/day (mean infection duration:",1/object@gamma,"days"),fill=60)
		gammadura<-paramduration(ville=i,nbVilles=nbVilles,object@gamma,object@duration)		
		vilgamma<-gammadura[[1]]
		vilduration<-round(gammadura[[2]]/unityear,digits=2)
		nbsimulgamma<-length(vilgamma)
		if(nbsimulgamma==1){
			cat(paste(format("Recovery rate gamma:",width=50),format(round(vilgamma,digits=5),width=10),
					"/day\n(Mean duration of infectious period: ",1/vilgamma," days)",sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulgamma){
			if(j==1) {
				    cat(paste(format("Recovery rate gamma:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(paste(round(vilgamma[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(paste(round(vilgamma[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			}
		}	
		#cat(paste("  death rate mu:",object@mu,"/day (mean life expectancy :",1/(object@mu *365),"years"),fill=60)
		mudura<-paramduration(ville=i,nbVilles=nbVilles,object@mu,object@duration)		
		vilmu<-mudura[[1]]
		vilduration<-round(mudura[[2]]/unityear,digits=2)
		nbsimulmu<-length(vilmu)
		if(nbsimulmu==1){
			cat(paste(format("Birth and death rate mu:",width=50),format(round(vilmu,digits=7),width=10),
					"/day\n(Mean life expectancy: ",1/(vilmu*365)," years)",sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulmu){
			if(j==1) {
				    cat(paste(format("Birth and death rate mu:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(paste(round(vilmu[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(paste(round(vilmu[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			}
		}	
		#synchrony phi:
		phidura<-paramduration(ville=i,nbVilles=nbVilles,object@phi,object@duration)		
		vilphi<-phidura[[1]]
		vilduration<-round(phidura[[2]]/unityear,digits=2)
		nbsimulphi<-length(vilphi)
		if(nbsimulphi==1){
			cat(paste(format("Synchrony phi:",width=50),format(round(vilphi,digits=5),width=10),
					sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulphi){
			if(j==1) {
				    cat(paste(format("Synchrony phi:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(paste(round(vilphi[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(paste(round(vilphi[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			}
		}	
		# State variables:
		cat("\nState variables:",fill=60)			
		popi <- pops.stoch[[i]]
		nbcol <- ncol(popi)
		popi <-popi[,c(2:nbcol)]
		print(apply(popi,2,summary))
	}	
}
#####################################################################################
# printing an object of seir.stoch
print.seir <- function(object,...) {
#printing the values of the populations as "time	S	E	P	R	N	I"
	object
   }

#####################################################################################
# printing coef of an object of seir.stoch
coef.seir<-function(object,...){
	#cat("\nPrinting the values of parameters of an object in the seir.stoch class\n")
	cat("@initial values of variables:",fill=60)
	varis <- c(N = as.numeric(object@N),
			S = as.numeric(object@S),
			E = as.numeric(object@E),
			I = as.numeric(object@I),
			R = as.numeric(object@R)
			)
	print(varis)
	cat("\n@values of parameters :",fill=60)
	parameters <- c(duration = object@duration,
			unitTIME = object@unitTIME,
			nbVilles = object@nbVilles,			
			T = object@T,
			mu = object@mu,
			beta0 = object@beta0,
			beta1 = object@beta1,
			sigma = object@sigma,
			gamma = object@gamma,
			epsilon = object@epsilon,
			rho = object@rho,
			phi = object@phi,
			seed = object@seed)
	print(parameters)
#	cat(paste("rng:",object@rng))
  }
########################################
str.seir<-function(object,...){
	cat("@populations:",fill=60)
	nbVilles <- object@nbVilles
	for(i in 1:nbVilles){	
		cat(paste("..$pop",i,": data frame"),fill=60)	
		cat("..$coef:",fill=60)	
		N<-object@N[i]
		S<-object@S[i]
		E<-object@E[i]
		I<-object@I[i]
		R<-object@R[i]
		duration<-sum(object@duration)
		unitTIME<-object@unitTIME
		#beta0		
		beta0dura<-paramduration(ville=i,nbVilles=nbVilles,object@beta0,object@duration)		
		vilbeta0<-beta0dura[[1]]
		vilduration<-beta0dura[[2]]
		nbsimulbeta0<-length(vilbeta0)
		if(nbsimulbeta0==1){
			beta0<-round(vilbeta0,digits=5)
		}
		else{
			beta0<-c()
			for(j in 1:nbsimulbeta0){
				beta0<-c(beta0,round(vilbeta0[j],digits=5))
			}
		}	
		#beta1
		beta1dura<-paramduration(ville=i,nbVilles=nbVilles,object@beta1,object@duration)		
		vilbeta1<-beta1dura[[1]]
		vilduration<-beta1dura[[2]]
		nbsimulbeta1<-length(vilbeta1)
		if(nbsimulbeta1==1){
			beta1<-round(vilbeta1,digits=5)
		}
		else{
			beta1<-c()
			for(j in 1:nbsimulbeta1){
				beta1<-c(beta1,round(vilbeta1[j],digits=5))
			}
		}
		#cat(paste("  latency rate sigma:",object@sigma,"/day (mean latency duration:",1/object@sigma,"days"),fill=60)
		sigmadura<-paramduration(ville=i,nbVilles=nbVilles,object@sigma,object@duration)		
		vilsigma<-sigmadura[[1]]
		vilduration<-sigmadura[[2]]
		nbsimulsigma<-length(vilsigma)
		if(nbsimulsigma==1){
			sigma<-round(vilsigma,digits=5)
		}
		else{
			sigma<-c()
			for(j in 1:nbsimulsigma){
				sigma<-c(sigma,round(vilsigma[j],digits=5))
			}
		}
		#gamma				
		gammadura<-paramduration(ville=i,nbVilles=nbVilles,object@gamma,object@duration)		
		vilgamma<-gammadura[[1]]
		vilduration<-gammadura[[2]]
		nbsimulgamma<-length(vilgamma)
		if(nbsimulgamma==1){
			gamma<-round(vilgamma,digits=5)
		}
		else{
			gamma<-c()
			for(j in 1:nbsimulgamma){
				gamma<-c(gamma,round(vilgamma[j],digits=5))
			}
		}
		#mu			
		mudura<-paramduration(ville=i,nbVilles=nbVilles,object@mu,object@duration)		
		vilmu<-mudura[[1]]
		vilduration<-mudura[[2]]
		nbsimulmu<-length(vilmu)
		if(nbsimulmu==1){
			mu<-round(vilmu,digits=5)
		}
		else{
			mu<-c()
			for(j in 1:nbsimulmu){
				mu<-c(mu,round(vilmu[j],digits=5))
			}
		}	
		#phi			
		phidura<-paramduration(ville=i,nbVilles=nbVilles,object@phi,object@duration)		
		vilphi<-phidura[[1]]
		vilduration<-phidura[[2]]
		nbsimulphi<-length(vilphi)
		if(nbsimulphi==1){
			phi<-round(vilphi,digits=5)
		}
		else{
			phi<-c()
			for(j in 1:nbsimulphi){
				phi<-c(phi,round(vilphi[j],digits=5))
			}
		}	
		######
		parami<-c(N=N,S=S,E=E,I=I,R=R,duration=duration,unitTIME=unitTIME,beta0=beta0,beta1=beta1,sigma=sigma,gamma=gamma,mu=mu,phi=phi)
		names(parami)[1:5]<-c("N","S","E","I","R")	
		print(parami)
		if(i!=nbVilles) cat("\n",fill=60)
	}

	
	#coef(object)
}
########################################

