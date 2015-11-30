#
##
# Library 'parallel' of R that helps us make simulation in parallel
library(parallel)
library("dizzysNewInfec")
library(MASS)
# 
#################
postExtGolVIL<-function(I=c(0,1,0,1))
{
  n=length(I)
  i=1
  post<-0
  if(I[n] ==0){
	  while(i<=n){
	    #print(paste("i=",i))
	    if(I[i]==0){  
	      j<-i+1
	      post<-i
	      while(j<=n){	
	        #   print(paste("j=",j))
	        if(I[j]!=I[i]){          
	          i<-j
		  post<-i
	          break
	        }
	        else{
	          j<-j+1			
		}
	      }
	      i<-j	
	    }
	    else
	      i<-i+1
	  }
	 return(post) 
   }
   else
	return(0)

}

####################################################
# This fucntion does just one simulation (one time) with the values of parameters and variables given
# parameter 'phi' here is a vector of two values, one is for phiMIN and other one is for phiMAX
# Here, we forcus on exploit the parameters:
# 	N: population size of subpopulation
#	nbVilles: number of subpopulation in a metapopulation
#	rho : coupling rate, coupling strength between any two subpopulations
#	phi: lower and upper limits of phase of forcing in the meatapopulation
#	epsilon: fixed is zero, infection rate from outside
# 	duration: time of simulation, 50 years 
# vector total
vecGRAIN <- c()
#
getValExtNEWYANN<- function(nRep=1,N=1e5,nbVilles=3,duration=36500,nbCONTACT0=100,nbCONTACT1=0.10,grain=1090,
	phiMAX=0,probVISITER=0.01,probINFECTER=0.01,statSTATE=TRUE){
	# vector contains time of local disease persistence of each subpopulation in the metapopulation
	# number 'seed'
	#grain=as.numeric(sample(0:100000,1,replace=F))
	#grain=as.numeric(Sys.time())
	grain=runif(1, 1.0, 1000000)
	phiPHASE=seq(0,phiMAX,length=nbVilles)
	#print(phiPHASE)
	#SEIR
	# calling the function 'simul' in the package 'dizzys' to do simulation
	obj<-stoSEIRNewInfec(N=N,nbVilles=nbVilles,duration=duration,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,statSTATE=statSTATE,#mu=1/(70*365),sigma=1/8,gamma=1/5,
phiPHASE=phiPHASE,probVISITER=probVISITER,probINFECTER=probINFECTER,seed=grain)
	#print(obj@N)
	#print(phiPHASE)

	# finding the time of local disease persistences in the metapopulation
	persTimeMETA<-0
	valPersi<-0
	sumInfected<-0
	for(ivil in 1:nbVilles){

		postExtGol<-postExtGolVIL(obj@pop[[ivil]][,4])
		if(postExtGol !=0)
		{
			valPersi<-(obj@pop[[ivil]][postExtGol,1])
			#persistence time
			persTimeMETA<-max(valPersi,persTimeMETA)	
		}
		sumInfected<-sumInfected+sum(obj@pop[[ivil]][,7])
	}
	#save persistence global
	# vecteur about the moment of the beginning point of local extinction
		veclocalExtMETA<-unlist(obj@localExtPOP)
	# vecteur about the moment of the ending point of local extinction
		vecRecolMETA<-unlist(obj@disRecolPOP)	
	res<-list(persTime=persTimeMETA,localEXT=veclocalExtMETA,recolDURA=vecRecolMETA,totalI=sumInfected,valGRAIN=grain)

	return(res)
}
######
######
# This function does many simulations in the parallel form.
# It means to repeat one simulation with M times and here M times of this simulation are parallel
# values=rep(1,M) helps us to makes M simulations in the parallel form
# numWorkers: Number of workers (R processes) used
###
paraGlobSimul<-function(values=rep(1,3),N=1e5,nbVilles=3,duration=36500,nbCONTACT0=100,nbCONTACT1=0.10,
	phiMAX=0,probVISITER=0.01,probINFECTER=0.01,statSTATE=TRUE,numWorkers=100){
	## Number of workers (R processes) to use:
	## Parallel calculation (mclapply):
	res<-mclapply(values, getValExtNEWYANN,
			#parameter
		N=N,nbVilles=nbVilles,duration=duration,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,
		phiMAX=phiMAX,probVISITER=probVISITER,probINFECTER=probINFECTER,statSTATE=statSTATE,
		#paralell
		mc.cores = numWorkers)
	#get the result
	#return(t(as.data.frame(res)))
	return(res)
}
###
###
# this function permits us to do many simulations with the different values of phiMAX
# 	\rho (coupling rate), fixed
#	\nbRep number of simulationis, fixed
#	\nbVilles, number of subpopulations in a metapopulation, fixed
#	\duration, time of each simulation, fixed
#	\epsilon, infection rate from outside, fixed
#	\N, population size of subpopulation, 
#	the same for all subpopulations
para.RATEEXT <-function(nbRep=10,nbPARA=10,N=3e5,nbVilles=3,duration=3*365,nbCONTACT0=100,nbCONTACT1=0.10,
	phiMAX=0,probINFECTER=0.01,probVISITER=0.01,statSTATE=TRUE){
	#simulation (parallel in R)	
	#data.frame  de résultats
	timeSIMUL<-system.time(
	for(i in 1:nbRep){
		resPara<-paraGlobSimul(values=rep(1,nbPARA),			
			#parameter
		N=N,nbVilles=nbVilles,duration=duration,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,
		phiMAX=phiMAX,probVISITER=probVISITER,probINFECTER=probINFECTER,statSTATE=statSTATE)
		if(i==1) dataVAL<-resPara
		else{ 
			lngls1<-length(dataVAL)
			for(j in (lngls1+1):(nbPARA+lngls1)) dataVAL[[j]]<-resPara[[j-lngls1]]
		}
		#print("KKKK")
		#print(resPara[[1]])

	})
	## Dem so phan tu in LIST
	nbElement<-length(dataVAL)
	## Tach du lieu cho time of persistence global
	### tach du lieu cho data.frame #local extinction #dureeRecol
	vecTemPERS<-c()	
	vecTotalLocalExtMETA<-c()
	vecTotalRecolMETA<-c()
	vecsumInfected<-c()
	vecGRAIN <- c()
	for(i in 1:nbElement){
		#si on a l'extinction globale
		if(dataVAL[[i]][[1]] >0) vecTemPERS<-c(vecTemPERS,dataVAL[[i]][[1]])
		vecTotalLocalExtMETA<-c(vecTotalLocalExtMETA,dataVAL[[i]][[2]])
		vecTotalRecolMETA<-c(vecTotalRecolMETA,dataVAL[[i]][[3]])
		#vec SUM INFECTED
		vecsumInfected<-c(vecsumInfected,dataVAL[[i]][[4]])
		#vec GRAIN
		vecGRAIN <- c(vecGRAIN,dataVAL[[i]][[5]])
	}

	#saving the histogram and the density of the "local extinction" and the "time of recolonisation"
	# estimer le taux du nombre d'extinction locale et de durée de recolonization
#peristance globale
	vecPers<-vecTemPERS
	#print("AAAAAAAAAAAAAAAAAAAAAA")
	#print(vecPers)
	if(length(vecPers) !=0){		
		sortVecPers<-sort(vecPers)	
		objEstGolEXT<-survreg(Surv(sortVecPers)~1, dist="exp")
		#confidence rate
		confRate<-confint(objEstGolEXT,level=0.95)	
		perGolLOWER<-confRate[1]
		extGolLOWER<-exp(-perGolLOWER)
		perGolUPPER<-confRate[2]
		extGolUPPER<-exp(-perGolUPPER)
		# estimated rate
		perGolRATE<-objEstGolEXT$coeff[1]
		extGolRATE<-exp(-perGolRATE)
	}
	else{
		#confidence rate
		perGolLOWER<-0
		extGolLOWER<-0
		perGolUPPER<-0
		extGolUPPER<-0
		# estimated rate
		perGolRATE<-0
		extGolRATE<-0
	}
	#result	
	#extinction local
	library(MASS)
	#vecTotalLocalExtMETA<-vecTotalLocalExtMETA
	if(length(vecTotalLocalExtMETA) !=0){
		estObjLocalEXT<-fitdistr(vecTotalLocalExtMETA,densfun="exponential")
		extLocalRATE<-estObjLocalEXT$estimate
		extLocalCONF<-confint(estObjLocalEXT)
		extLocalLOWER<-extLocalCONF[1]
		extLocalUPPER<-extLocalCONF[2]
	}
	else{
		extLocalRATE<-0
		extLocalLOWER<-0
		extLocalUPPER<-0
	}
	##show histogram
	#x <- vecTotalLocalExtMETA
	#h<-hist(x, breaks=10, col="red", xlab="Miles Per Gallon",main=paste("nbVil",nbVilles," rho=",probVISITER," N=",N))
	#xfit<-seq(min(x),max(x),length=40)
	#yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
	#yfit <- yfit*diff(h$mids[1:2])*length(x)
	#lines(xfit, yfit, col="blue", lwd=2) 
	## density
	# Kernel Density Plot
	#d <- density(x) # returns the density data
	#plot(d) # plots the results 
	
	# de durée
	if(length(vecTotalRecolMETA) !=0){
		estObjDurRECL<-fitdistr(vecTotalRecolMETA,densfun="exponential")
		durRecolRATE<-estObjDurRECL$estimate
		estDurRECLCONF<-confint(estObjDurRECL)
		durRecolLOWER<-estDurRECLCONF[1]
		durRecolUPPER<-estDurRECLCONF[2]
	}
	else{
		durRecolRATE<-0
		durRecolLOWER<-0
		durRecolUPPER<-0
	}

	##resultas

	#result
	resline<-c(nbVilles,phiMAX,probVISITER,N,extGolLOWER,extGolRATE,extGolUPPER,perGolLOWER,perGolRATE,perGolUPPER,extLocalLOWER,extLocalRATE,extLocalUPPER,
		durRecolLOWER,durRecolRATE,durRecolUPPER,sum(vecsumInfected),sum(vecTotalLocalExtMETA),
		sum(vecTotalRecolMETA),timeSIMUL[3])
	resRATE<-data.frame(t(resline))
	names(resRATE)<-c("nbVilles","phiMAX","rho","Nmeta","extGolLOWER","extGolRATE","extGolUPPER","perGolLOWER","perGolRATE","perGolUPPER",
		"extLocalLOWER","extLocalRATE","extLocalUPPER","durRecolLOWER","durRecolRATE","durRecolUPPER",
		"sumInfMETA","sumLocalEXT","sumRECOL","timeELAP")
	##rates of metapopuatiopn
	#persistence of time
	#locale extinction et recolonisation
	resRATE<-list(resRATE,persTIME=vecPers,vecGRAIN=vecGRAIN)
	save(resRATE,file=paste("DATAMETAPOP/resVil",nbVilles,"N",N,"probVIS",probVISITER,"phiMAX",phiMAX,".Rda",sep=""))
	return(resRATE)
}
##
##############
##############
# exploiting multi-N and multi-rho for the relation between "estimated persistence rate" and "phiMAX"
# nbVilles, fixed
para_expNmeta<-function(nbRep=10,nbPARA=10,duration=50*365,nbVilles=3,
	nbCONTACT0=100,nbCONTACT1=0.10,probVISITER=0.0,probINFECTER=0.01,statSTATE=TRUE,
	phiMAX=0,vecN=c(1e4,5e4,1e5,3e5,5e5,7e5,1e6)){	
	# for each pair
	print(paste("So lan lap=",nbRep*nbPARA))
	nameHispdf<-paste("figHISnbVil",nbVilles,"probVIS",probVISITER,"phiM",phiMAX,".pdf",sep="")
	pdf(nameHispdf)
		
	for(i in 1:length(vecN)){
		#print(paste("VAL of N=",vecN[i]))
		# step 1:  doing simulations (100 fois)
		# step 2:  finding the extinct rate
		rateDATA<-para.RATEEXT(nbRep=nbRep,nbPARA=nbPARA,duration=duration,
		nbVilles=nbVilles,N=vecN[i],nbCONTACT0=nbCONTACT0,			nbCONTACT1=nbCONTACT1,probVISITER=probVISITER,probINFECTER=probINFECTER,phiMAX=phiMAX,statSTATE=statSTATE)
		if(i==1) resRATE<-rateDATA[[1]]
		else 	resRATE<-rbind(resRATE,rateDATA[[1]])
		
			
	}
	dev.off()
	# table of #nbVilles #N #nbCONTACT0 #nbCONTACT1 #nbMulCONTACT #probVISITER #phiMAX #rateEXT #timeSIMU
	#total	
	nameRes<-paste("DATAMETAPOP/resMETAVil",nbVilles,"probVIS",probVISITER,"phiMAX",phiMAX,".Rda",sep="")
	save(resRATE,file=nameRes)
	return(resRATE)
}
################
nofiltrerRES<-function(nbVilles=3,probVISITER=0.1,vecN=c(1e5,3e5,5e5,7e5,1e6),duration=30*365){
	#
	idctr=12
	idbf=11
	idaft=13
	idctrX=4

	resRATE<-para_expNmeta(nbVilles=nbVilles,probVISITER=0.0,vecN=vecN,duration=duration,phiMAX=0)
	resrho0phi0<-resRATE
	
	#
	
	resRATE<-para_expNmeta(nbVilles=nbVilles,probVISITER=probVISITER,vecN=vecN,duration=duration,phiMAX=0)
	resrhoDIFphi0<-resRATE
	
	#
	
	resRATE<-para_expNmeta(nbVilles=nbVilles,probVISITER=probVISITER,vecN=vecN,duration=duration,phiMAX=pi)
	resrhoDIFphiDIF<-resRATE
	nfile<-paste("fignbVil",nbVilles,"probVIS",probVISITER,".pdf",sep="")

	pdf(nfile)
	maxy<-max(resrhoDIFphi0[,idaft])
	miny<-min(resrho0phi0[,idbf])
	plot(resrhoDIFphi0[,idctrX],resrhoDIFphi0[,idctr],type="b",col="red",lwd=3,ylim=range(miny,maxy),xlab="Size of Metapopulation (#individual)",ylab="Local extinction rate during 50 years simulation time",main="rho=0.1, nbVilles=5, 100 simulations")
	arrows(resrhoDIFphi0[,idctrX],resrhoDIFphi0[,idbf],resrhoDIFphi0[,idctrX],resrhoDIFphi0[,idaft],code=3,angle=90,length=0.05,col="red",lwd=3)

	lines(resrhoDIFphiDIF[,idctrX],resrhoDIFphiDIF[,idctr],type="b",col="blue",lwd=3)
	arrows(resrhoDIFphiDIF[,idctrX],resrhoDIFphiDIF[,idbf],resrhoDIFphiDIF[,idctrX],resrhoDIFphiDIF[,idaft],code=3,angle=90,length=0.05,col="blue",lwd=3)

	lines(resrho0phi0[,idctrX],resrho0phi0[,idctr],type="b",col="black",lwd=3)
	arrows(resrho0phi0[,idctrX],resrho0phi0[,idbf],resrho0phi0[,idctrX],resrho0phi0[,idaft],code=3,angle=90,length=0.05,col="black",lwd=3)
	legend(5e5,maxy/2,c("rho=0,phiMAX=0","rho=0.1,phiMAX=0","rho=0.1,phiMAX=pi"),lty=c(1,1,1), lwd=c(2.5,2.5),col=c("black","red","blue"))

	dev.off()
}

#####

filtrerRES <- function(nbVilles=8,probVISITER=0.1,vecphiMAX=c(0,pi),vecN=c(5e4,1e5,3e5,5e5,7e5),duration=50*365)
{
	#minimum
	resrhoDIFphiPI<-para_expNmeta(nbVilles=nbVilles,probVISITER=probVISITER,phiMAX=pi,vecN=vecN)
	vecRATErhoDIFphiPI<-resrhoDIFphiPI[,12]
	### medium
	dtrhoDIFphi0<-data.frame()
	print(length(vecRATErhoDIFphiPI))
	for(i in 1:length(vecRATErhoDIFphiPI)){
		tp=-0.1
		while(tp <vecRATErhoDIFphiPI[i]){
		 resrhoDIFphi0<-para_expNmeta(nbVilles=nbVilles,probVISITER=probVISITER,phiMAX=0,vecN=vecN[i],
			duration=duration)
			tp<-resrhoDIFphi0[12]
		}
		print(paste("pi=",vecRATErhoDIFphiPI[i],"  phi0=",resrhoDIFphi0[12]))
		print(resrhoDIFphi0)
		if(i==1) dtrhoDIFphi0<-as.data.frame(dtrhoDIFphi0)
		else	dtrhoDIFphi0<-rbind(dtrhoDIFphi0,resrhoDIFphi0)

	}
	vecRATErhoDIFphi0<-dtrhoDIFphi0[,12]
	##maximum
	print(dtrhoDIFphi0)
	dtrho0phi0<-data.frame()
	print(length(vecRATErhoDIFphi0))
	for(j in 1:length(vecRATErhoDIFphi0)){
		tp=-0.1
		while(tp <vecRATErhoDIFphi0[i]){
		 resrho0phi0<-para_expNmeta(nbVilles=nbVilles,probVISITER=probVISITER,phiMAX=0,vecN=vecN[j],
			duration=duration)
			tp<-resrho0phi0[12]
		}
		print(paste("rhoDIFphiDIF=",vecRATErhoDIFphiPI[i],"  rhoDIFphi0=",resrhoDIFphi0[12],"  rho0phi0= ",tp))
		print(resrho0phi0)
		if(i==1) dtrho0phi0<-as.data.frame(resrho0phi0)
		else	dtrho0phi0<-rbind(dtrho0phi0,resrho0phi0)

	}




	resSAVE<-list(resrhoDIFphiPI,dtrhoDIFphi0,resrho0phi0)
	namefile=paste("resMETAvil",nbVilles,"rho",probVISITER,".Rda",sep="")
	save(resSAVE,file=namefile)
	return(resSAVE)
	


}
####

