source("expTroisSupHIS.R")
###
##############
# exploiting multi-N and multi-rho for the relation between "estimated persistence rate" and "phiMAX"
# nbVilles, fixed
para_expvecVIL<-function(nbRep=10,nbPARA=10,duration=50*365,
	nbCONTACT0=100,nbCONTACT1=0.10,probVISITER=0.0,probINFECTER=0.01,statSTATE=TRUE,
	phiMAX=0,N=1e6,vecVil=c(19,20,25,30,35,40)){	
	# for each pair	
	for(i in 1:length(vecVil)){
		print(paste("nbVil=",vecVil[i], "phiMAX=",phiMAX, "probVIS=",probVISITER))
		print(paste("So lan lap=",nbRep*nbPARA))
		
		# step 1:  doing simulations (100 fois)
		# step 2:  finding the extinct rate
		rateDATA<-para.RATEEXT(nbRep=nbRep,nbPARA=nbPARA,duration=duration,
		nbVilles=vecVil[i],N=N,nbCONTACT0=nbCONTACT0,			nbCONTACT1=nbCONTACT1,probVISITER=probVISITER,probINFECTER=probINFECTER,phiMAX=phiMAX,statSTATE=statSTATE)
		if(i==1) resRATE<-rateDATA[[1]]
		else 	resRATE<-rbind(resRATE,rateDATA[[1]])
		
			
	}
	# table of #nbVilles #N #nbCONTACT0 #nbCONTACT1 #nbMulCONTACT #probVISITER #phiMAX #rateEXT #timeSIMU
	#total	
	nameRes<-paste("DATAMETAPOP/resMETAN",N,"probVIS",probVISITER,"phiMAX",phiMAX,".Rda",sep="")
	save(resRATE,file=nameRes)
	return(resRATE)
}
################
filtrer2Phase <- function(nbVilles=8,probVISITER=0.1,vecphiMAX=c(0,pi),N=1e6,duration=50*365)
{
	#minimum
	resrhoDIFphiPI<-para.RATEEXT(nbVilles=nbVilles,probVISITER=probVISITER,phiMAX=pi,N=N)
	raterhoDIFphiPI<-resrhoDIFphiPI[[1]][,12]
	### medium
	tp=-0.1
	while(tp <raterhoDIFphiPI){
	resrhoDIFphi0<-para.RATEEXT(nbVilles=nbVilles,probVISITER=probVISITER,phiMAX=0,N=N,
			duration=duration)
	tp<-resrhoDIFphi0[[1]][,12]
	}

	resSAVE<-list(resrhoDIFphi0,resrhoDIFphiPI)
	namefile=paste("resMETAvil",nbVilles,"rho",probVISITER,"N",N,".Rda",sep="")
	save(resSAVE,file=namefile)
	print(resSAVE)
	return(resSAVE)
}
####
