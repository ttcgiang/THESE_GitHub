#
##
# Library 'parallel' of R that helps us make simulation in parallel
library(parallel)
# our package 'dizzys'
library("dizzysNEWYANN")
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
getValPersNEWYANN<- function(nRep=1,N=1e5,nbVilles=3,duration=100*365,nbCONTACT0=100,nbCONTACT1=0.10,grain=1090,
			 nbMulCONTACT=1,phiMIN=0.0,phiMAX=0.0,probVISITER=0.0,probINFECTER=0.01,namefileMETA="meta_DONNEES.txt"){
	# vector contains time of local disease persistence of each subpopulation in the metapopulation
	# number 'seed'
	#grain=as.numeric(sample(0:100000,1,replace=F))
	#grain=as.numeric(Sys.time())
	grain=runif(1, 1.0, 1000000)
	#SEIR
	# calling the function 'simul' in the package 'dizzys' to do simulation
	obj<-seir.stochYANN(N=N,nbVilles=nbVilles,duration=duration,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,disSTAT=rep(1,nbVilles),#mu=1/(70*365),sigma=1/8,gamma=1/5,
			nbMulCONTACT=nbMulCONTACT,phiMIN=phiMIN,phiMAX=phiMAX,probVISITER=probVISITER,probINFECTER=probINFECTER,seed=grain)
	#plot.seir(obj)
	# finding the time of local disease persistences in the metapopulation
	vecpersTimeMETA<-c()
	for(ivil in 1:nbVilles){
		#postExtGol<-postExtGolVIL(obj@pop[[ivil]][,4])
		postExtGol <- tail(subset(obj@pop[[ivil]],P>0),n=1)[,1]
		vecpersTimeMETA<-c(vecpersTimeMETA,postExtGol)
	}
	##estimer l'extinction locale
	sortVecPers<-sort(vecpersTimeMETA)
	#print(paste("nbVilles = ", nbVilles, " taille de POP=",N, " phiMAX= ",phiMAX, "probVISISIT=", probVISITER))
	#print("PERSISTENCE")
	#print(vecpersTimeMETA)	
	#print(sortVecPers)
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

	##res
	resline<-c(nbVilles,N,phiMAX,extGolLOWER,extGolRATE,extGolUPPER,perGolLOWER,perGolRATE,perGolUPPER)
	#resRATE<-data.frame(t(resline))
	#names(resRATE)<-c("nbVilles","Nmeta","phiMAX","extGolLOWER","extGolRATE","extGolUPPER","perGolLOWER","perGolRATE","perGolUPPER")

	return(resline)
}
######
testCCS<-function(vecN=c(1e5,2e5,3e5,4e5,5e5),vecVIL=c(2,5,10,15),vecphiMAX=c(0,pi),probVISITER=0.01){
	for(phitp in 1:length(vecphiMAX)){	
	for(vil in 1:length(vecVIL)){
	for(i in 1:length(vecN)){
		restp<-getValPersNEWYANN(N=vecN[i],nbVilles=vecVIL[vil],phiMAX=vecphiMAX[phitp],probVISITER=probVISITER)
		#print(restp)
		if(i==1){
			resSAVE<-data.frame(t(restp))
		}
		else
			resSAVE<-rbind(resSAVE,restp)
		#print(resSAVE)
	}
	#print("RESULT")
	#print(resSAVE)
	plot(resSAVE[,2],resSAVE[,5],type="b",col="red",lwd=3,ylab="Size of subPOP",xlab="Extinction rate",main=paste("nbVilles",vecVIL[vil], "phiMAX=",vecphiMAX[phitp]),ylim=range(min(resSAVE[,6]),max(resSAVE[,4])))
	arrows(resSAVE[,2],resSAVE[,4],resSAVE[,2],resSAVE[,6],code=3,angle=90,length=0.05,col="black",lwd=3)
		
	}
	namesRES<-paste("kquanbVilles",vecVIL[vil],"phiMAX",vecphiMAX[phitp],".Rda")
	save(resSAVE,file=namesRES)
	}
	#names(resSAVE)<-c("nbVilles","Nmeta","phiMAX","extGolLOWER","extGolRATE","extGolUPPER","perGolLOWER","perGolRATE","perGolUPPER")
	#save(resSAVE,file="kqua.Rda")
	return(resSAVE)
}
