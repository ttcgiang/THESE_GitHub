demo.simul.sto <- function(){
	#STO
	objSTO <- stoSEIRNewInfec(N=10e5,nbVilles=2)
	plot(objSTO,add=F,col="red")
		pause()
	plot(globSEIRSimulNewInfec(objSTO,duration=10*365,nbVilles=1),add=T,col=c("black"))
		pause()
	plot(globSEIRSimulNewInfec(objSTO,duration=10*365,nbVilles=3),add=F)
		pause()
	#STO, DET
	#objSTO <- stoSEIRNewInfec(nbVilles=2)
	plot(objSTO,add=F)
	plot(globSEIRSimulNewInfec(objSTO,typeSIMU="deter",duration=10*365),add=T,col="black",lwd=3)
	pause()

	objSTO <- stoSEIRNewInfec(N=10e5,nbVilles=1,duration=50*365)
	plot(objSTO,add=F)
	plot(globSEIRSimulNewInfec(objSTO,typeSIMU="deter",duration=50*365),add=T,col="black",lwd=3)
}
demo.simul.sto()


demo.simul.det<- function(){
	#DET
	objDET <- detSEIRNewInfec(N=10e5,duration=30*365)
	plot(objDET,add=F)
		pause()
	plot(globSEIRSimulNewInfec(objDET,duration=60*365,phiPHASE=pi/2),add=T,col=c("black"))
		pause()

	#DET, STO
	objDET <- detSEIRNewInfec(duration=20*365,S=741559,E=2794,I=1675,N=10e6)
	plot(objDET,add=F,col="red")
	plot(globSEIRSimulNewInfec(objDET,typeSIMU="stoch",typeRNG="good",nbVilles=1,duration=10*365),add=T,col="black",lwd=2)
	pause()
	plot(globSEIRSimulNewInfec(objDET,typeSIMU="stoch",typeRNG="fast",nbVilles=2,duration=10*365),add=F,col="black",lwd=1)
}
demo.simul.det()



demo.continue <- function(){
	#STO, STO
	sto<-stoSEIRNewInfec(N=10e5,duration=5*365,nbCONTACT1=0.1,nbVilles=2)
	plot(globSEIRSimulNewInfec(sto,typeSIMU="stoch",continue=T,duration=5*365,nbCONTACT1=0,phiPHASE=c(pi/2,0)))
		pause()
	#DET, DET
	det<-detSEIRNewInfec(N=10e5,duration=50*365)
	plot(globSEIRSimulNewInfec(det,typeSIMU="deter",continue=T,duration=50*365,nbCONTACT1=0,phiPHASE=pi/2))
		pause()

	#DET, STO
	plot(globSEIRSimulNewInfec(detSEIRNewInfec(duration=5*365,S=741559,E=2794,I=1675,N=10e6),typeRNG="fast",typeSIMU="stoch",continue=T,nbCONTACT1=0,duration=20*365))
		pause()
	plot(globSEIRSimulNewInfec(detSEIRNewInfec(duration=50*365,N=10e5),typeRNG="fast",continue=T,duration=10*365,typeSIMU="stoch"))
		pause()
	#STO, DET
	plot(globSEIRSimulNewInfec(stoSEIRNewInfec(N=10e5,duration=5*365,nbVilles=2),typeSIMU="deter",continue=T,duration=20*365))
		pause()
	plot(globSEIRSimulNewInfec(stoSEIRNewInfec(N=10e5,duration=5*365,nbVilles=2),typeSIMU="deter",continue=T,duration=10*365,nbCONTACT1=0,phiPHASE=pi/2))	
}
demo.continue()


