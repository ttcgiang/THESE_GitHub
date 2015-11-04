demo.globaldizzysNewInfec<-function(){
	globSEIRSimulNewInfec(N=1e7,nbVilles=3,typeSIMU="stoc")->a
	plot(a,col=c(1,2,5))
	pause()
	plot(a,z="S",col=c(1,2,5))
	pause()
	plot(a,z="S",col=c(1,2,5),proj=list(c("time","P")),box=T)
	pause()
	plot(a,z="S",col="red")
	pause()

	globSEIRSimulNewInfec(N=1e7)->a
	plot(a)
	pause()
	b<-globSEIRSimulNewInfec(a,continue=T,append=T,t0=1000)
	plot(b)
	pause()
	plot(a,add=T,col="red")
	pause()

	b<-globSEIRSimulNewInfec(a,continue=T,append=T,t0=1000,duration=20*365)
	plot(b)
	pause()

	b<-globSEIRSimulNewInfec(a,continue=T,append=T,t0=NULL)
	plot(b)
	pause()

	b<-globSEIRSimulNewInfec(a,typeSIMU="sto",continue=T,append=F,t0=1000)
	plot(b)
	pause()
	b<-globSEIRSimulNewInfec(a,typeSIMU="sto",nbVilles=4,continue=T,append=T,t0=1000)
	plot(b,col=c(1,2)) 
	pause()
	b<-globSEIRSimulNewInfec(a,typeSIMU="sto",nbVilles=4,continue=T,append=F,t0=800)
	plot(b)
	pause()

	globSEIRSimulNewInfec(typeSIMU="sto",N=1e7)->a
	plot(a)
	pause()
	globSEIRSimulNewInfec(a,typeSIMU="det",continue=T,append=F)->b
	plot(b)
	pause()
	globSEIRSimulNewInfec(typeSIMU="sto",nbVilles=3,N=c(1e7,1e6))->a
	globSEIRSimulNewInfec(a,typeSIMU="det",continue=T,append=T)->b
	plot(b)
	pause()
	summary(b)
	str(b)
	coef(b)
}
demo.globaldizzysNewInfec()
