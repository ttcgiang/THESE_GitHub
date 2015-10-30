demo.simul.sto <- function(){
	#STO
	objSTO <- seir.stoch(N=10e5,nbVilles=2)
	plot(objSTO,add=F,col="red")
		pause()
	plot(simul(objSTO,duration=10*365,nbVilles=1),add=T,col=c("black"))
		pause()
	plot(simul(objSTO,duration=10*365,nbVilles=3),add=F)
		pause()
	#STO, DET
	#objSTO <- seir.stoch(nbVilles=2)
	plot(objSTO,add=F)
	plot(simul(objSTO,type="deter",duration=10*365),add=T,col="black",lwd=3)
	pause()

	objSTO <- seir.stoch(N=10e5,nbVilles=1,duration=50*365)
	plot(objSTO,add=F)
	plot(simul(objSTO,type="deter",duration=50*365),add=T,col="black",lwd=3)
}
demo.simul.sto()


demo.simul.det<- function(){
	#DET
	objDET <- seir.det(N=10e5,duration=30*365)
	plot(objDET,add=F)
		pause()
	plot(simul(objDET,duration=60*365,phi=pi/2),add=T,col=c("black"))
		pause()

	#DET, STO
	objDET <- seir.det(duration=20*365,S=741559,E=2794,I=1675,N=10e6)
	plot(objDET,add=F,col="red")
	plot(simul(objDET,type="stoch",rng="good",nbVilles=1,duration=10*365),add=T,col="black",lwd=2)
	pause()
	plot(simul(objDET,type="stoch",rng="fast",nbVilles=2,duration=10*365),add=F,col="black",lwd=1)
}
demo.simul.det()



demo.continue <- function(){
	#STO, STO
	sto<-seir.stoch(N=10e5,duration=5*365,beta1=0.1,nbVilles=2)
	plot(simul(sto,type="stoch",continue=T,duration=5*365,beta1=0,phi=c(pi/2,0)))
		pause()
	#DET, DET
	det<-seir.det(N=10e5,duration=50*365)
	plot(simul(det,type="deter",continue=T,duration=50*365,beta1=0,phi=pi/2))
		pause()

	#DET, STO
	plot(simul(seir.det(duration=5*365,S=741559,E=2794,I=1675,N=10e6),rng="fast",type="stoch",continue=T,beta1=0,duration=20*365))
		pause()
	plot(simul(seir.det(duration=50*365,N=10e5),rng="fast",continue=T,duration=10*365,type="stoch"))
		pause()
	#STO, DET
	plot(simul(seir.stoch(N=10e5,duration=5*365,nbVilles=2),type="deter",continue=T,duration=20*365))
		pause()
	plot(simul(seir.stoch(N=10e5,duration=5*365,nbVilles=2),type="deter",continue=T,duration=10*365,beta1=0,phi=pi/2))	
}
demo.continue()


