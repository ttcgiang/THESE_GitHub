demo.dizzysgene<-function(){
	simul(N=1e7,nbVilles=3,type="stoc")->a
	plot(a,col=c(1,2,5))
	pause()
	plot(a,z="S",col=c(1,2,5))
	pause()
	plot(a,z="S",col=c(1,2,5),proj=list(c("time","P")),box=T)
	pause()
	plot(a,z="S",col="red")
	pause()

	simul(N=1e7)->a
	plot(a)
	pause()
	b<-simul(a,continue=T,append=T,t0=1000)
	plot(b)
	pause()
	plot(a,add=T,col="red")
	pause()

	b<-simul(a,continue=T,append=T,t0=1000,duration=20*365)
	plot(b)
	pause()

	b<-simul(a,continue=T,append=T,t0=NULL)
	plot(b)
	pause()

	b<-simul(a,type="sto",continue=T,append=F,t0=1000)
	plot(b)
	pause()
	b<-simul(a,type="sto",nbVilles=4,continue=T,append=T,t0=1000)
	plot(b,col=c(1,2)) 
	pause()
	b<-simul(a,type="sto",nbVilles=4,continue=T,append=F,t0=800)
	plot(b)
	pause()

	simul(type="sto",N=1e7)->a
	plot(a)
	pause()
	simul(a,type="det",continue=T,append=F)->b
	plot(b)
	pause()
	simul(type="sto",nbVilles=3,N=c(1e7,1e6))->a
	simul(a,type="det",continue=T,append=T)->b
	plot(b)
	pause()
	summary(b)
	str(b)
	coef(b)
}
demo.dizzysgene()
