demo.sir <- function(){
	#plotting deterministic SIR model with sigma=Inf and other correspondent parameters
	# refer to the book "THE STOCHASTIC DYNAMICS OF EPIDEMIC MODELS"  of Andrew James Black in 2010	
	plot(detSEIRNewInfec(duration=100*365,sigma=Inf,N=10e6))
	pause()

	#Comparing deterministic SIR model with stochastic SIR model
	S=700000; E=0; I=100; R=9299900;

	plot(globSEIRNewInfec(duration=5*365,type="stoch",S=S,E=E,I=I,R=R,nbCONTACT0=150, nbCONTACT1=0.0,sigma=Inf,gamma=0.077),col="red")
	plot(globSEIRNewInfec(duration=5*365,type="deter",S=S,E=E,I=I,R=R,nbCONTACT0=150, nbCONTACT1=0.00,sigma=Inf,gamma=0.077),col="blue",add=T,lty=3,lwd=2)
	pause()	

	plot(globSEIRNewInfec(duration=5*365,type="stoch",nbVilles=2,S=S,E=E,I=I,R=R,nbCONTACT0=150, nbCONTACT1=0.0,sigma=Inf,gamma=0.077),col="red")
	plot(globSEIRNewInfec(duration=5*365,type="deter",S=S,E=E,I=I,R=R,nbCONTACT0=150, nbCONTACT1=0.00,sigma=Inf,gamma=0.077),col="blue",add=T,lty=3,lwd=2)
	pause()	

	#Comparing SIR model direct/adaptivetau methods
		#with "direct"/"adaptivetau" method
	plot(globSEIRNewInfec(duration=5*365,type="stoch",method="direct",nbVilles=1,S=S,E=E,I=I,R=R,nbCONTACT0=150, nbCONTACT1=0.0,sigma=Inf,gamma=0.077),col="red")
	plot(globSEIRNewInfec(duration=5*365,type="stoch",method="adaptivetau",nbVilles=1,S=S,E=E,I=I,R=R,nbCONTACT0=150, nbCONTACT1=0.0,sigma=Inf,gamma=0.077),add=T,col="blue")
}
demo.sir()

