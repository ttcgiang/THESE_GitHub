pkgname <- "dizzysNewInfec"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('dizzysNewInfec')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("coef.seirNewInfec")
### * coef.seirNewInfec

flush(stderr()); flush(stdout())

### Name: coef.seirNewInfec
### Title: Coefficient of seirNewInfec Object
### Aliases: coef
### Keywords: coefficients seirNewInfec object

### ** Examples

	seirobj1<-globSEIRSimulNewInfec(N=1e7)
	coef(seirobj1)

	seirobj2<-globSEIRSimulNewInfec(nbVilles=3, N=c(1e7,1e6))
	coef(seirobj2)



cleanEx()
nameEx("detSEIRNewInfec")
### * detSEIRNewInfec

flush(stderr()); flush(stdout())

### Name: detSEIRNewInfec
### Title: Creat a seir Object
### Aliases: detSEIRNewInfec
### Keywords: deterministic model

### ** Examples

	obj<-globSEIRNewInfec(typeSIMU="deter",nbVilles=3,N=c(1e7,1e6))
	plot(obj)	



cleanEx()
nameEx("equiNewInfec")
### * equiNewInfec

flush(stderr()); flush(stdout())

### Name: equiNewInfec
### Title: Finding Endemic Equilibrium of a seasonally-forced SEIR/SIR
###   model
### Aliases: equiNewInfec
### Keywords: SEIR/SIR model limit cycle equilibrium

### ** Examples

## The point on the limit cycle depends on the input phase value 'phi':
	res<-equiNewInfec(duration=100*365,unitTIME=1,N=10e6,mu=1/(70*365),
    nbCONTACT0=300,nbCONTACT1=.1,probINFECTER=0.1,sigma=1/7,gamma=1/7,phiPHASE=c(0),periDISE=365)	
	print(res)



cleanEx()
nameEx("globSEIRNewInfec")
### * globSEIRNewInfec

flush(stderr()); flush(stdout())

### Name: globSEIRNewInfec
### Title: Creat a seir Object
### Aliases: globSEIRNewInfec
### Keywords: deterministic model stochastic model

### ** Examples

	obj<-globSEIRNewInfec()
	plot(obj)
	obj<-globSEIRNewInfec(typeSIMU="sto",nbVilles=3,N=c(1e7,1e6))
	plot(obj)
	obj<-globSEIRNewInfec(typeSIMU="deter",nbVilles=3,N=c(1e7,1e6))
	plot(obj)	



cleanEx()
nameEx("globSEIRSimulNewInfec")
### * globSEIRSimulNewInfec

flush(stderr()); flush(stdout())

### Name: globSEIRSimulNewInfec
### Title: Redoing or Continuing a Simulation.
### Aliases: globSEIRSimulNewInfec
### Keywords: seir model R package

### ** Examples

	#STO, STO
	sto<-globSEIRNewInfec(N=1e6,typeSIMU="stoch",duration=5*365,nbVilles=2)
	plot(globSEIRSimulNewInfec(sto,typeSIMU="stoch",continue=TRUE,duration=5*365,nbCONTACT1=0,phiPHASE=c(pi/2,0)),col=c(1,2))
	#DET, DET
	det<-globSEIRNewInfec(N=10e4,typeSIMU="deter",duration=50*365)
	plot(globSEIRSimulNewInfec(det,typeSIMU="deter",continue=TRUE,duration=5*365,nbCONTACT1=0.1,phiPHASE=pi))



cleanEx()
nameEx("lines.seirNewInfec")
### * lines.seirNewInfec

flush(stderr()); flush(stdout())

### Name: lines.seirNewInfec
### Title: Add Connected Line Segments of an seir object to a Plot 2D/3D
### Aliases: lines
### Keywords: lines projection on plane

### ** Examples

#creating a plot
#adding a line to the plot
	globSEIRSimulNewInfec(nbVilles=2)->obj
	globSEIRSimulNewInfec(nbVilles=1)->obj1
	#2D
	plot(obj,col="red")
	lines(obj1,col="blue",lwd=2)
	#3D
	plot(obj,z="S",col="red",proj=list(c("time","P")))
	lines(obj1,z="S",col="blue",proj=list(c("time","P")))



cleanEx()
nameEx("persNewInfec")
### * persNewInfec

flush(stderr()); flush(stdout())

### Name: persNewInfec
### Title: persNewInfec in a Metapopulation
### Aliases: persNewInfec
### Keywords: metapopulation persNewInfec

### ** Examples

	obj1<-obj1<-globSEIRSimulNewInfec(nbVilles=5,N=1e5,nbCONTACT0=100,duration=365*30)
	objper<-persNewInfec(obj1)
	objper@persistence



cleanEx()
nameEx("plot.persNewInfec")
### * plot.persNewInfec

flush(stderr()); flush(stdout())

### Name: plot.persNewInfec
### Title: Plotting Kaplan Meier Survival Curve
### Aliases: plot
### Keywords: Kaplan–Meier curve Kaplan–Meier estimator

### ** Examples

	p<-persNewInfec(globSEIRSimulNewInfec(nbVilles=5,N=1e5,nbCONTACT0=100,duration=365*30))
	plot.persNewInfec(p)
	x11()
	plot.persNewInfec(p,curvetype="pop",col=c("green","blue"),vilabline=c(1,3))



cleanEx()
nameEx("plot.seirNewInfec")
### * plot.seirNewInfec

flush(stderr()); flush(stdout())

### Name: plot.seirNewInfec
### Title: Plotting 2D/3D a seir Object
### Aliases: plot
### Keywords: seir model

### ** Examples

	obj<-globSEIRSimulNewInfec(nbVilles=3, N=1e6, nbCONTACT0=100)
	plot(obj,col=c("red","blue"),lwd=2,xlab="time (day)", ylab="number of infectives")
	pause()
	plot(obj,z="S",col=c("red","blue"),lwd=2,xlab="time (day)", ylab="number of infectives",zlab="number of susceptible")
	pause()
	#plot(obj,z="S",col=c("red","blue"),lwd=2,proj=list(c("time","P"),c("time","S")),box=F,xlab="time (day)", ylab="number of infectives",zlab="number of susceptible")



cleanEx()
nameEx("pop.seirNewInfec")
### * pop.seirNewInfec

flush(stderr()); flush(stdout())

### Name: pop.seirNewInfec
### Title: Extract Values of State Variables of each City according to
###   Time.
### Aliases: pop.seirNewInfec

### ** Examples

	obj<-globSEIRSimulNewInfec(nbVilles=3, N=1e6)
	tpobj<-pop(obj) 
	class(tpobj)
	tpobj<-pop(obj,fct="sum")
	class(tpobj)
	tpobj<-pop(obj,subset=c(1,2),fct="sum")	
	class(tpobj)



cleanEx()
nameEx("stoSEIRNewInfec")
### * stoSEIRNewInfec

flush(stderr()); flush(stdout())

### Name: stoSEIRNewInfec
### Title: Creat a seir stochastic Object
### Aliases: stoSEIRNewInfec
### Keywords: stochastic model

### ** Examples

	obj<-globSEIRNewInfec(typeSIMU="sto",nbVilles=3,N=c(1e7,1e6))
	plot(obj)	



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
