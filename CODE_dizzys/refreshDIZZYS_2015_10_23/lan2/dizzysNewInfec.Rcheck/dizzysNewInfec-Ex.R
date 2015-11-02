pkgname <- "dizzysNewInfec"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('dizzysNewInfec')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("coef.seir")
### * coef.seir

flush(stderr()); flush(stdout())

### Name: coef.seir
### Title: Coefficient of seir Object
### Aliases: coef.seir
### Keywords: coefficients seir object

### ** Examples

	seirobj1<-seir(N=1e7)
	coef(seirobj1)

	seirobj2<-seir(nbVilles=3, N=c(1e7,1e6))
	coef(seirobj2)



cleanEx()
nameEx("confint.pers.rate")
### * confint.pers.rate

flush(stderr()); flush(stdout())

### Name: confint.pers.rate
### Title: Confidence interval for the estimated global disease persistence
###   rate in a metapopulation.
### Aliases: confint.pers.rate
### Keywords: seir model R package

### ** Examples

	#STO, STO
	sto<-seir(N=10e5,type="stoch",duration=30*365,beta1=0.1,nbVilles=10)
	objper<- persistence(sto)
	pers.obj<-pers.rate.obj(objper)
	confint.pers.rate(pers.obj)
	pause()



cleanEx()
nameEx("dizzysNewInfec-package")
### * dizzysNewInfec-package

flush(stderr()); flush(stdout())

### Name: dizzysNewInfec-package
### Title: Simulating SEIR/SIR models
### Aliases: dizzysNewInfec-package dizzysNewInfec
### Keywords: package

### ** Examples

~~ simple examples of the most important functions ~~



cleanEx()
nameEx("equilibrium")
### * equilibrium

flush(stderr()); flush(stdout())

### Name: equilibrium
### Title: Finding Endemic Equilibrium of a seasonally-forced SEIR/SIR
###   model
### Aliases: equilibrium
### Keywords: SEIR/SIR model limit cycle equilibrium

### ** Examples

## The point on the limit cycle depends on the input phase value 'phi':
	equilibrium(mu=.000456, beta0=1000/365, beta1=.1, sigma=1/7, gamma=1/7, phi=0, T=365,duration=100,unitTIME=1)
	equilibrium(mu=.000456, beta0=1000/365, beta1=.1, sigma=1/7, gamma=1/7, phi=pi/2, T=365, duration=100,unitTIME=7)



cleanEx()
nameEx("lines.seir")
### * lines.seir

flush(stderr()); flush(stdout())

### Name: lines.seir
### Title: Add Connected Line Segments of an seir object to a Plot 2D/3D
### Aliases: lines.seir
### Keywords: lines projection on plane

### ** Examples

#creating a plot
#adding a line to the plot
	seir(nbVilles=2)->obj
	seir(nbVilles=1)->obj1
	#2D
	plot(obj,col="red")
	lines(obj1,col="blue",lwd=2)
	#3D
	plot(obj,z="S",col="red",proj=list(c("time","P")))
	lines(obj1,z="S",col="blue",proj=list(c("time","P")))



cleanEx()
nameEx("pers.rate.obj")
### * pers.rate.obj

flush(stderr()); flush(stdout())

### Name: pers.rate.obj
### Title: Calculating the global disease persistence rate in a
###   metapopulation.
### Aliases: pers.rate.obj
### Keywords: seir model R package

### ** Examples

	#STO, STO
	sto<-seir(N=10e5,type="stoch",duration=30*365,beta1=0.1,nbVilles=10)
	objper<- persistence(sto)
	pers.obj<-pers.rate.obj(objper)
		pause()



cleanEx()
nameEx("persistence")
### * persistence

flush(stderr()); flush(stdout())

### Name: persistence
### Title: Persistence in a Metapopulation
### Aliases: persistence
### Keywords: metapopulation persistence

### ** Examples

	obj1<-seir.stoch(nbVilles=10,N=1e5)
	objper<-persistence(obj1)
	objper@persistence



cleanEx()
nameEx("plot.pers")
### * plot.pers

flush(stderr()); flush(stdout())

### Name: plot.pers
### Title: Plotting Kaplan Meier Survival Curve
### Aliases: plot.pers
### Keywords: Kaplan–Meier curve Kaplan–Meier estimator

### ** Examples

	p<-persistence(seir(type="sto",nbVilles=15,N=1e5))	
	plot.pers(p)
	x11()
	plot.pers(p,curvetype="pop",col=c("green","blue"),vilabline=c(1,3))



cleanEx()
nameEx("plot.seir")
### * plot.seir

flush(stderr()); flush(stdout())

### Name: plot.seir
### Title: Plotting 2D/3D a seir Object
### Aliases: plot.seir
### Keywords: seir model

### ** Examples

	 obj<-seir(nbVilles=3, N=5e5)
	plot(obj,col=c("red","blue"),lwd=2,xlab="time (day)", ylab="number of infectives")
	plot(obj,z="S",col=c("red","blue"),lwd=2,xlab="time (day)", ylab="number of infectives",zlab="number of susceptible")
	plot(obj,z="S",col=c("red","blue"),lwd=2,proj=list(c("time","P"),c("time","S")),box=F,xlab="time (day)", ylab="number of infectives",zlab="number of susceptible")



cleanEx()
nameEx("pop.seir")
### * pop.seir

flush(stderr()); flush(stdout())

### Name: pop.seir
### Title: Extract Values of State Variables of each City according to
###   Time.
### Aliases: pop.seir

### ** Examples

	obj<-seir(nbVilles=3, N=1e7)
	tpobj<-pop(obj) 
	class(tpobj)
	tpobj<-pop(obj,fct="sum")
	class(tpobj)
	tpobj<-pop(obj,subset=c(1,2),fct="sum")	
	class(tpobj)



cleanEx()
nameEx("print.seir")
### * print.seir

flush(stderr()); flush(stdout())

### Name: print.seir
### Title: Printing seir Object
### Aliases: print.seir

### ** Examples

	obj<-seir(nbVilles=3, N=1e7)
	print(obj)



cleanEx()
nameEx("seir-class")
### * seir-class

flush(stderr()); flush(stdout())

### Name: seir-class
### Title: seir Class '"seir"'
### Aliases: seir-class
### Keywords: classes

### ** Examples

showClass("seir")



cleanEx()
nameEx("simul")
### * simul

flush(stderr()); flush(stdout())

### Name: simul
### Title: Redoing or Continuing a Simulation.
### Aliases: simul
### Keywords: seir model R package

### ** Examples

	#STO, STO
	sto<-seir(N=10e5,type="stoch",duration=5*365,beta1=0.1,nbVilles=2)
	plot(simul(sto,type="stoch",continue=TRUE,duration=5*365,beta1=0,phi=c(pi/2,0)))
		pause()
	#DET, DET
	det<-seir(N=10e5,type="deter",duration=50*365)
	plot(simul(det,type="deter",continue=TRUE,duration=50*365,beta1=0,phi=pi/2))
		pause()



cleanEx()
nameEx("str.seir")
### * str.seir

flush(stderr()); flush(stdout())

### Name: str.seir
### Title: Describe the Structure of a seir Object.
### Aliases: str.seir
### Keywords: seir model R object

### ** Examples

	obj<-seir()
	str(obj)



cleanEx()
nameEx("summary.seir")
### * summary.seir

flush(stderr()); flush(stdout())

### Name: summary.seir
### Title: Object Summaries
### Aliases: summary.seir
### Keywords: summary seir class

### ** Examples

	obj<-seir()
	summary(obj)



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
