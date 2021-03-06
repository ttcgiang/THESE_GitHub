#demo.globSEIRNewInfec
#We have two object typeSIMU that can be created. One with the direct algorithm of Gillespie. One other with the adaptivetau algorithm of Gillespie.
#We can give the initial values of variables for each ville by using the villeglobSEIRNewInfec0 parameter 
#	ville	S	E	I	R
#	0	1000	23	10	4566
#	..........
#If ville is't in villeglobSEIRNewInfec0, then this ville will use the values of S,E,I,R by defaults
demo.globSEIRNewInfec<-function(){
	#defaults
	directObj <- globSEIRNewInfec(N=1e7,method="direct")
	plot(directObj)
	pause()
		#changing the values of parameters or variables
		#equi parameter is "TRUE"
	directObj <- globSEIRNewInfec(N=1e7,typeSIMU="sto",method="direct",duration=5*365,nbVilles=3,equi=TRUE)
	plot(directObj)
	pause()
		#equi parameter is "FALSE"
	directObj <- globSEIRNewInfec(N=1e7,typeSIMU="stoch",equi=FALSE,typeRNG="good")
	plot(directObj)
	pause()
		#typeRNG parameter is "good" with equi=TRUE
	directObj <- globSEIRNewInfec(N=1e7,typeSIMU="stoch",nbVilles=3, typeRNG="good")
	plot(directObj)
	pause()
		#typeRNG parameter is "fast" with equi=TRUE
	directObj <- globSEIRNewInfec(N=1e7,typeSIMU="stoch",nbVilles=3, typeRNG="fast")
	plot(directObj)
	pause()
	
	#globSEIRNewInfec.adaptivetau object pour n population of a metapopulation. 
		#1ville	
	#adaptObj <- globSEIRNewInfec(method="adaptivetau")
	#plot(adaptObj)
	#pause()	
		#changing the values of parameters or variables
		#3villes
	#adaptObj <- globSEIRNewInfec(N=1e7,typeSIMU="stoch",method="adaptivetau",duration=10*365,nbVilles=3,equi=TRUE)
	#plot(adaptObj)
	#pause()
	#adaptObj <- globSEIRNewInfec(N=1e7,typeSIMU="stoch",method="adaptivetau",duration=10*365,nbVilles=3,equi=FALSE)
	#plot(adaptObj)
	#pause()
}
demo.globSEIRNewInfec()


#comparing the rand() function of C++ and the ran1() function of Pr.Yann Paris 13 
comp.typeRNGGood.typeRNGFast<-function(){
	nbRep = 10	
	good <- replicate(nbRep,system.time(globSEIRNewInfec(N=1e7,typeSIMU="stoch",duration=20*365,typeRNG="good"))[3])	
 	fast <- replicate(nbRep,system.time(globSEIRNewInfec(N=1e7,typeSIMU="stoch",duration=20*365,typeRNG="fast"))[3])
	plot(c(1:nbRep),good,ylim=range(0,good,fast),lwd=5,col="blue", ylab="time of simulation", xlab="number of repeations")
	points(c(1:nbRep),fast,lwd=5,col="black")
	title(main="the run-time of the good/fast random function",cex=0.6)
	box()
	legend(nbRep/2,max(good,fast)/2, c("good","fast"),pch=19,cex=0.5,col=c("blue","black"),lty=c(1,1))
	
	
}
comp.typeRNGGood.typeRNGFast()



#testing the generic functions
demo.print<- function(){
	#print
	onecity <- globSEIRNewInfec(N=1e7,typeSIMU="stoch")
	threecities <- globSEIRNewInfec(N=1e7,typeSIMU="stoch",nbVilles=3)
	#print(onecity)
	pop(onecity)
	pop(threecities)
	pause()
	#coef
	coef(onecity)
	coef(threecities)
	pause()
	summary(onecity)	
	pause()
	summary(threecities)
	pause()
}
demo.print()

#testing plot
#add parameter
#2d and 3d for plotting
demo.plot <- function(){
	#2d with add=F
	plot(globSEIRNewInfec(N=1e7,typeSIMU="stoch",nbVilles=2),add=F,ylab="P")
		pause()
	plot(globSEIRNewInfec(N=1e7,typeSIMU="stoch",nbVilles=1),add=T,ylab="P",colPops=c("black"))
		pause()
	#3d
	plot(globSEIRNewInfec(N=1e7,typeSIMU="stoch",nbVilles=1,duration=10*365),zlab="S",add=F,ylab="P",colPops=c("green"))
		pause()
	plot(globSEIRNewInfec(N=1e7,typeSIMU="stoch",nbVilles=3,duration=10*365),zlab="S",add=F,ylab="P",colPops=c("green"))
		pause()
}
demo.plot()


#testing for persistence
demo.stoch.persistence<-function(){
	#curvetypeSIMU="KM"
	p<-perNewInfec(globSEIRNewInfec(typeSIMU="sto",nbVilles=15,N=1e5))	
	plot.pers(p)
	pause()
	plot.comp.surv.estim(p)
	pause()
	#curvetypeSIMU="pops"
	plot.pers(p,curvetypeSIMU="pop",col=c("green","blue"),vilabline=c(1,3))
}
demo.stoch.persistence()



