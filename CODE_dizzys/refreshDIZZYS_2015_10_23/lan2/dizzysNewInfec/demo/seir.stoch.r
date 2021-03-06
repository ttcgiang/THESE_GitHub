#demo.seir
#We have two object type that can be created. One with the direct algorithm of Gillespie. One other with the adaptivetau algorithm of Gillespie.
#We can give the initial values of variables for each ville by using the villeSEIR0 parameter 
#	ville	S	E	I	R
#	0	1000	23	10	4566
#	..........
#If ville is't in villeSEIR0, then this ville will use the values of S,E,I,R by defaults
demo.seir<-function(){
	#defaults
	directObj <- seir(N=1e7,method="direct")
	plot(directObj)
	pause()
		#changing the values of parameters or variables
		#equi parameter is "TRUE"
	directObj <- seir(N=1e7,type="sto",method="direct",duration=5*365,nbVilles=3,equi=TRUE)
	plot(directObj)
	pause()
		#equi parameter is "FALSE"
	directObj <- seir(N=1e7,type="stoch",equi=FALSE,rng="good")
	plot(directObj)
	pause()
		#rng parameter is "good" with equi=TRUE
	directObj <- seir(N=1e7,type="stoch",nbVilles=3, rng="good")
	plot(directObj)
	pause()
		#rng parameter is "fast" with equi=TRUE
	directObj <- seir(N=1e7,type="stoch",nbVilles=3, rng="fast")
	plot(directObj)
	pause()
	
	#seir.adaptivetau object pour n population of a metapopulation. 
		#1ville	
	adaptObj <- seir(method="adaptivetau")
	plot(adaptObj)
	pause()	
		#changing the values of parameters or variables
		#3villes
	adaptObj <- seir(N=1e7,type="stoch",method="adaptivetau",duration=10*365,nbVilles=3,equi=TRUE)
	plot(adaptObj)
	pause()
	adaptObj <- seir(N=1e7,type="stoch",method="adaptivetau",duration=10*365,nbVilles=3,equi=FALSE)
	plot(adaptObj)
	pause()
}
demo.seir()


#comparing the rand() function of C++ and the ran1() function of Pr.Yann Paris 13 
comp.rngGood.rngFast<-function(){
	nbRep = 10	
	good <- replicate(nbRep,system.time(seir(N=1e7,type="stoch",duration=20*365,rng="good"))[3])	
 	fast <- replicate(nbRep,system.time(seir(N=1e7,type="stoch",duration=20*365,rng="fast"))[3])
	plot(c(1:nbRep),good,ylim=range(0,good,fast),lwd=5,col="blue", ylab="time of simulation", xlab="number of repeations")
	points(c(1:nbRep),fast,lwd=5,col="black")
	title(main="the run-time of the good/fast random function",cex=0.6)
	box()
	legend(nbRep/2,max(good,fast)/2, c("good","fast"),pch=19,cex=0.5,col=c("blue","black"),lty=c(1,1))
	
	
}
comp.rngGood.rngFast()



#testing the generic functions
demo.print<- function(){
	#print
	onecity <- seir(N=1e7,type="stoch")
	threecities <- seir(N=1e7,type="stoch",nbVilles=3)
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
	plot(seir(N=1e7,type="stoch",nbVilles=2),add=F,ylab="P")
		pause()
	plot(seir(N=1e7,type="stoch",nbVilles=1),add=T,ylab="P",colPops=c("black"))
		pause()
	#3d
	plot(seir(N=1e7,type="stoch",nbVilles=1,duration=10*365),zlab="S",add=F,ylab="P",colPops=c("green"))
		pause()
	plot(seir(N=1e7,type="stoch",nbVilles=3,duration=10*365),zlab="S",add=F,ylab="P",colPops=c("green"))
		pause()
}
demo.plot()


#testing for persistence
demo.stoch.persistence<-function(){
	#curvetype="KM"
	p<-persistence(seir(type="sto",nbVilles=15,N=1e5))	
	plot.pers(p)
	pause()
	plot.comp.surv.estim(p)
	pause()
	#curvetype="pops"
	plot.pers(p,curvetype="pop",col=c("green","blue"),vilabline=c(1,3))
}
demo.stoch.persistence()



