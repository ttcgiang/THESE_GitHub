lvrates <- function(x, params, t) {
nu<-params["nu"]
mu<-params["mu"]
beta<-params["beta"]
sigma<-params["sigma"]
gamma<-params["gamma"]
S<-x["S"]
E<-x["E"]
I<-x["I"]
R<-x["R"]

return(c(nu,mu*(S+E+I+R),
## prey growth rate
mu*S,mu*E,mu*I,mu*R,
## prey death / predator growth rate
beta*S*I/(S+E+I+R),sigma*E,gamma*I))
}

seir_stochAdaptivetau <- function(tmax=10*365,N=1e+05,S=7304,E=29,I=18,mu=3.913894e-05,
	beta=1000/365,sigma=0.125,gamma=0.2,nu=0.0,...) {
	library("adaptivetau")
	S <- round(S)
	E <- round(E)
	I <- round(I)
	R <- N-S-E-I
# S, E, and I are numbers of susceptibles, exposed and infectious
# respectively. All the rates are defined per day.
#       mu is the death (and birth) rate
#       beta is the contact rate
#       sigma is the inverse of the average duration of exposed
#       gamma is the recovery rate
# The model is defined by the following 3 differential equations:
	params <- c(mu=mu,beta=beta,sigma=sigma,gamma=gamma,nu=nu)
	x <- c(S=S,E=E,I=I,R=R)
# The model is defined by the following state-change matrix:
	model <- t(matrix(c(
#		 S   E   I   R
		-1,  0,  1,  0,  # infection from outside
		 1,  0,  0,  0,  # birth
		-1,  0,  0,  0,  # susceptible death
		 0, -1,  0,  0,  # exposed death
		 0,  0, -1,  0,  # infectious death
		 0,  0,  0, -1,  # recovered death
		-1,  1,  0,  0,  # infection
		 0, -1,  1,  0,  # becoming infectious
		 0,  0, -1,  1), # recovery
		ncol=4,byrow=T))		
# Propensity vector:
	a <- c(nu,mu*(S+E+I+R),mu*S,mu*E,mu*I,mu*R,
		beta*S*I/(S+E+I+R),sigma*E,gamma*I)
# Stochastic simulations by direct method (with GillespieSSA package):
	print(x)
	print(a)
	print(model)
	print(params)
	print("fkfkkffkfkfkfkfk")
	kq<-lvrates(x,params,0)
	print(kq)
	out <- ssa.exact(x,model,lvrates,params,tf=tmax)
	print("ffffffj")
#	print(out)
# Putting the output of stochastic simulations in good shape:
	#out <- out$data
	rownames(out) <- paste(1:nrow(out))
	out <- as.data.frame(out)
	names(out) <- c("time","S","E","I","R")
# Ploting the prevalence as a function of time:
	#plot(out$time/365,out$I,type="l",col="red",
		#xlab="time (years)",ylab="prevalence",...)
# Optionally giving the data frame as an output:
	#invisible(out)
	return(out)
}

figureAdaptivetau<-function(nb=10){
	for(i in 1:nb) {
		out<- seir_stochAdaptivetau()
	if(i==1) {
		plot(out$time/365,out$I,type="l",col="red",xlab="time (years)",ylab="prevalence per year",ylim=range(1,200),main="SEIR model with package adapivetau")
	}
	else{
		lines(out$time/365,out$I,type="l",col="red",xlab="time (years)",ylab="prevalence")
	}
	}
}

