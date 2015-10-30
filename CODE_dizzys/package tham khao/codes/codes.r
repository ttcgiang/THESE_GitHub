equilibrium <- function(duration=10,N=1e6,mu=1/(70*365),
    beta0=1000/365,beta1=.1,sigma=1/7,gamma=1/7,phi=0,T=365) {
# Finding the equilibrium on a long duration and a small population size,
# starting with random initial values for the variables:
	Nsmall <- 1000
	S <- .1
	E <- 0
	I <- .0001
	output <- seir_det(times=seq(0,100*365,le=5000),Nsmall,S,E,I,T,mu,beta0,
		beta1,sigma,gamma,phi,F)
	return(tail(output,1)[-1]/Nsmall)
}



deterministic_dynamics <- function(duration=10,N=1e6,S,E,I,mu=1/(70*365),
    beta0=1000/365,beta1=.1,sigma=1/7,gamma=1/7,phi=0,T=365) {
# Finding the equilibrium on a long duration and a small population size,
# starting with random initial values for the variables:
	equil <- equilibrium(duration,N,mu,beta0,beta1,sigma,gamma,phi,T)
# Finding the maximum I on this equilibrium
#	S <- as.numeric(equil["S"])
#	E <- as.numeric(equil["E"])
#	I <- as.numeric(equil["I"])
	output <- seir_det(times=seq(0,duration*365,le=50*duration),N,S,E,
		I,T,mu,beta0,beta1,sigma,gamma,phi,F)
	return(output)
}



comparison_det_stoch0 <- function(duration=10,N=1e6,mu=1/(70*365),
    beta0=1000/365,beta1=.1,sigma=1/7,gamma=1/7,phi=0,T=365) {
	equil <- equilibrium(duration,N,mu,beta0,beta1,sigma,gamma,phi,T)
	S <- as.numeric(equil["S"])
	E <- as.numeric(equil["E"])
	I <- as.numeric(equil["I"])
	det <- deterministic_dynamics(duration,N,S,E,I,mu,beta0,beta1,sigma,gamma,phi,T)
	S0 <- round(S*N)
	E0 <- round(E*N)
	I0 <- round(I*N)
	R0 <- N - S0 - E0 - I0
	system("rm Output/*")
	system(paste("code_ttcgiang_12_2012/Projet_SEIR_ver01 -nbVilles 1 -tmax",duration*365,"-sigma",sigma,"-gamma",gamma,"-mu",mu,"-epsilon 0 -topology 0 -rho 0 -unitTemps 1 -graine",as.numeric(Sys.time()),"-S0",S0,"-E0",E0,"-I0",I0,"-R0",R0,"-N0",N,"-beta0",beta0,"-beta1",beta1,"-phiMIN 0 -phiMAX 0 -nbSimu",1,"-typeFonc 0"))
	system("mv Output/Vil1* Output/Vil1.csv")
	simulations <- read.table("Output/Vil1.csv",header=T)
	with(det,plot(time/365,I,type="l",lwd=1,col="red",ylim=range(c(det$I,simulations$I)),xlab="time (year)",ylab="number of infectives"))
	with(simulations,lines(t.jours./365,I))
	with(det,lines(time/365,I,lwd=1,col="red"))
}




comparison_det_stoch <- function(duration=10,N=1e6,mu=1/(70*365),
    beta0=1000/365,beta1=.1,sigma=1/7,gamma=1/7,phi=0,T=365) {
	equil <- equilibrium(duration,N,mu,beta0,beta1,sigma,gamma,phi,T)
	S <- as.numeric(equil["S"])
	E <- as.numeric(equil["E"])
	I <- as.numeric(equil["I"])
	det <- deterministic_dynamics(duration,N,S,E,I,mu,beta0,beta1,sigma,gamma,phi,T)
	S0 <- round(S*N)
	E0 <- round(E*N)
	I0 <- round(I*N)
	R0 <- N - S0 - E0 - I0
	system("rm Output/*")
	system(paste("code_ttcgiang_12_2012/Projet_SEIR_ver01 -nbVilles 2 -tmax",duration*365,"-sigma",sigma,"-gamma",gamma,"-mu",mu,"-epsilon 0 -topology 0 -rho 0 -unitTemps 1 -graine",as.numeric(Sys.time()),"-S0",S0,"-E0",E0,"-I0",I0,"-R0",R0,"-N0",N,"-beta0",beta0,"-beta1",beta1,"-phiMIN 0 -phiMAX 0 -nbSimu",1,"-typeFonc 0"))
	system("mv Output/Vil1* Output/Vil1.csv")
	system("mv Output/Vil2* Output/Vil2.csv")
	simulations1 <- read.table("Output/Vil1.csv",header=T)
	simulations2 <- read.table("Output/Vil2.csv",header=T)
	with(det,plot(time/365,I,type="l",lwd=1,col="red",ylim=range(c(det$I,simulations1$I)),xlab="time (year)",ylab="number of infectives"))
	with(simulations1,lines(t.jours./365,I))
	with(simulations2,lines(t.jours./365,I,col="grey"))
	with(det,lines(time/365,I,lwd=1,col="red"))
}
