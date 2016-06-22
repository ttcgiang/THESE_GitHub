equilibrium <- function(duration=100,unitTIME=1,N=10e6,mu=1/(70*365),
    beta0=1000/365,beta1=.1,sigma=1/7,gamma=1/7,phi=0,T=365) {
# Finding the equilibrium on a long duration and a small population size,
# starting with random initial values for the variables:
	sei_eq <- .Call("getEquilibrium",as.numeric(mu),as.numeric(beta0),as.numeric(sigma),as.numeric(gamma))	
	sei_eq <- sei_eq*N
	obj_det <- seir.det(duration=duration,unitTIME=unitTIME,N=N,S=sei_eq[1],E=sei_eq[2],I=sei_eq[3],T=T,mu=mu,beta0=beta0,
		beta1=beta1,sigma=sigma,gamma=gamma,phi=phi)
	#plot(obj_det)
	output <- obj_det@pop[[1]]
#returning a vector of the equilibrium values of S, E, I
	return(unlist(tail(output,1)[-1]))
}

