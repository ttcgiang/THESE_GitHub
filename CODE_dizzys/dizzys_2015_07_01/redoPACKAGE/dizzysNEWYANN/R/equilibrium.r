#
# Created by TRAN Thi Cam Giang
# Function permit us to transmit the values of parameters,
# Then, call the function C++ in R, to caculating the equilibrium variables of the SEIR model
#
equilibrium <- function(duration=100,unitTIME=1,N=10e6,mu=1/(70*365),
    beta0=1000/365,beta1=.1,sigma=1/7,gamma=1/7,phi=0,T=365) {

	#Call the function C++ in R and at the same time, transmit the values of parameters
	# return S*, E*, I*
	sei_eq <- .Call("getEquilibrium",as.numeric(mu),as.numeric(beta0),as.numeric(sigma),as.numeric(gamma))	
	#S*, E*, I* (*N)
	sei_eq <- sei_eq*N
	#Doing one deterministic simulation in 100 years
	obj_det <- seir.det(duration=duration,unitTIME=unitTIME,N=N,S=sei_eq[1],E=sei_eq[2],I=sei_eq[3],T=T,mu=mu,beta0=beta0, beta1=beta1,sigma=sigma,gamma=gamma,phi=phi)
	
	#getting the values of variables of the population 1
	output <- obj_det@pop[[1]]
        #returning a vector of the equilibrium values of S, E, I
	return(unlist(tail(output,1)[-1]))
}

#################The end #######################

