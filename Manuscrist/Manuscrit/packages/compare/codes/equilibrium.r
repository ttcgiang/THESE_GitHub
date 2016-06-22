equilibrium <- function(duration=100,mu=1/(70*365),
    beta0=1000/365,beta1=.1,sigma=1/7,gamma=1/7,phi=0,T=365) {
# Finding the equilibrium on a long duration and a small population size,
# starting with random initial values for the variables:
	Seq <- (gamma+mu)*(sigma+mu)/(beta0 * sigma)
	Eeq <- mu*((1/(sigma+mu)) - ((gamma+mu)/(beta0 * sigma)))
	Ieq <- mu*((beta0*sigma - (gamma+mu)*(sigma+mu))/(beta0 *(gamma+mu)*(sigma+mu)))
	Req <- 1- Seq - Eeq - Ieq
	output <- seir_det(duration=duration,N=1,S=Seq,E=Eeq,I=Ieq,T=T,mu=mu,beta0=beta0,
		beta1=beta1,sigma=sigma,gamma=gamma,phi=phi)
	equal <- tail(output,1)[-1]
        vecEqual <- rbind(equal[1,1],equal[1,2],equal[1,3])
	vecEqual <- as.vector(vecEqual)
	return(vecEqual)
}

