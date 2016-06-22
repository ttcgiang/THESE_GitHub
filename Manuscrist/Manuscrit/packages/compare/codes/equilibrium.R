equilibrium <-
function(mu=1/(70*365),
    beta0=1000/365,beta1=.1,sigma=1/7,gamma=1/7,phi=0,T=365) {
# Finding the equilibrium on a long duration and a small population size,
# starting with random initial values for the variables:
	Seq <- (gamma+mu)*(sigma+mu)/(beta0 * sigma)
	Eeq <- mu*((1/(sigma+mu)) - ((gamma+mu)/(beta0 * sigma)))
	Ieq <- mu*((beta0*sigma - (gamma+mu)*(sigma+mu))/(beta0 *(gamma+mu)*(sigma+mu)))
	Req <- 1- Seq - Eeq - Ieq
	output <- seir_det(times=seq(0,100*365,le=5000),N=1,Seq,Eeq,Ieq,T,mu,beta0,
		beta1,sigma,gamma,phi,F)
	return(tail(output,1)[-1])
}
