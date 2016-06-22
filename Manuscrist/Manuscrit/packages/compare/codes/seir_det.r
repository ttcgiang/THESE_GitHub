seir_det <- function(times=seq(0,100*365,le=5000),N=10000,S=.1,E=0,I=.0001,T=365,
	mu=1/(70*365),beta0=1000/365,beta1=.1,sigma=1/7,gamma=1/7,phi=0,plot=T,...) {
	require(deSolve)
# N is the (constant) total population size, S, E, and I are numbers of
# susceptibles, exposed and infectious respectively. All the rates are
# defined per day.
# 	mu is the death (and birth) rate
# 	beta is the contact rate
# 	sigma is the inverse of the average duration of exposed
# 	gamma is the recovery rate
# The model is defined by the following 3 differential equations:
	model <- function(t,states,parameters) {
		with(as.list(c(states,parameters)),{
			beta <- beta0*(1+beta1*cos(2*pi*t/T+phi))
			dS <- N*mu - beta*S*I/N - mu*S
			dE <- beta*S*I/N - sigma*E - mu*E
			dI <- sigma*E - gamma*I - mu*I
			list(c(dS,dE,dI))
		})
	}
# Solving the system of differential equations (with deSolve package):
	states <- c(S=S,E=E,I=I)*N
	parameters <- c(N=N,mu=mu,beta0=beta0,beta1=beta1,sigma=sigma,gamma=gamma,phi=phi)
	out <- as.data.frame(ode(states,times,model,parameters))
# Ploting the prevalence as a function of time:
	if(plot) with(subset(out,time>90*365),plot(time/365,I,type="l",col="red",
		xlab="time (years)",ylab="prevalence",...))
# Optionally giving the data frame as an output:
	invisible(out)
}
