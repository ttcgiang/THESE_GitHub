# Parameters and initial conditions.
parameters <- list(beta0=17/13, beta1=0.1, gamma=1/13, omega=2*pi/365, mu=1/(50 * 365))
initials <- c(S=1/17, I=1e-4, R=1-1/17-1e-4)
# Solve the system.
sir.sinusoidal.forcing <- SIRSinusoidalForcing(pars = parameters,init = initials,time = 0:(60 * 365))
PlotMods(sir.sinusoidal.forcing)

# Solve bifurcation dynamics for 20 years.
# If max(time) < 3650, bifurcation dynamics are solved for 3650 time-steps.
parameters2 <- list(beta0=17/13, beta1=seq(0.001,0.251,by = 0.001), gamma=1/13, omega=2*pi/365, mu=1/(50 * 365))
# Uncomment the following lines:
# bifur <- SIRSinusoidalForcing(pars = parameters2,
#  				init = initials,
#				time = 0:(20 * 365))
# PlotMods(bifur, bifur = TRUE)

##SEIR model
# Parameters and initial conditions.
parameters <- c(mu = 1/(70 * 365), beta=520/365, sigma=1/14, gamma=1/7)
initials <- c(S = 0.1, E = 1e-04, I = 1e-04, R = 1 - 0.1 - 1e-4 - 1e-4)
# Solve and plot.
seir <- SEIR(pars = parameters, init = initials, time = 0:(60 * 365))

PlotMods(seir)

###SEIR4AgeClasses
# Parameters and initial conditions.
parameters <- list(beta = matrix(c(2.089, 2.089, 2.086, 2.037,
				2.089, 9.336, 2.086, 2.037,
				2.086, 2.086, 2.086, 2.037,
				2.037, 2.037, 2.037,2.037),
				nrow = 4, ncol = 4),
				sigma = 0.125, gamma = 0.2,
				mu = c(0, 0, 0, 1) / (55 * 365),
				nu = c(1 / (55 * 365), 0, 0, 0),
				n = c(6, 4, 10, 55) / 75)
initials <- c(S=c(0.05, 0.01, 0.01, 0.008),
		E=c(0.0001, 0.0001, 0.0001, 0.0001),
		I=c(0.0001, 0.0001, 0.0001, 0.0001),
		R=c(0.0298, 0.04313333, 0.12313333, 0.72513333))
#Solve and plot.
#Uncomment the following lines (running it
seir4.age.classes <- SEIR4AgeClasses(pars= parameters,
					init= initials,
					time= 0:36500)
PlotMods(seir4.age.classes,variables = c('I1' , 'I2' , 'I3', 'I4'), grid=F)

###SIR
# Parameters and initial conditions.
parameters <- c(beta = 1.4247, gamma = 0.14286)
initials <- c(S = 1 - 1e-06, I = 1e-06, R = 1 - 1e-06 - 1e-06)
# Solve and plot.
sir <- SIR(pars = parameters, init = initials, time = 0:70)
PlotMods(sir)
## natalité/mortalité
# Parameters and initial conditions.
parameters <- c(mu = 1/(70 * 365),
beta = 520/365, gamma = 1/7)
initials <- c(S=0.1, I=1e-4, R=1 - 0.1 - 1e-4)
# Solve and plot.
sir.birth.death <- SIRBirthDeath(pars = parameters, init = initials,time = 0:(60 * 365))
PlotMods(sir.birth.death)

##
# Parameters and initial conditions.
parameters <- c(beta = 1, gamma = 1 / 10, mu = 5e-4)
initials <- c(X = 500, Y = 25, N = 5e3)
# Solve and plot.
system.time(sir.demog.stoch <- SIRDemogStoch(pars = parameters,init = initials, time = 2 * 365))
PlotMods(sir.demog.stoch)



