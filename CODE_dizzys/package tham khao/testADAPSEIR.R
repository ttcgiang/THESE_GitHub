library(adaptivetau)
transitions = list(c(prey = +1),# trans 1: prey grows
	c(prey = -2, pred = +1), # trans 2: predation
	c(pred = -1))# trans 3: predator dies

lvRateF <- function(x, params, t) {
	return(c(params$r * x["prey"],# rate of prey growing
	params$beta * x["prey"]*x["pred"] * # rate of predation
	(x["prey"] >= 2),
	params$delta * x["pred"])) # rate of predators dying
}

set.seed(4) # set random number generator seed to be reproducible
simResults = ssa.adaptivetau(init.values = c(prey = 1000, pred = 500),
	transitions, lvRateF,
	params = list(r=10, beta=0.01, delta=10),
	tf=12)

matplot(simResults[,"time"], simResults[,c("prey","pred")], type= "l",
xlab="Time", ylab= "Counts (log scale)", log="y" )

legend("bottomleft", legend=c("prey", "predator"), lty=1:2, col=1:2)
