library("dizzys")
for(i in 1:10) {
if(i==1) plot(obj<-seir(duration=10*365,N=1e5,beta0=1000/365,beta1=0.0,sigma=0.2,gamma=0.125,mu=3.913894e-05,seed=as.numeric(sample(1:1000,1,replace=F)),method="direct"),unitTIME=365,ylim=range(1,200),xlab="time (years)",ylab="prevalence per year", main="SEIR model with package dizzys")
else
lines(obj<-seir(duration=10*365,N=1e5,beta0=1000/365,beta1=0.0,sigma=0.2,gamma=0.125,mu=3.913894e-05,seed=as.numeric(sample(1:1000,1,replace=F)),method="direct"),unitTIME=365)
}

