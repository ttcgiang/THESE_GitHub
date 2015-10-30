/* file seir_det().cpp */
#include <R.h>
#include "seir_det.h"

static double parms[8];
#define beta0 parms[0]
#define beta1 parms[1]
#define T parms[2]
#define phi parms[3]
#define mu parms[4]
#define sigma parms[5]
#define gamma parms[6]
#define N parms[7]


/* initializer */
void initparms(void (* odeparms)(int *, double *))
{
	int nparms=8;
	odeparms(&nparms, parms);
}

/* Derivatives */
#define time y[0]
#define S y[1]
#define E y[2]
#define I y[3]
#define dtime  ydot[0]
#define dS ydot[1]
#define dE ydot[2]
#define dI ydot[3]
void derivs (int *neq, double *t, double *y, double *ydot,double *yout, int *ip)
{
	double pi=3.141593;
	dtime = 1;
	double beta=beta0*(1+beta1*cos(2*pi*time/T+phi));
	dS =  N*mu - beta*S*I/N - mu*S;	
	if(sigma==INFINITY){
		//dS =  mu - beta*S*I - mu*S;
		dE = 0.0;
		dI = beta*S*I/N - gamma*I - mu*I;
	}
	else{
		//dS =  N*mu - beta*S*I/N - mu*S;
		dE = beta*S*I/N  - sigma*E - mu*E;
		dI = sigma*E - gamma*I - mu*I;
	}
}
//The end

