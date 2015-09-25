/*
 *  Created by TRAN Thi Cam Giang on 03/2013.
 *  Goal: calculate the deterministic equations of the SEIR/SIR model.
 *  By associating C++ with R.
 *  we dip R in C++
 */

//Library of R
#include <R.h>
#include "seirDetGENE.h"

//Define the station of each parameter in the vector 'parms'
static double parms[8];
#define beta0 parms[0]
#define beta1 parms[1]
#define T parms[2]
#define phi parms[3]
#define mu parms[4]
#define sigma parms[5]
#define gamma parms[6]
#define N parms[7]


/* initialize */
void initparmsGENE(void (* odeparms)(int *, double *))
{
	int nparms=8;
	odeparms(&nparms, parms);
}

/* Derivatives*/
/* There are four ordinal different equations*/
#define time y[0]
#define S y[1]
#define E y[2]
#define I y[3]
#define dtime  ydot[0]
#define dS ydot[1]
#define dE ydot[2]
#define dI ydot[3]

void derivsGENE (int *neq, double *t, double *y, double *ydot,double *yout, int *ip)
{
    //Initialize the values of parameters
	double pi=3.141593;
	dtime = 1;

    // the sinusoidal function of the contact parameter \beta
	double beta=beta0*(1+beta1*cos(2*pi*time/T+phi));

    // the ODE of \dS
	dS =  N*mu - beta*S*I/N - mu*S;	

    if(sigma==INFINITY){// for the SIR model
        // The ODE of \dE
		dE = 0.0;
        // The ODE of \dI
		dI = beta*S*I/N - gamma*I - mu*I;
	}
    else{// for the SEIR model
        // The ODE of \dE
		dE = beta*S*I/N  - sigma*E - mu*E;
        // The ODE of \dI
        dI = sigma*E - gamma*I - mu*I;
	}
}
//The end

