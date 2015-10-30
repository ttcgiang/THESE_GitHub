/* equilibrium */ 
#include "equilibrium.h"
#include <stdexcept>
  
SEXP getEquilibrium(SEXP mu0,SEXP beta0,SEXP sigma0,SEXP gamma0){
	double mu, beta, gamma, sigma;
	mu = NUMERIC_VALUE(mu0);
	beta = NUMERIC_VALUE(beta0);
	sigma = NUMERIC_VALUE(sigma0);
	gamma = NUMERIC_VALUE(gamma0);
 	
	SEXP rsei_eq;
   	double *csei_eq;
   	int len = 3;

   	// Allocating storage space:
   	PROTECT(rsei_eq = NEW_NUMERIC(len));
	csei_eq = NUMERIC_POINTER(rsei_eq);
	double sir_R0 = beta/(gamma+mu);
	if(sigma==INFINITY){
		if(sir_R0>1.0){
			csei_eq[0] = 1/sir_R0;
			csei_eq[1] = 0.0;
			csei_eq[2] =  mu*(sir_R0 - 1)/beta;
		}
		else{
			 error("The equilibrium value R0 is less than 1, R0 = beta/(gamma+mu)");
		}
			
	}
	else{
		csei_eq[0] = (gamma+mu)*(sigma+mu)/(beta * sigma);
		csei_eq[1] = mu*((1/(sigma+mu)) - ((gamma+mu)/(beta * sigma)));
		csei_eq[2] = mu*((beta*sigma - (gamma+mu)*(sigma+mu))/(beta *(gamma+mu)*(sigma+mu)));
	} 
   	

   UNPROTECT(1);
   return rsei_eq;

}

//the end
