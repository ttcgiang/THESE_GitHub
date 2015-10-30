#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>  
using namespace std;

extern "C" {

/* equilibrium */
SEXP getEquilibrium(SEXP mu,SEXP beta,SEXP sigma,SEXP gamma);

/***********************************************************************************************/
}


