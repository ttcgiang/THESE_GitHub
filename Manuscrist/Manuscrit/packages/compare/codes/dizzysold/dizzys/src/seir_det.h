#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

extern "C" {

/* initializer */
void initparms(void (* odeparms)(int *, double *));

/* Derivatives */
void derivs (int *neq, double *t, double *y, double *ydot,double *yout, int *ip);

/***********************************************************************************************/
}


