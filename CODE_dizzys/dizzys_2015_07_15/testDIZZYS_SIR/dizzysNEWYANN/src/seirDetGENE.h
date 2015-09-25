#include <stdio.h>
#include <stdlib.h>
using namespace std;

extern "C" {

/* initializer */
void initparmsGENE(void (* odeparms)(int *, double *));

/* Derivatives */
void derivsGENE(int *neq, double *t, double *y, double *ydot,double *yout, int *ip);

/***********************************************************************************************/
}


