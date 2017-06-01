#ifndef VARSTATIONAIRE_H
#define VARSTATIONAIRE_H
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <fstream>
#include <sstream>
#include <limits>
#include <stdexcept>
#include <random>
#include <iostream>
#include "initialFUNC.h"

using namespace std;

/*using a Random number generator (uniforme distribution)
Goal: generating a random matrix n*n (transition matrix)
Input:
    nrow : both number row and number column of the matrix
    seuil: condition, a random number in [0,seuil]
    min, max: [min, max] is the interval of the generator.

Output:
    a random matrix nrow*nrow
*/
double **createUnifNumb(int nrow, double seuil, double min, double max );
/**********************************/
/*
 Goal: finding a stationary distribution from the transition matrix by by just multiplying the transition
matrix by itself over and over again until it finally converges.
 Input:
    matrix: the transition matrix
    nrox : both number row and number column of the matrix
    power : number time we multiply the transition matrix by itself
 Output:
    a matrix for the stationary distribution. (all rows are the same values)

 */
double** exponentMatrix(double** matrix,int nrow, int power);
/**********************************/
/*
 Goal: finding the vecteur of the stationary distribution
 Input:
    matrix: matrix of the stationary distribution
    nrow : both number row and number column of the matrix

 Output:
    vecteur of the stationary distribution

 */
double* stationarydist(double** matrix,int nrow);
/**********************************/
/*
 Goal: calculating the contact rate for each city following the fomula YANN et GIANG

 Input:
    nbVilles : number of cities in the metapopulation
    prbVisit : \rho_{i,j}, the probability that an individual from subpopulation i
  visits subpopulation j.
    mtVarb : matrix of all variables S, E, I, R at the time t.
    nbCont : \kappa_{j} is the average number of contacts per unit of time a susceptible will have when visiting city j.
    prbInf : c_{i,k} is the probability that a susceptible individual native from i being in contact with another infected individual native from k
  gets infected.
    prbFromk: \xi_{jk} refers to the probability that an individual y meeting x in C_{j} comes from C_{k}.

Output:
    a vector of the contact rates for all cities in the metapopulation.
   */

double* calLamdaYANNGIANG(int nbVilles,int nbVAR, double** prbVisit,unsigned long** mtVarb,
                          double* nbCont, double** prbInf, double** prbFromk);
/**********************************/

#endif // VARSTATIONAIRE_H
