/*
1. Faire la simulation stochastique SEIR et SIR en utilisant l'algorithme Gillespie en 1967.
2. Modèle SIR est un cas spécifique du modèle SEIR. 
Si le taux de guérison gamma est infinitif,
alors, modèle SEIR devient SIR.
*/
//bibliothèque de C++
#include <ctime>
#include <stdio.h>
#include <vector>
#include <map>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <stdexcept>
#include "index.h"
#include "rangenyan.h"

using namespace std;

/******************************************************************************************************/
/*
Objectif:
	Fonction est utilisée dans le cas, les valeurs des parameters, et des variables sont différents pour chaque ville.
	Simule le modele SEIR avec plusieurs populations et beta sinusoidal.
        Chaque sous-population est identique. Les valeurs des variables sont calcules  intervalles de temps fixes (unitTemps).

Entrée: 
	nbVilles: nombre de sous-populations.
        gamma: taux de guérison (/jour).
        sigma: inverse de la durée moyenne d'infection (/jour).
        beta0: valeur moyenne du taux de contact (/ind/jour).
        beta1: phase du forçage
        tmax: durée maximale de la simulation (jours).
        unitTemps: pas de temps auquel les valeurs des variables sont enregistrées 				(jours).
        S0: valeur initiale de susceptibles dans une sous-population.
        E0: valeur initiale d'exposé dans une sous-population.
        I0: valeur initiale d'infectieux dans une sous-population.
        R0: valeur initiale de guéris dans une sous-population.
        N0: valeur initiale de la taille d'une sous-population.
        rmu: taux de mortalité égal au taux de natalité (/jour).
        epsilon: [0;1] taux d'infection de l'extrieur.
        rho: [0;1] taux de couplage entre les sous-populations
        graine: graine du générateur de nombres aléatoires.
        phiMIN: phase minimum dans la métapopulation.
        phiMAX: phase maximun dans la métapopulation.

Sortie:
	Un vecteur en trois dimensions:
		la dimension 1 pour temps
		la dimension 2 pour valeurs des variables
		la dimension finalles pour le nombre de villes
*/
vector<vector<vector<unsigned long> > > simulerSEIRGENE(int nbVilles,double sigma,double gamma,double mu,
                                             double beta0,double beta1,double *phi,
                                             double** arr_rho,
                                             double seed,double unitTIME, double tmax, double typeRNG,
                                             double T, unsigned long **valeursSEIR0);


/******************************************************************************************************/
/*
Objectif:
	Calculer les taux pour les fonctions de propensités. 
	Cette fonction établie pour aider la simulation stochastique qui utilise l'algorithme adaptive-tau.
*/
double* seirratefunc(double* x,int nbVilles,double beta0,double beta1,double mu, double sigma, double gamma,double *phi,double** arr_rho, double T, double m_T);
//LA FIN


