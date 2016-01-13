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
#include "evolVac.h"


using namespace std;
extern "C" {

/******************************************************************************************************/
/*
Objectif:
    Initialiser les valeurs initiales pour les parametres et pour les variables de SIMULATION
  par défaut ou  partir de terminal.

Entrée:
    Trouver valeurs à partir de terminal.

Sortie:
    Un vecteur en une dimension.
*/
double *initialervalParSIM(int argc, char** argv);
unsigned long **initialerVarSEIR(int nbVilles,int nvar,
                                 unsigned long S0, unsigned long E0, unsigned long I0, unsigned long R0);

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
vector<vector<vector<unsigned long> > > simulerSEIR(int nbVilles, double sigma, double gamma, double mu,
                                             double beta0, double beta1, double *phi,
                                             double** arr_rho,
                                             double seed, double unitTIME, double tmax, double typeRNG,
                                             double T, unsigned long **valeursSEIR0,
                                             bool valParVAC, double tauxVaccToutesVilles,
                                             vector<double>*politiqueEVOL);

/***********************************************************************************************************/
/*
Objectif:
    Fonction est utilisée dans le cas, les valeurs des parameters, et des variables sont identiques pour toutes villes.
    Fonction appellera la fonction "simulerSIER".

Entrée:
    valParSIM: contient les valeurs des parameters et des variables

Sortie:
    Un vecteur en trois dimensions:
        la dimension 1 pour temps
        la dimension 2 pour valeurs des variables
        la dimension finalles pour le nombre de villes

*/


/***********************************************************************************************************/
/*
Objectif:
    Calcule le taux de contact beta qui suit un forçage saisonier, pour toutes les sous-populations
 avec beta_i(t) = beta0 *(1+ beta1*cos(2*pi*t + phi_i).

Entrée:
    nbVilles: nombre de subpopulations.
        beta0 : valeur moyenne du taux de contact (/ind/jours).
        beta1 : l'amplitude du taux de contact (/ind/jours).
        phi_i : le phase du forçage pour la sous-population i, [0, 2*pi].
    t: au moment t de simulation.
Sortie :
        Un vecteur en 1 dimension contient les noeuds. Le nombre de noeuds est égal au nombre du sous-populations. La valeur de chaque noeud est la valeur de beta pour chaque sous-population.
 */
double *calculerBeta(int nbVilles, double beta0, double beta1,
                     double periBETA,double t, double*phi);

/******************************************************************************************************/
/*
Object:
    Faire l'événement m dans la sous-population "nextville".

Entrée:
    m: événement m
    nextville: sous-population où l'événement m se produit.
        y: table qui enregistre les nombres d'individus dans chaque groupe.
    sigma: c'est la condition pour le modèle SEIR ou le SIR
    matrCumulI: c'est la matrice qui contient le nombre d'infectés pendant une unité de temps
Sortie:
    table qui enregistre les nombres d'individus dans chaque groupe après que l'événement m a fini.
  */
void fairEvenementM(int m, double sigma, int nextville, unsigned long **y, int *matrCumulI);

/******************************************************************************************************/
/*
Object:
    Calculer les valeurs de Lamda (c'est le taux d'infection S -> E) pour toutes les sous-populations avec beta sinusoidal, le nouveau modle SEIR au moment t.

Entrée:
          nbVill: nombre de sous-populations, chiffre entier positif.
          beta: matrix en 1 dimension qui contient les valeurs de contact des toutes les sous-populations, (/ind/jours).
          rho: matrice en 2 dimension de taux de couplage entre deux sous-populations i et j, [0,1].
          y: table servant à sauvegarder le nombre de groupe S, E, I, R au moment t.
          lamda: table pour sauvegarder les valeurs lamda au moment t, (/ind/jours).
          epsilon: taux d'infection de l'exterieur, [0,1].

 Sortie:
          valeus de taux d'infection des toutes les sous-populations.
  */

void calculerLamda(int nbVilles, double **arr_xi, double *beta, double **arr_rho,
                   unsigned long **y, double *lamda);


/******************************************************************************************************/
/*
Object:
    calculer les valeurs de la fonction de propensité, dans ce cas, beta est sinusoidal.

Entrée:
        nbVill: nombre de sous-populations, chiffre entier.
        rmu: taux de naissance et taux de mortalité. (/jours).
        sigma: taux de transmission E -> I (/jours).
        gamma: taux de transmission I -> R (/jours).
        lamda: table pour sauvegarder les valeur de fonction de propensité (/ind/jours).
        y:  table pour  sauvegarder les nombres de personnes dans chaque groupe.
        f: table pour sauvegarder les valeurs de la fonction de propensité.

Sortie:
        valeurs de fonctions de propensités des toutes les sous-populations.
  */
void calculerF(int nbVilles,double rmu, double sigma, double gamma,
                  double *lamda, unsigned long **y,double **f);

/******************************************************************************************************/
 /*
Object:
    Fonction crée un table en deux dimensions qui contient les taux de couplage entre ville i et j.

Entrée:
    rho0: taux de couplage initiale, [0,1]

*/

double** calculerProbVisiter(int nbVilles, double rho0);


 /******************************************************************************************************/
 /*
Object:
    Fonction  calculer les valeurs de epsilon (taux d'infection de l'extérieur)pour entre deux population i et j

Entrée:
    epsilon0: taux infection de l'extérieur initiale,[0,1]

*/

 // double **calculerProbVoir(int nbVilles, unsigned long **y, double **arr_rho);
void calculerProbVoir(int nbVilles, unsigned long **y, double **arr_rho, double**arr_xi);


/******************************************************************************************************/
/*
 Objectif:
    Calculer le phase pour chaque la population avec phi_i = i*(phiMAX - phiMIN)/(nbVille-1), avec i de 0 à
(nbVilles-1)

 Entrée :
        nbVille : nombre de sous-populations, chiffre entier.
        phiMAX : phase maximun,[0,2*pi].
        phiMIN : phase minimum, [0, 2*pi].

 Sortie :
       Un vecteur en 1 dimension contient les noeuds.
       Le nombre de noeuds est gal au nombre du sous-populations. La valeur de chaque noeud est le phase pour 				chaque sous-population.
  */
 double *calculerPhases(int nbVilles, double phiMAX, double phiMIN);
/******************************************************************************************************/
/*
Objectif:
    Fonction est d'initialiser l'adresse pour un pointeur 1 dimention.
Entrée :
        nbVill : nombre de cases memoires qu'il faut initialiser.
        Matrice : nom de pointeur  1 dimention.
Sortie :
        un pointeur en 1 dimention.
*/
void create1D(int nbVilles,double *&vect);


 /******************************************************************************************************/
 /*
   Objectif :
    Fonction est d'initialiser l'adresse pour un pointeur à 2 dimentions.
   Entrée :
        nbVill : nombre de cases memoires qu'il faut initialiser.
        Matrice : un pointeur à deux dimentions.
   Sortie :
        un pointeur à 2 dimentions.
   */
 void create2D(int nbVilles, int nevent, double **&vect);

/******************************************************************************************************/
 /*
Objectif:
    Fonction choisit une ville où un événement va se produire.
Entrée:
    nbVille : nombre de sous-populations
        nevent : nombre d'événement qui peut se produire.
        p2 : probabilite de propensité
        f : fonctions de propensité
Sortie:
       ville choisie.
   */
int choisirVille(int nbVilles, int nevent, double p2, double **f);

/******************************************************************************************************/
 /*
Objectif:
    Fonction choisit l'événement qui se produit dans la sous-population choisie "nextville"
Entrée:
       nevent : nombre d'événement qui peut se produire dans cette sous-population.
       nextville: sous-population choisie
       p3 : probabilité s'établit pour qu'un événement peut se produire.
       f : fonctions de pensilite de la souspopulation choisie.
Sortie:
       événement choisi.
   */

int choisirEvenement(int nevent, int nextville, double p3, double **f);

/******************************************************************************************************/
/*
Objectif:
    Redimensionner un vecteur à une dimension avec la longueur donnée.
Entrée:
    V: vecteur redimensionnée
    n: longueur à la quelle le vecteur V doit arriver
Sortie:
    Vecteur à une dimension

*/
vector<double> resizeVector1D(vector<double> V,int n);

/******************************************************************************************************/
/*
Objectif:
    Stocker les resultats dans les fichier de text, un fichier par sous-population.
Entrée:
    valParSIM: vecteur qui contient les valeurs de parameters.
    tabVille: table qui stocke les résultats des villes après avoir fait la simulation.
Sortie:
    Fichiers qui retiennent les résulats des villes. Leurs noms viennent des noms des paramèters.
*/
 void writeInFile(double*valParSIM,vector<vector<vector<unsigned long> > > tabVille);


/*
Objectif:
    Créer une table à deux dimensions qui stocke le résultat total des toutes villes selon unité de temps.
    Après:
        enregistrer le résultat total à un fichier en text.
*/
  vector<vector<double> >  tableauToTal_00(vector<vector<vector<double> > > tabVille);
  void writeInFileToTal_00(string nameFile,vector<vector<double> >  tableauToTal);

/******************************************************************************************************/
  /******************************************************************************************************/
  void verifierEvolVac(vector<double> valParSIM, vector<double> valParVAC,
                       int nbEXPOL,int nbSTRAT,int nbPOL);

  /******************************************************************************************************/
  //Pour Vaccination
  int vaccinerToutesVilles_V1(double tauxVaccToutesVilles,int nbVilles, unsigned long **y,
                              double_vector *politique);

  /******************************************************************************************************/
  //Pour EvolVac
  double eval_policy(vector<double> valParSIM,
                     char SAUVER[5],vector<nodeSEIR>*valeursSEIR0,
                     vector<double>*valParVAC,double_vector* policy);

  /******************************************************************************************************/
//LA FIN

