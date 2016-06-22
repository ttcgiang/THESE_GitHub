//#ifndef __00_SIMULATIONSEIR_H__
//#define __00_SIMULATIONSEIR_H__
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
#include <stdexcept>
#include "index.h"
#include "rangenyan.h"
using namespace std;

extern "C" {


/******************************************************************************************************/
//Initialiser les valeurs initiales pour les parametres de SIMULATION
//par dfaut ou  parti du terminal
vector<double> initialervalParSIM_00(int, char**);

/******************************************************************************************************/
/*
Simule le modele SEIR avec plusieurs populations et beta sinusoidal.
                Chaque sous-population est identique. Le taux de contact beta est constant. Les valeurs
                des variables sont calcules  intervalles de temps fixes (unitTemps).
                INPUT:
                        nbVilles: nombre de sous-populations.
                        gamma: taux de gurison (/jour).
                        sigma: inverse de la dure moyenne d'infection (/jour).
                        beta0: valeur moyenne du taux de contact (/ind/jour).
                        beta1: phase du forage
                        tmax: dure maximale de la simulation (jours).
                        unitTemps: pas de temps auquel les valeurs des variables sont enregistres 				(jours).
                        S0: valeur initiale de susceptibles dans une sous-population.
                        E0: valeur initiale d'exposs dans une sous-population.
                        I0: valeur initiale d'infectieux dans une sous-population.
                        R0: valeur initiale de guris dans une sous-population.
                        N0: valeur initiale de la taille d'une sous-population.
                        rmu: taux de mortalit gal au taux de natalit (/jour).
                        topology: type de modle
                                0: pas de structure spatiale explicite. .
                        epsilon: [0;1] taux d'infection de l'extrieur.
                        rho0: [0;1] taux de couplage entre les sous-populations
                        graine: graine du gnrateur de nombres alatoires.
                        phiMIN: phase minimum dans la mtapopulation.
                        phiMAX: phase maximun dans la mtapopulation.

                OUTPUT:
            Output/BatVACVilX_NbVilXTmaxXSigmaXGammaXRmuXEpsilonXTopXRhoXUnittempsXGraineXSXEXIXRXNXbeta0_Xbeta1_XphiMINXphiMAXXnbSimuX_SimSYN_00.csv:
                                fichier de rsultat par sous-population. Dans chaque fichier, il y a les
                                valeurs de temps, S, E, I, R, N dans 6 colonnes diffrentes.

            Output/BatchtolSEIRNbVilXTmaxXSigmaXGammaXRmuXEpsilonXTopXRhoXUnittempsXGraineXSXEXIXRXNXbeta0_Xbeta1_XphiMINXphiMAXXnbSimuX_SimSYN_00.csv:
                                fichier de rsultat pour toute la mtapopulation. Dans chaque fichier,
                                il y a les valeurs de temps, S, E, I, R, N dans 6 colonnes diffrentes.
  */
//vector<vector<vector<double> > > simulerSEIR(vector<double> valParSIM, vector<vector<double> >valeursSEIR0);
vector<vector<vector<double> > > simulerSEIR(int nbVilles,vector<double>sigma,vector<double>gamma,vector<double>mu,
                                             vector<double>beta0,vector<double>beta1,vector<double>phi,
                                             double** arr_rho, double** epsilon,
                                             double seed,double unitTIME, double tmax, double typeRNG,
                                             double T, unsigned long S,  unsigned long E,  unsigned long I,  unsigned long R,
                                             vector<vector<double> >valeursSEIR0);
 vector<vector<vector<double> > > resSimSEIR(vector<double> valParSIM,vector<vector<double> >valSEIR0);

/***********************************************************************************************************/
/*
 initialise le nombre d'individus dans chaque groupe (S, E, I, R, N) pour chaque sous-population.
 Chaque sous-population est identique.
 INPUT :
      nbVill: nombre de sous-populations, chiffre entier positif.
        S0: valeur initiale de susceptibles dans une sous-population.
        E0: valeur initiale d'exposs dans une sous-population.
        I0: valeur initiale d'infectieux dans une sous-population.
        R0: valeur initiale de guris dans une sous-population.
        N0: valeur initiale de la taille d'une sous-population.
 OUTPUT :
      Un matrice en 2 dimension pour toutes sous-populations.
      Le nombre de lignes est	gale au nombre de sous-populations.
      Ensuite, il y a les valeurs de S, E, I, R, N dans 5 colonnes diffrentes.
      */
void initialerY_00(int nbVill,int neq, unsigned long **&y, unsigned long S0, unsigned long E0, unsigned long I0,
                                                        unsigned long R0, unsigned long N0);


/***********************************************************************************************************/
/*
Calcule le taux de contact beta qui subit un forcage saisonier, pour toutes les sous-populations
 avec beta_i(t) = beta0 *(1+ beta1*cos(2*pi*t + phi_i).
 INPUT :
        beta0 : valeur moyenne du taux de contact (/ind/jours).
        beta1 : l'amplitude du taux de contact (/ind/jours).
        phi_i : le phase du forage pour la sous-population i, [0, 2*pi].
 OUTPUT :
        un vecteur en 1 dimension contient les noeuds. Le nombre de noeuds est gal au nombre du sous-populations. La valeur de chaque noeud est la valeur de beta pour chaque sous-population.
 */
//void calculerBeta_00(double *beta, int nbVill, double beta0, double beta1, double periodeBETA,double t, double *phi);
vector<double>  calculerBeta_00(int nbVilles, vector<double> beta0, vector<double> beta1,
                     double periodeBETA,double t, vector<double> phi);

/******************************************************************************************************/
/*
faire l'vnement m dans la sous-population nextville.
        INPUT:
                m: vnement m
                nextville: sous-population o l'vnement m se produit.
                y: table   sauvegarder les nombres de personnes dans chaque groupe.
        OUTPUT:
                table   sauvegarder les nombres de personnes dans chaque groupe aprs que l'vnement m a fini.
  */
//void fairEvenementM_00(int m, int nextville, unsigned long **y, double *matrCumulI);
//void fairEvenementM_00(int m, double sigma, int nextville, unsigned long **y, double *matrCumulI);
void fairEvenementM_00(int m, vector<double> sigma, int nextville, unsigned long **y, double *matrCumulI);

/******************************************************************************************************/
/*
calculer le somme fsum de toutes les sous-populations dans la mtapopulation.

  */
double calculerFsumToutesVilles_00(int nbVill, int nevent, double **f);

/******************************************************************************************************/
/*calculer les valeurs de Lamda (c'est le taux de contact S  E) pour toutes les sous-populations avec beta sinusoidal, le nouveau modle SEIR au moment t.
  INPUT:
          nbVill: nombre de sous-populations, chiffre entier positif.
          beta: matrix en 1 dimension qui contient les valeurs de contact des toutes les sous-populations, (/ind/jours).
          rho: matrice en 2 dimension de taux de couplage entre deux sous-populations i et j, [0,1].
          y: table servant  sauvegarder le nombre de groupe S, E, I, R au moment t.
          lamda: table pour sauvegarder les valeurs lamda au moment t, (/ind/jours).
          epsilon: taux d'infection de l'exterieur, [0,1].

  OUTPUT:
          valeus de taux de contact des toutes les sous-populations.
  */
//void calculerLamda_00(int nbVilles, double **epsilon, double *beta, double **rho,
  //                   unsigned long **y, double *lamda);
  void calculerLamda_00(int nbVilles, double **epsilon, vector<double>beta, double **rho,
                     unsigned long **y, double *lamda);


/******************************************************************************************************/
/*
calculer les valeurs de la fonction de propensit, dans le cas, beta est sinusoidal.
        INPUT:
                nbVill: nombre de sous-populations, chiffre entier.
                rmu: taux de naissance et taux de mortalit. (/jours).
                sigma: taux de transmission E  I (/jours).
                gamma: taux de transmission I  R (/jours).
                lamda: table  sauvegarder les valeur de lamda (/ind/jours).
                y:  table   sauvegarder les nombres de personnes dans chaque groupe.
                f: table  sauvegarder les valeurs de la fonction de propensit.

        OUTPUT:
                valeurs de fonctions de propensits des toutes les sous-populations.
  *///
//void calculerF_00(int nbVill,double nu,
  //             double sigma, double gamma, double *lamda,
    //           unsigned long **y,double **f);
  void calculerF_00(int nbVilles,vector<double> rmu, vector<double> sigma, vector<double> gamma,
                    double *lamda, unsigned long **y,double **f);

/******************************************************************************************************/
/*
 calculer le phase pour chaque la population avec phi_i = i*(phiMAX  phiMIN)/nbVille.
 INPUT :
        nbVille : nombre de sous-populations, chiffre entier.
        phiMAX : phase maximun,[0,2*pi].
        phiMIN : phase minimum, [0, 2*pi].
 OUTPUT :
        un vecteur en 1 dimension contient les noeuds.
       Le nombre de noeuds est gal au nombre du sous-populations. La valeur de chaque noeud est le phase pour 				chaque sous-population.
  */
//double *calculerPhases_00(int nbVill, double phiMAX, double phiMIN);
 // double calculerTauxCouplage_00(int nbVill, double rho0);
double** calculerTauxCouplage_00(int nbVilles, double rho0);


 /******************************************************************************************************/
 // Fonction  calculer les valeurs de epsilon (taux d'infection de l'extérieur)pour entre deux population i et j
  double **calculerTauxExterieur_00(int nbVill, double epsilon0);


/******************************************************************************************************/
/*Calcule les valeurs de taux de couplage entre deux population i et j.
De nos jours, chaque valeur est identique.
    INPUT :
        nbVille : nombre de sous-populations, chiffre entier.
        rho0 : valeur initiale pour le taux de couplage [0,1].
    OUTPUT :
        un vecteur est en 2 dimensions.
        Le nombre de lignes est gal au nombre de colonnes et mme au nombre de sous-populations.
        Chaque case contient la valeur e taux de couplage entre la sous-population j et la sous-population j correspondantes.
*/

//double **calculerTauxCouplage_00(int nbVill, double rho0);
    vector<double> calculerPhases_00(int nbVilles, double phiMAX, double phiMIN);
/******************************************************************************************************/
/*
  Fonction est de creer un fichier pour sauvegarder le nombre total de chaque groupe S, E, I, R au moment t.
    INPUT:
        nameFileTotal : nombre de fichier.
    OUTPUT :
        un fichier.
  */
 FILE * creerFileToTal_00(char nameFileTotal[255]);

 /******************************************************************************************************/
 /*
   Fonction est d'initialiser l'adresse pour un pointeur  1 dimention.
   INPUT :
        nbVill : nombre de cases memoires qu'il faut initialiser.
        Matrice : un pointeur  1 dimention.
   OUTPUT
        un pointeur  1 dimention.
   */

 void creerAdresse1Point_00(int nbVill,double *&Matrice);

 /******************************************************************************************************/
 /*
   Fonction est d'initialiser l'adresse pour un pointeur  2 dimention.
   INPUT :
        nbVill : nombre de cases memoires qu'il faut initialiser.
        Matrice : un pointeur  1 dimention.
   OUTPUT
        un pointeur  2 dimention.
   */
void creerAdresse2Point_00(int nbVill, int nevent, double **&Matrice);

/******************************************************************************************************/
/*
  Fonction est de sauvegarder les valeurs des variables globales au moment t dans un fichier.
  INPUT:
       pFile : nombre de fichier.
       nbVill : nombre de sous-populations.
       t : au moment t.

  OUTPUT
       enregistrer les valeurs des variables S, E, I, R au moment t sur le fichier pFile.
  */
void sauvegarderVariablesGlobales_00(FILE* pFile, int nbVill, double t, unsigned long table[100][5000][7],unsigned long indexTable[100], unsigned long indexTblDer[100]);

/******************************************************************************************************/
/*
Fonction : creer une liste de nom de fichier de resultats.
Ces fichiers sont dans le rpertoire Output.
*/
//string *creerListeNomFichier(int nbVill, char nameFilePAM[255]);
string *creerListeNomFichier_00(int nbVill,unsigned long **y, char nameFilePAM[255]);

/******************************************************************************************************/
 /*
   Fonction chosit une ville qui fait des evenements.
   INPUT:
       nbVille : nombre de sous-populations
       nevent : nombre d'evenement qui peut se produire.
       p2 : probabilite.
       f : fonctions de pensilite
    OUTPUT:
       ville choisie.
   */

int choisirVille(int nbVilles, int nevent, double p2, double **f);

/******************************************************************************************************/
 /*
   Fonction chosit l'evement qui se produit dans la sous-population "nextville"
   INPUT:
       nevent : nombre d'evenement qui peut se produire dans cette sous-population.
       nextville: sous-population o on est l:
       p3 : probabilite.
       f : fonctions de pensilite
    OUTPUT:
       evenement choisi.
   */

int choisirEvenement(int nevent, int nextville, double p3, double **f);

/******************************************************************************************************/
/*
  Fonction cre un chiffre entre 0 et 1 au hasard.
  */
double creerGraine_00();

vector<double> resizeVector1D(vector<double> V,int n);

/******************************************************************************************************/
/*
sauvegarder les resultats dans les fichier de text, un fichier par sous-population.
*/
  void writeInFile_00(vector<double> valParSIM,vector<vector<vector<double> > > tabVille);

  vector<vector<double> >  tableauToTal_00(vector<vector<vector<double> > > tabVille);
  void writeInFileToTal_00(string nameFile,vector<vector<double> >  tableauToTal);

/******************************************************************************************************/
double* seirratefunc(double* x,int nbVilles,vector<double>beta0,vector<double> beta1,
                         vector<double> mu, vector<double> sigma, vector<double> gamma,vector<double>phi,
                         double** arr_rho, double** epsilon, double T, double m_T);
//#endif // !__00_SIMULATIONSEIR_H_


