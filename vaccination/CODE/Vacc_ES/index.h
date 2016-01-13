#ifndef INDEX_H
#define INDEX_H



// Situer les valeurs des situations
const int it = 0, iS = 1, iE = 2, iI = 3, iR = 4, iN = 5, ievent = 6, iV = 7; // Nhan ajouté iV pour compartiment V

//Situer les index des parametres à simuler

// noms pour les index des parametres

// noms pour les index des parametres
const int itmax = 0, inbVilles = 1,
          isigma = 2, igamma = 3,imu = 4,iepsilon = 5,
          irho = 6,iunitTIME = 7, iseed = 8,
          iS0 = 9, iI0 = 10, iE0 = 11, iR0 = 12, iN0 = 13,
          ibeta0 = 14, ibeta1 = 15, iphiMIN = 16, iphiMAX =17, inbSimu =18,
          iPeriBETA = 19, itypeRNG=20,iV0 = 19; // Nhan ajouté iV pour compartiment V


//Situer les index des parametres à vaccination
const int itstartVAC = 0, itfinalVAC = 1, iPeriode = 2, itauxglobVAC = 3, inbMaxVaccins = 4;

//Les parametres d'exporer
const int inbEXPLOR = 0, inbSTRAT = 1,//nombre de strategies sont expolres dans une fois d'exploration
          inbPOLIT = 2;//nombre de politiques essayees dans une strategie

//Creer un type de node nodeSEIR
//Creer une structure en Node
struct nodeSEIR{
    double t;
    unsigned long S,E,I,R,N,V; // Nhan ajouté compartiment V
    int event;
};

#define DEUX_PI ( 2.0 * 3.141592653589793238462643383279502884197169399375 ) // PI x 2
// définir le numéro pi
#define Pi 3.141592653589

// noms pour des générateurs de nombres aléatoires
const int GOOD=0; // utiluser le générateur de C++
const int FAST=1; // utiliser le générateur de Yann
const int GENERATE=0; // parameter dans la fonction de Yann

#endif // INDEX_H

