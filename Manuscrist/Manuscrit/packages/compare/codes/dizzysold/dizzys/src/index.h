#ifndef INDEX_H
#define INDEX_H

extern "C" {

// DESCRIRE LES INDEX DES PARAMETRES
/**********************************************************************************************/
// Situer les valeurs des situations
const int it=0,iS=1,iE=2,iI=3,iR=4,iN=5,iP=6;

//Situer les index des parametres  simuler
const int itmax = 0, inbVilles = 1,
          isigma = 2, igamma = 3,imu = 4,iepsilon = 5,
          itopology = 6, irho = 7,
          iunitTIME = 8, iseed = 9,
          iS0 = 10, iI0 = 11, iE0 = 12, iR0 = 13, iN0 = 14,
          ibeta0 = 15, ibeta1 = 16, iphiMIN = 17, iphiMAX =18, iNbSimu =19,
          iPeriodeBETA = 20, itypeRNG=21;


#define Pi 3.141592653589
const int GOOD=0; //using the random generator of YANN
const int FAST=1;
const int GENERATE=0;

}
/**********************************************************************************************/
#endif // INDEX_H
