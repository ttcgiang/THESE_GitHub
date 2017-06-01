#ifndef VACOPTIMAL_H
#define VACOPTIMAL_H
//Initialiser les bibliotheques utilisees
#include "vacOptimal.h"

using namespace std;



// Vaccination
// but est de vacciner p% de la population S->R
// Avec tous les T jours,
// La même chose pour toutes les villes
int vaccinerToutesVilles_V1(double tauxVaccToutesVilles,int nbVilles, unsigned long **y,double_vector *politique)
{
    unsigned long nbPersonneVaccine=0;
    int i;

    if (politique == NULL)
        for(int i=0;i<nbVilles;i++){
            nbPersonneVaccine=(y[i][iS]*tauxVaccToutesVilles)/100;
            y[i][iS] = y[i][iS]-nbPersonneVaccine;
            y[i][iR] = y[i][iR]+nbPersonneVaccine;
        }
    else
        for(i=0;i<nbVilles;i++)
            if (y[i][iS]*(*politique)[0]+y[i][iE]*(*politique)[1]+y[i][iI]*(*politique)[2]+y[i][iR]*(*politique)[3]+100.0*(*politique)[4] > 0) {
//			if (y[i][iS]*(*politique)[0]+100.0*(*politique)[1] > 0) {
                    nbPersonneVaccine=(y[i][iS]*tauxVaccToutesVilles)/100;
                    y[i][iS] = y[i][iS]-nbPersonneVaccine;
                    y[i][iR] = y[i][iR]+nbPersonneVaccine;
            }

    return nbPersonneVaccine;
}
#endif // VACOPTIMAL_H

