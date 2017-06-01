#ifndef VACOPTIMAL_H
#define VACOPTIMAL_H
//Initialiser les bibliotheques utilisees
#include <ctime>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <list>
#include <map>
#include <stdlib.h>
#include <math.h>
#include<sys/stat.h>
#include<sys/types.h>
#include <time.h>
#include <stdlib.h>

using namespace std;



// Vaccination
// but est de vacciner p% de la population S->R
// Avec tous les T jours,
// La même chose pour toutes les villes
int vaccinerToutesVilles_V1(double tauxVaccToutesVilles,int nbVilles, unsigned long **y,double_vector *politique);
#endif // VACOPTIMAL_H
