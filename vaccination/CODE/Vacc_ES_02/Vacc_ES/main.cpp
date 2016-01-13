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
#include "stdio.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include<sys/stat.h>
#include<sys/types.h>
#include <time.h>
#include <stdlib.h>
#include "evostrat.h"


using namespace std;

int main (int argc, char* argv [])
{
    double eval_policy(double_vector* policy);
    char SAUVER[5];
    SAUVER[0]='N';
    int i,j;

    cout << "----\n";
    // mon objectif: reussir a diminuer le nombre d'infectes en utilisant 2000 vaccins au plus.
    int nombreMaxVaccine = 2000;
    int tauxVaccToutesVilles = 80;

    // evaluate policy
    double_vector politique(5),fitness_sequence;

    for (i=0; i < 20; i++) {
        fitness_sequence.clear();
        politique = evo_strat(politique,&eval_policy,100,1.0,true,&fitness_sequence);

        cout << "fitness sequence = (";
        for (j=0; j < fitness_sequence.size(); j++)
            cout << fitness_sequence[j] << ",";
        cout << ")\n";
    }

    cout << "that's it !!!\n";
    return 0;
}
// avec les parametres actuels, je vaccine 200 fois, c'est a dire une fois tous les 10 jours, et il y a 2000 jours


