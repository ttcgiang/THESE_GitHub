#include "batchVac_01.h"
#include "index.h"
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>


int main(int argc, char *argv[])
{
    //Variable
    char SAUVER[5] = "";
    double *valParSIM, valParVAC;
    SAUVER[0] = 'N';
    valParSIM = initialervalParSIM(argc, argv);


    //cin >> politique[0] >> endl;

    double_vector politique(5);
    //cout << politique.size() << endl;

    for(int i = 0; i < politique.size(); i ++) {
        cin >> politique[i];
    }
    //for(int i = 0; i < politique.size(); i ++) {
    //    politique[i] = (double)politique[i] / (double)10.0;
    //}
     srand(time(NULL));
    //cout<<"\nSimulation avec vaccination"<<endl;
    cout << eval_policy(valParSIM,SAUVER,NULL,&valParVAC,&politique) <<endl;
    //verifierEvolVac(valParSIM,valParVAC,2,5,100);

}
