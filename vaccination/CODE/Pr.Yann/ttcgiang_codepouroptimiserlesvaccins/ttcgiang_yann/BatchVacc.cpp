/***************************************************************************************************************************/
// But de programme est de faire le simulateur d'epidemie SEIR en utilisant la méthode directe de Gillespie en 1977 
// Etape1: choisir le temps ou la reaction se produit avec la probabilite p1
// Etape2: choisir une ville ou la reaction se produit avec la probabilite p2
// Etape3: choisir une reaction des 8 reactions produite avec la ville choisie dans l'etape 2 avec la probabilite p3
// Par : TRAN Thi Cam Giang, P15 -IFI
/***************************************************************************************************************************/


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
using namespace std;

// Situer les valeurs des situations
const int iS =1, iE = 2, iI = 3, iR=4, iN=5;
#define DEUX_PI ( 2.0 * 3.141592653589793238462643383279502884197169399375 ) // PI x 2
#define PI 3.141592653589
// Fonction d'initialisation la table qui contient les nombres de personnes dans chaque groupe
void InitialerY(int nbVill,int neq, unsigned long **&y, unsigned long rnitN, unsigned long rnitI);

// Fonction à calculer les valeurs de Lamda pour chaque ville
void CalculerLamda(int nbVill, int topology, double beta, double epsilon, unsigned long **y, double *lamda);

// Fonction à calculer les fonctions de propensité
//void CalculerF(int nbVill,double rmu, double sigma, double gamma, double *lamda, unsigned long **y,double **f);
void CalculerF(int nbVill,double nu,double rmu, double sigma, double gamma, double *lamda, unsigned long **y,double **f);

// Fonction à calculer le total des fonctions de propensités
double CalculerFsum(int nevent, int nextville,double **f);

// Fonction à faire le événement m
void FairEvenementM(int m, int nextville, unsigned long **y);

// Fonction à calculer le somme fsum de toutes les villes
double TotalVilles(int nbVill, int nevent, double **f);

// Fonction à sauvegarder les resultats dans les fichiers, un fichier par ville
 void WriteInFile(string nameFile,int ivill, unsigned long table[100][5000][7],unsigned long indexTable[100]);

// Fonction à afficher les dernieres lignes de resultats de chaque ville
 void AfficheDerniereValeur(int nbVill,unsigned long rinitN, double epsilon, double tmax, unsigned long table[100][5000][7],unsigned long indexTabDer[100]);

// Fonction à afficher les arguments de l'appel
void AfficheArgument(int nbVill,unsigned long rinitN, double epsilon, double tmax,int argc, char **argv);
/**
 * Retourne un nombre pseudo-aléatoire selon une loi normale de paramètres mu et sigma
 * @param mu moyenne (espérance mathématique) de la distribution
 * @param sigma écart-type de la distribution (doit être strictement positif)
 */
double creerBetaAleatoire(double beta0, double sigma2) {

        // On récupère deux nombres pseudo-aléatoires indépendants selon une loi uniforme sur l'intervalle [0;1]
        double randNumUni = ((double) rand())/((double) RAND_MAX);
        double randNumBi = ((double) rand())/((double) RAND_MAX);

        // On récupère un nombre pseudo-aléatoire selon une loi normale centrée réduite
        // (Paramètres : moyenne = 0, écart-type = 1)
        // Utilisation de l'algorithme de Box-Muller
        double cos2Pi = cos(DEUX_PI*randNumBi);
        double randNumNorm = 0.0;
        if(cos2Pi>= 0.0)
        {
            randNumNorm = sqrt(-2.0*log(randNumUni))*cos2Pi;
        }
        else{
            randNumNorm = sqrt(-2.0*log(randNumUni))*(-cos2Pi);
        }
        return (beta0 + sigma2 * randNumNorm);
}


// C'est la fonction principale
int main(int argc, char ** argv)
{

    // Initialer le nombre de cas = 5 (S, E, I, R, N); le nombre d'evenements = 9
    int neq=6, nevent=9;

    // le nombre de villes utilisé
    int nbVill=1;


    // Type de Topology :
    // topology = 0, modèle en ^iles
    // topology = 1, modèle en cercle
    // topology = 2 , modèle en réseaux

    int topology = 0;

    //Taux de E->I
    double sigma = (double) 1/8;

    //Taux de I->R
    double gamma= (double) 1/5;

    //Taux de transmission
    double beta0 = (double)1250/365;
    double beta=0, sigma2=0;

    //Taux de natalité
    double rmu = (double)1/(70*365);

    //Taux de transmission entre deux villes
    double epsilon =0.1;
    double nu = 0.01;

    //Parametres de probabilités, p1: probabilité de temps, p2: probabilité de ville, p3: probabilité de réaction qui se produit
    double p1, p2, p3;
    double tstep;

    // Initialisation les nombres de personnes dans le groupe I et le groupe N
    unsigned long rinitN=10000;
    unsigned long rinitI=100;

    // Initialisation
    double tmax =(double)500;
    double fsum =0.0, fsumVilles =0.0;
    double t = 0.0;
    int m=0, nextville=0;
    double fprop, fvill, fvillTemps;
    double unitTemps=1.0;
    double graine=10;
    char SAUVER[5]="Y";
    int sizeTable,sizeTblDer;
    int MAXSIZEEVENTTABLE=5000;



    //Parameters pour conserver les entrees
    char str[255]="";
    for (int i=1; i<argc; i++)
    {
        // Recuperer le premier point de la fonction de transformation
          strcpy(str,argv[i]);
          if ((strcmp(str,"-sigma") == 0))
          {
              i++;
              sigma =(double) atof(argv[i]);
          }
          else if((strcmp(str,"-gamma")) == 0)
          {
              i++;
              gamma = atof(argv[i]);
          }
          else if((strcmp(str,"-beta"))==0)
          {
              i++;
              beta0 = atof(argv[i]);
          }
          else if((strcmp(str,"-sigma2"))==0)
          {
              i++;
              sigma2 = atof(argv[i]);
          }
          else if((strcmp(str,"-rmu"))==0)
          {
              i++;
              rmu = atof(argv[i]);
          }
          else if((strcmp(str,"-initI"))==0)
          {
              i++;
              rinitI = atoi(argv[i]);
          }
          else if((strcmp(str,"-initN"))==0)
          {
              i++;
              rinitN = atoi(argv[i]);
          }
          else if((strcmp(str, "-nbVilles"))==0)
          {
              i++;
              nbVill= atoi(argv[i]);
          }
          else if((strcmp(str,"-topology"))==0)
          {
              i++;
              topology = atoi(argv[i]);
          }
          else if((strcmp(str,"-epsilon"))==0){
              i++;
              epsilon = atof(argv[i]);
          }
          else if((strcmp(str,"-nu"))==0){
              i++;
              nu = atof(argv[i]);
          }
          else if((strcmp(str,"-tmax"))==0){
              i++;
              tmax = atof(argv[i]);
          }
          else if((strcmp(str,"-unitTemps"))==0){
              i++;
              unitTemps = atof(argv[i]);

          }
          else if((strcmp(str,"-SAUVER"))==0){
              i++;
              strcpy(SAUVER,argv[i]);
          }

          else if((strcmp(str,"-graine"))==0){
              i++;
              graine = atof(argv[i]);
          }

        }



    cout<<"rmu = "<<rmu<<",  sigma= "<<sigma<<",gamma="<<gamma<<", beta0= "<<beta0<<", epsilon="<<epsilon<<",nu="<<nu<<endl;
    cout<<"initI="<<rinitI<<",  initN="<<rinitN<<",  nbVill="<<nbVill<<endl;
    cout<<"unitTemps="<<unitTemps<<",  tmax="<<tmax<<endl;
    cout<<"sauvergarder les resultats: "<<SAUVER<<endl;
    cout<<"topology="<<topology<<endl;
    cout<<"sigma2="<<sigma2<<endl;
    cout<<"graine="<<graine<<endl;

    // La programme marche :
    //Initialisation
    double *lamda = new double[nbVill];

    //Sauvegarder les valeurs dans chaque cas et dans chaque evenement
    unsigned long **y;
    double **f;
    f = new double*[nbVill];
    for (int nextville = 0; nextville <  nbVill; nextville++)
    {
        f[nextville] = new double[nevent];
        for(int i = 0; i< nevent ; i++) f[nextville][i] = 0.0;
    }
    //Initialiser les fichiers pour sauvegarder ler resultal par ville
    //ville i: t    S    E   I   R    N
    //        0.001 10   20  20  20   70
    //        ...    ..  ..  ..  ..   ..
    
    string listeNomFichier[nbVill];
    int nbNomVilles=0;
    if((strcmp(SAUVER,"Y"))==0){
        char nameFile[255];
        FILE * pFile;
        for(int i=1; i<= nbVill;i++){
            sprintf(nameFile,"SimVil%dEps%.3fTop%d.csv",i,epsilon,topology);
            pFile = fopen (nameFile,"w+");
            listeNomFichier[nbNomVilles]= string(nameFile);
            fclose(pFile);
            nbNomVilles++;

        }
    }


    //Creer un fichier total qui sauvegarde le total de S, E, I, R au moment t
    FILE *totalFile;
    char nameFileTotal[255];
    sprintf(nameFileTotal,"tolSEIRNbVil%dEps%.3fTop%dGam%.2fSig%.2f.csv",nbVill,epsilon,topology,gamma,sigma);
    totalFile = fopen(nameFileTotal,"w+");

    //Creer les tables qui sauvegardent des resultats quand le programme marche
    unsigned long table[nbVill][5000][7];
    unsigned long indexTable[nbVill];
    //
    unsigned long tableSauverDernierLigne[nbVill][5000][7];
    unsigned long indexTblDer[nbVill];

    for(int i=0; i<nbVill;i++){
        indexTable[i]=0;
        indexTblDer[i]=0;

    }

    //Creer une variable de temps
    double pointTemps = unitTemps;

    //Calculer le temps où la programme marche selon le temps de CPU
    clock_t start, end;
    start = clock();

    //Initialiser les nombres de personnes de chaque groupe
    srand(graine);
    InitialerY( nbVill, neq, y, rinitN,  rinitI);

    //Initialiser les valeurs des tables
     for (int i = 0; i <  nbVill; i++)
     {
            table[i][0][0] =0;
            table[i][0][1]= y[i][iN]-y[i][iI];
            table[i][0][2] = 0;
            table[i][0][3] = y[i][iI];
            table[i][0][4]= 0;
            table[i][0][5] = y[i][iN];
            table[i][0][6] = 0;
            indexTable[i]++;

            tableSauverDernierLigne[i][0][0] =0;
            tableSauverDernierLigne[i][0][1]= y[i][iN]-y[i][iI];
            tableSauverDernierLigne[i][0][2] = 0;
            tableSauverDernierLigne[i][0][3] = y[i][iI];
            tableSauverDernierLigne[i][0][4]= 0;
            tableSauverDernierLigne[i][0][5] = y[i][iN];
            tableSauverDernierLigne[i][0][6] = 0;
            indexTblDer[i]++;

     }

    //Initialer les valeurs de Lamda pour chaque ville
    CalculerLamda(nbVill, topology,beta0,epsilon,y,lamda);

    //Initialiser les valeurs pour la fonction de propensité
    CalculerF(nbVill,nu,rmu,sigma,gamma,lamda,y,f);

    //Caculer le temps où la réaction se produit
    while(t<tmax){

        //Generate uniform random numbers
        p1 = rand()/(double)RAND_MAX;
        p2 = rand()/(double)RAND_MAX;
        p3 = rand()/(double)RAND_MAX;

        //Caculer le total fsum de toutes les villes
        fsumVilles = TotalVilles(nbVill,nevent,f);

        //Determine time interval and update time
        if(fsumVilles > 0.0) {
            //tstep = -log(p1/fsumVilles)/fsumVilles;
            tstep = -log(p1)/fsumVilles;
            t = t + tstep;
        }
        else {
            cout<<"no event"<<endl;
            cout<<"parce que: fsumVilles = "<<fsumVilles<<endl;
            break;
        }

        //Selectionne aleatoirement une ville ou se produira l'evenement
        fvillTemps =0;
        for(int i=0; i<nbVill; i++){
            fvill = CalculerFsum(nevent,i,f);
            fvillTemps = fvillTemps+fvill;
            if(p2 < fvillTemps/fsumVilles) {
                nextville = i;
                break;
            }
        }


        // Selectionne aleatoirement la nature du prochain evenement dans nextville
        fprop = 0.0;
        m = -1;
        fsum = CalculerFsum(nevent,nextville,f);
        for(int j=0; j< nevent; j++){
            fprop = fprop+f[nextville][j];
            if(p3 < (fprop/fsum)) {
                m = j;
                break;
            }
        }

        // Faire l'événement m
        if(m !=-1){
                //Faire l'événement m
                FairEvenementM(m,nextville,y);

                //Recalculer les valeurs de Lamda
                beta = creerBetaAleatoire(beta0,sigma2);
                CalculerLamda(nbVill,topology,beta,epsilon,y,lamda);

                //Recalculer les valleurs de la foction de propensité
                CalculerF(nbVill,nu,rmu,sigma,gamma,lamda,y,f);

                //Recalculer la valeur de fsum
                fsum= CalculerFsum(nevent,nextville,f);

                //Caculer le total fsum de toutes les villes
                fsumVilles = TotalVilles(nbVill,nevent,f);

                //Sauvegarder la valeur de t, le nombre de personnes dans le groupe I et l'événement qui se produit 
                if(t> pointTemps){
                    sizeTable = indexTable[nextville];

                    if(sizeTable< 5000){
                        table[nextville][sizeTable][0] =pointTemps;
                        table[nextville][sizeTable][1] = y[nextville][iS];
                        table[nextville][sizeTable][2] = y[nextville][iE];
                        table[nextville][sizeTable][3] = y[nextville][iI];
                        table[nextville][sizeTable][4] = y[nextville][iR];
                        table[nextville][sizeTable][5] = y[nextville][iN];
                        table[nextville][sizeTable][6]= m;
                        indexTable[nextville]++;
                    }
                    else{
                        if((strcmp(SAUVER,"Y"))==0) {
                            WriteInFile(listeNomFichier[nextville],nextville,table,indexTable);
                        }
                        sizeTblDer = indexTblDer[nextville];
                        for(int i=0; i<7; i++){
                            tableSauverDernierLigne[nextville][sizeTblDer][i]= table[nextville][sizeTable][i];
                        }
                        indexTable[nextville]=0;
                    }

                    pointTemps +=unitTemps;

                    //Sauvegarder dans le fichier total
                    unsigned long totalS=0, totalI=0, totalE=0, totalR=0;
                    char line[255];
                    for(int k=0; k<nbVill; k++){
                        unsigned long dernierObject = indexTable[k];
                        if(dernierObject>0){
                            dernierObject = dernierObject-1;
                            totalS += table[k][dernierObject][iS];
                            totalE += table[k][dernierObject][iE];
                            totalI += table[k][dernierObject][iI];
                            totalR += table[k][dernierObject][iR];
                    }
                        else{
                            dernierObject = indexTblDer[k]-1;
                            totalS += table[k][dernierObject][iS];
                            totalE += table[k][dernierObject][iE];
                            totalI += table[k][dernierObject][iI];
                            totalR += table[k][dernierObject][iR];
                        }
                    }
                    sprintf(line, "%ld    %ld     %ld     %ld       %.0f  \n",totalS, totalE, totalI, totalR,t);
                    fputs (line,totalFile);
                }


            }

        else
        {
             // PB !!
               cout<<"Giang a pas reussi car il y a un m=-1!"<<endl;
        }

        //Le programme fini!
    }

    end = clock();

    // Le temps CPU utilisé:
    double totalDure = (double)(end - start)/CLOCKS_PER_SEC;
    cout<<"Le temps CPU utilisé est de :  "<<totalDure<<" s"<<endl;

    //Afficher les resultats
    for(int nextville=0; nextville<nbVill; nextville++){
        sizeTable = indexTable[nextville];
        sizeTblDer = indexTblDer[nextville];
        for(int i=0; i<7; i++){
            tableSauverDernierLigne[nextville][sizeTblDer][i]= table[nextville][sizeTable][i];
        }
    }
   // AfficheDerniereValeur(nbVill,rinitN,epsilon,tmax,tableSauverDernierLigne,indexTblDer);


    //Afficher tous les resultats de chaque ville
    if((strcmp(SAUVER,"Y"))==0){
        for(int nextville=0; nextville<nbVill; nextville++){
            WriteInFile(listeNomFichier[nextville],nextville,table,indexTable);
        }
    }

    //Afficher tous les arguments de l'appel
    //AfficheArgument(nbVill,rinitN,epsilon,tmax,argc, argv);

    //Libérer la mémoire
    delete []lamda;
    for (int nextville = 0; nextville <  nbVill; nextville++)
    {
        delete []y[nextville];
        delete []f[nextville];

    }
    delete []y;
    delete[]f;
   // delete[]table;
    fclose(totalFile);
    cout<<"Giang a reussi!"<<endl;
    cout<<"***********************************************************************"<<endl;

    return(0);
}

/********************************************************************************************************************************/
//Le but de cette fonction est d'initialiser les nombres de prosonnes dans chaque groupe (S, E, I, R, N) pour chaque ville
//  nbVill: nombre de villes
//  neq: nombre de groupes
//  y: table pour sauvegarder les nombres de personnes dans chaque groupe
//  rnitN: nombre de population initial
//  rnitI: nombre de personnes infectées initial
void InitialerY(int nbVill,int neq, unsigned long **&y, unsigned long rnitN, unsigned long rnitI){
    y = new unsigned long*[nbVill];
    for (int nextville = 0; nextville <  nbVill; nextville++)
    {
        y[nextville] = new unsigned long[neq];
        y[nextville][iI] = rnitI;//rnitN*0.01;
        y[nextville][iN] = rnitN;
        y[nextville][iS] = rnitN - rnitI;//y[nextville][iN]*0.2;//-y[nextville][iI];
        y[nextville][iE] = 0;
        y[nextville][iR] = 0;//y[nextville][iN]-y[nextville][iS]-y[nextville][iI];

    }
    return;
}


/*********************************************************************************************************************************/
//Le but de cette fonction est de calculer les valeurs de Lamda (c'est le taux de transmission S->E) pour chaque ville
//  nbVill: nombre de villes
//  topology: type de topology
//  beta: valeur de trasmission
//  epsilon: taux de transmission entre deux villes
//  y: table pour sauvegarder les nombres de personnes dans chaque groupe
//  lamda: table pour sauvegarder les valeurs de lamda
void CalculerLamda(int nbVill, int topology, double beta, double epsilon, unsigned long **y, double *lamda){
    double sumIj;
    double tp,tpi;
    double tp1,tp2;

    for(int v=0; v<nbVill; v++){
        if(topology == 0){// en forme en ^iles
            sumIj = 0;
        }

        if(topology==1) {//en forme de cycle
            sumIj=0.0;
            int villefinal =nbVill-1;
            if(villefinal==0){//il y a seulement une ville
                sumIj = y[v][iI];
            }
            else if(villefinal==1)//il y a deux villes
            {
                if(v==0) sumIj=y[v+1][iI];
                else sumIj =y[v-1][iI];
            }
            else //il y a plus de deux villes
            {
            if(v==0){
                tp1 = (double)y[v+1][iI];
                tp2 = (double)y[villefinal][iI];
                sumIj= sumIj+tp1+tp2;
            }
            else if(v==villefinal){
                tp1 = (double)y[v-1][iI];
                tp2 = (double)y[0][iI];
                sumIj=sumIj+tp1+tp2;
            }
            else{
                tp1 = (double)y[v-1][iI];
                tp2 = (double)y[v+1][iI];
                sumIj=sumIj+tp1+tp2;
            }
          }
        }

        if(topology ==2){//en forme d'un graphe complet
            sumIj=0.0;//
            for(int l=0; l<nbVill;l++){
                if(l !=v){
                    tp = (double)y[l][iI];
                    sumIj= sumIj+tp;
                }
            }
        }

        tpi =(double)y[v][iI];
        lamda[v] = (double)beta*(tpi+ epsilon*sumIj)/y[v][iN];
    }
}

/**********************************************************************************************************************/
//Le but de cette fonction est de calculer les valeurs de la fonction de propensité
//  nbVill: nombre de villes
//  rmu: taux de naissance et taux de mortalité
//  sigma: taux de transmission E->I
//  gamma: taux de transmission I->R
//  lamda: table à sauvegarder les valeur de lamda
//  y:  table à  sauvegarder les nombres de personnes dans chaque groupe
//  f: table pour sauvegarder les valeurs de la fonction de propensité
void CalculerF(int nbVill,double nu,double rmu, double sigma, double gamma, double *lamda, unsigned long **y,double **f){
    for(int nextville=0; nextville< nbVill;nextville++){

         //infection from outside
        //f[nextville][0]= nu;//*y[nextville][iS];


        //Birth
        f[nextville][1] = rmu*y[nextville][iN];

        //Susceptible death
        f[nextville][2] = rmu*y[nextville][iS];

        //Infection
        f[nextville][3] = lamda[nextville]*y[nextville][iS] + nu*y[nextville][iS];

        //Exposed death
        f[nextville][4] = rmu*y[nextville][iE];

        //Movement into next class
        f[nextville][5] = sigma*y[nextville][iE];

        //Infected death
        f[nextville][6] = rmu*y[nextville][iI];

        //Movement into next class
        f[nextville][7] = gamma*y[nextville][iI];

        //Recovered death
        f[nextville][8] = rmu*y[nextville][iR];
    }
}

/*************************************************************************************************************************/
//Le but de cette fonction est de calculer le total des fonctions de propensités pour une ville
//  nevent: nombre d'événements
//  nextville: ville i
//  f: table àsauvegarder les valeurs des fonctions de propensités
//  nevent: nombre d'événements
double CalculerFsum(int nevent, int ville,double **f){
    double fsum=0.0;
    for(int i=0; i< nevent;i++){
        fsum=fsum + f[ville][i];
    }
    return fsum;
}



/***********************************************************************************************************************/
//Le but de cette fonction est de calculer le somme fsum de toutes les villes

double TotalVilles(int nbVill, int nevent, double **f){
    double fsumVilles =0.0;
    for(int i = 0; i<nbVill; i++){
        for(int j = 0; j<nevent; j++){
            fsumVilles = fsumVilles + f[i][j];
        }
    }
    return fsumVilles;
}

/*********************************************************************************************************************/
//Le but de cette fonction est de faire l'événement m
//  m: événement m
//  nextville: ville nextville
//  y: table à  sauvegarder les nombres de personnes dans chaque groupe
void FairEvenementM(int m, int nextville, unsigned long **y){

  //m=0 : infection de dehors
    /*

    if((m ==0)&&(y[nextville][iS]>0)){
         y[nextville][iI]=y[nextville][iI]+1;
         y[nextville][iS]=y[nextville][iS]-1;
         nombreTotalInfectes++;

    }
    */


    // événement : une personne est née
    if(m ==1){
        y[nextville][iS]=y[nextville][iS]+1;
        y[nextville][iN]=y[nextville][iN]+1;
    }

    // événement : une personne est morte
    if(m ==2){
        y[nextville][iS] = y[nextville][iS]-1;
        y[nextville][iN]=y[nextville][iN]-1;
    }

    // événement : une personne est exposée S->E
    if(m == 3){
        y[nextville][iS]=y[nextville][iS]-1;
        y[nextville][iE]=y[nextville][iE]+1;        
    }

    // événement : une personne exposée est morte
    if(m==4){
        y[nextville][iE]=y[nextville][iE]-1;
        y[nextville][iN]=y[nextville][iN]-1;
    }

    // événement : une personne exposée est infecté E->I
    if(m == 5){
        y[nextville][iE]=y[nextville][iE]-1;
        y[nextville][iI]=y[nextville][iI]+1;
    }

    // événement : une personne infectée est morte
    if(m == 6){
        y[nextville][iI]=y[nextville][iI]-1;
        y[nextville][iN]=y[nextville][iN]-1;
    }

    // événement : une personne infecté est guérie
    if(m == 7){
        y[nextville][iI]=y[nextville][iI]-1;
        y[nextville][iR]=y[nextville][iR]+1;
    }

    // événement : une personne guérie est morte
    if(m==8){
        y[nextville][iR]=y[nextville][iR]-1;
        y[nextville][iN]=y[nextville][iN]-1;
    }
}


//***********************************************************************************************/
// Le but de cette fonction est de sauvegarder les resultats dans les fichier de text, un fichier par ville
// nameFile: le nom de fichier qui sauvegarde les resultats
// iVill: ville i
// table: table contient les resultat

 void WriteInFile(string nameFile,int ivill, unsigned long table[100][5000][7],unsigned long indexTable[100]){

    int sizeTable;
    char line[200];
    unsigned long S, E, I, R, N;
    int event;
    int t;
    FILE * pFile;
    sizeTable = indexTable[ivill];
    int ville = ivill+1;
    pFile = fopen (nameFile.c_str(),"a+");
    for(int j=0; j<sizeTable;j++){

        t = table[ivill][j][0];
        S = table[ivill][j][1];
        E = table[ivill][j][2];
        I = table[ivill][j][3];
        R = table[ivill][j][4];
        N = table[ivill][j][5];
        event = table[ivill][j][6];
        sprintf(line, "%d         %ld    %ld     %ld     %ld     %ld    %d\n",t ,S, E, I, R, N, event);
        fputs (line,pFile);
    }
        fclose (pFile);
    return;
}

//***********************************************************************************************/
// Le but de cette fonction est de sauvegarder les dernieres lignes des resultats par ville
// Le nom de fichier a la forme suivante :
//    SimuNbVille<nombre des villes>Nini<nombre de population initial>Epsion<valeur de epsilon>Tmax<valeur de tmax>_SUMMARY.csv
// nbVill: nombre de villes
// rinitI: le nombre de personnes infectees initiale
// epsilon: le taux de contact entre deux villes
// tmax: le temps maximal ou la simulation marche
// table: table contient les resultat

 void AfficheDerniereValeur(int nbVill,unsigned long rinitN, double epsilon, double tmax, unsigned long table[100][5000][7],unsigned long indexTabDer[100]){
    unsigned long S, E, I, R, N,event;
    int t;
    unsigned long sizeTable;
    char line[200];
    FILE * pFile;
    int ville;
    unsigned long SIZE;
    char nameFile[100];
    sprintf(nameFile,"%s%d%s%ld%s%.3f%s%.1f%s","Output/SimuNbVille",nbVill,"Ninit",rinitN,"Epsilon",epsilon,"Tmax",tmax,"_SUMMARY.csv");
    pFile = fopen (nameFile,"w");
    sprintf(line,"%s", "Villes;  S;    E;    I;     R;   N;   t\n");
    fputs (line,pFile);
    for(int ivill=0; ivill<nbVill; ivill++){
        SIZE = indexTabDer[ivill];
        sizeTable = SIZE-1;
        ville = ivill+1;
        S = table[ivill][sizeTable][iS];
        E = table[ivill][sizeTable][iE];
        I = table[ivill][sizeTable][iI];
        R = table[ivill][sizeTable][iR];
        N = table[ivill][sizeTable][iN];
        t = table[ivill][sizeTable][0];
        sprintf(line, "%d;      %ld;    %ld;     %ld;     %ld;     %ld;   %d  \n",ville,S, E, I, R, N,t);
        fputs (line,pFile);
    }
    fclose (pFile);
    return;
}

/*****************************************************************/
// Le but de cette fonction cree un fichier qui contient les valeurs de tous les paramètres en argument de l'appel
// Le nom de fichier a la forme suivante :
//          SimuNbVille<nombre des villes>Nini<nombre de population initial>Epsion<valeur de epsilon>Tmax<valeur de tmax>_PARAMS.csv
// Le contenu du fichier a la forme suivante:
//              Arguments;  Nom;        Valeur
//              1;          Nbville;    1
//              2;          Nini;       5000
//              .......
void AfficheArgument(int nbVill,unsigned long rinitN, double epsilon, double tmax,int argc, char **argv){
    FILE *pFile;
    char nameFile[100];
    char line[200];
    sprintf(nameFile,"%s%d%s%ld%s%.3f%s%.1f%s","SimuNbVille",nbVill,"Ninit",rinitN,"Epsilon",epsilon,"Tmax",tmax,"_PARAMS.csv");
    pFile = fopen (nameFile,"w");
    sprintf(line,"%s", "Argument;  Nom;    Valeur\n");
    fputs (line,pFile);

    int nbArgument=0;
    char str[50]="";
    for (int i=1; i<argc; i++)
    {
        // Recuperer le premier point de la fonction de transformation
          strcpy(str,argv[i]);          
          nbArgument++;
          if ((strcmp(str,"-sigma") == 0))
          {
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"sigma", argv[i]);
              fputs (line,pFile);
          }
          else if((strcmp(str,"-gamma")) == 0)
          {
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"gamma",argv[i]);
              fputs (line,pFile);
          }
          else if((strcmp(str,"-beta"))==0)
          {
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"beta", argv[i]);
              fputs (line,pFile);
          }
          else if((strcmp(str,"-rmu"))==0)
          {
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"rmu", argv[i]);
              fputs (line,pFile);
          }
          else if((strcmp(str,"-initI"))==0)
          {
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"rInitI", argv[i]);
              fputs (line,pFile);
          }
          else if((strcmp(str,"-initN"))==0)
          {
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"rInitN", argv[i]);
              fputs (line,pFile);
          }
          else if((strcmp(str, "-nbVill"))==0)
          {
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"nbVilles", argv[i]);
              fputs (line,pFile);
          }
          else if((strcmp(str,"-topology"))==0)
          {
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"topology", argv[i]);
              fputs (line,pFile);
          }
          else if((strcmp(str,"-epsilon"))==0){
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"epsilon", argv[i]);
              fputs (line,pFile);
          }
          else if((strcmp(str,"-tmax"))==0){
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"tmax", argv[i]);
              fputs (line,pFile);
          }
          else if((strcmp(str,"-unitTemps"))==0){
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"unitTemps",argv[i]);
              fputs (line,pFile);

          }
          else if((strcmp(str,"-SAUVER"))==0){
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"SAUVER",argv[i]);
              fputs (line,pFile);
          }
          else if((strcmp(str,"-graine"))==0){
              i++;
              sprintf(line,"%d;  %s;  %s\n",nbArgument,"graine",argv[i]);
              fputs (line,pFile);
          }
        }
    fclose(pFile);
    return;
}
/********************************************************************************************/
