/***************************************************************************************************************************/
// But de programme est de faire le simulateur de vaccinations ne pas utiliser l'algorithme de l'appentissage par renforcement.
// Vaccination mormale, c'est-à-dire que donner %S vaccinés à chaque point de vaccination.
// Etape1: choisir le temps ou la reaction se produit avec la probabilite p1
// Etape2: faire vaccination avec les memes choses pour toutes les villes
// Etape3: choisir une ville ou la reaction se produit avec la probabilite p2
// Etape4: choisir une reaction des 8 reactions produite avec la ville choisie dans l'etape 2 avec la probabilite p3
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
#include "evostrat.h"
#include "VaccSEIR_EvoStrat.h"
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

char SAUVER[5]="Y";
double tauxVaccToutesVilles = 0; // entre 0 et 100
int nombreMaxVaccine = 10000000; // nombre max de vaccins utilises
int nombreTotalInfectes;
int nombreTotalVaccine = 0;

/*
 *
 *
 *
 */




/*
 sequence obtenue avec les variables I et Constante:
 fitness sequence = (-6341,-4461,-4461,-4461,-4306,-4306,-4306,-4306,-4306,-4306,-4295,-4295,-4295,-4295,-4295,-4295,-4295,-4291,
 -4291,-4291,-4291,-4291,-4291,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,
 -4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4282,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,
 -4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4266,-4194,-4194,
 -4194,-4194,-4194,-4194,-4194,-4194,-4194,-4054,-4054,-4054,-4054,-4054,-4054,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,
 -4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4140,-4108,-4108,-4108,-4108,-4108,-4108,
 -4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4046,-4046,
 -4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,-4046,
 -4046,-4015,-4015,-4015,-4015,-4015,-4015,-4015,-4015,-4015,-4015,-4002,-4002,-4002,-4002,-4002,-4002,-4002,-4002,-4002,-4002,-3999
 
avec SEIR:
-4377,-4377,-4377,-4377,-4375,-4362,-4362,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4302,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4288,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,-4257,
-4393,-4370,-4323,-4323,-4323,-4323,-4323,-4301,-4301,-4301,-4301,-4301,-4301,-4301,-4265,-4265,-4265,-4265,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,-4263,
-4332,-4332,-4332,-4332,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4285,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4284,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4203,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,-4108,
-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4172,-4162,-4162,-4162,-4162,-4162,-4162,-4162,-4162,-4162,-4162,-4162,-4162,-4162,-4110,-4110,-4110,-4110,-4110,-4110,-4110,-4110,-4110,-4107,-4107,-4107,-4107,-4107,-4107,-4107,-4107,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,-3995,
-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,-4042,
-4087,-4087,-4087,-4087,-4087,-4087,-4087,-4087,-4087,-4087,-4087,-4087,-4087,-4087,-4087,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,-4052,
-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,-3952,
 
 */
double eval_policy(double_vector* policy)
{
	int trial;
	//int ntrials=100;
	int ntrials=100;
	int ev=0;
	for (trial=0; trial < ntrials; trial++) {
		//-- Make mymain be silent
		std::stringstream   redirectStream;
		std::streambuf*     old_cout  = std::cout.rdbuf( redirectStream.rdbuf() );		
		mymain(0,NULL,policy);
		//--- start talking again
		std::cout.rdbuf(old_cout);

		ev += nombreTotalInfectes;
	}
	return - ev / ntrials;
}
/*
 *
 *
 *
 */


// C'est la fonction principale
int mymain(int argc, char ** argv,double_vector* policy=NULL)
{	
    // Initialer le nombre de cas = 5 (S, E, I, R, N); le nombre d'evenements = 8
    int neq=6, nevent=9;

    // le nombre de villes utilisé
    int nbVill=1;

    // Type de Topology, si topology = 1 , alors c'est la métapopulation sans structure; si topology != 1, alors c'est la structure en cercle
    int topology = 1;

    //Taux de E->I
    //double sigma = (double) 1/8;
    //YANN MODIF
	double sigma = 0.143;
	
    //Taux de I->R
    //YANN MODIF
    //double gamma= (double) 1/5;
	double gamma = 0.143;

    //Taux de transmission
    //double beta0 = (double)1250/365;
    //YANN MODIF
	double beta0 = 2.740;

   // double beta=0, sigma2=0;

    //Taux de naissance
    //double rmu = (double)1/(70*365);
    //YANN MODIF
	double rmu = 0.00027;
	
    //Taux de transmission entre deux villes
    double epsilon =0.1;
    double nu = 0.01; //infection from outside

    //Parametres de probabilités, p1: probabilité de temps, p2: probabilité de ville, p3: probabilité de réaction qui se produit
    double p1, p2, p3;
    double tstep;

    // Initialisation les nombres de personnes dans le groupe I et le groupe N
    unsigned long rinitN=10000;
    unsigned long rinitI=100;

    // Initialisation
    //double tmax =(double)500;
	//YANN MODIF
	double tmax = 2000.0;
	
    double fsum =0.0, fsumVilles =0.0;
    double t = 0.0;
    int m=0, nextville=0;
    double fprop, fvill, fvillTemps;
    double unitTemps=1.0;
    double graine=10;
    int sizeTable=0,sizeTblDer=0;
  //  int MAXSIZEEVENTTABLE=5000;

    //Les parametres pour la vaccination
    // Pendant nombreJourVaccination, on va vacciner combien de phase, c'est à dire combien de jour par fois de vaccination?
    int periode=10;
    double tempsPolitiqueVaccination =0.0;
    int kMax=0;
    //Le temps où on commence à vacciner
    double tstart=0.0;
    //Le temps où on finit à vaccinaer
    double tfinal=tmax;

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
          else if((strcmp(str,"-beta0"))==0)
          {
              i++;
              beta0 = atof(argv[i]);
          }
          else if((strcmp(str,"-sigma"))==0)
          {
              i++;
              sigma = atof(argv[i]);
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

          else if((strcmp(str,"-periode"))==0){
               i++;
               periode = atoi(argv[i]);
           }
           else if((strcmp(str,"-tstart"))==0){
               i++;
               tstart = atof(argv[i]);
           }
           else if((strcmp(str,"-tfinal"))==0){
               i++;
               tfinal = atof(argv[i]);
           }
          else if((strcmp(str,"-tauxVaccToutesVilles"))==0){
                  i++;
                  tauxVaccToutesVilles = atoi(argv[i]);
          }

        }
    tempsPolitiqueVaccination  = tfinal-tstart;
    kMax = (int)(tempsPolitiqueVaccination /periode) +1;
    cout<<endl<<"Parametres de SIMULATION SEIR!"<<endl;
    cout<<"rmu = "<<rmu<<",  sigma= "<<sigma<<",gamma="<<gamma<<", beta0= "<<beta0<<", epsilon="<<epsilon<<" ,nu="<<nu<<endl;
    cout<<"initI="<<rinitI<<",  initN="<<rinitN<<",  nbVill="<<nbVill<<endl;
    cout<<"unitTemps="<<unitTemps<<",  tmax="<<tmax<<endl;
    cout<<"sauvergarder les resultats: "<<SAUVER<<endl;
    cout<<"topology="<<topology<<endl;
    cout<<"sigma="<<sigma<<endl;
    cout<<"graine="<<graine<<endl;

    cout<<endl<<"Parametres de VACCINATION!"<<endl;
    cout<<"Periode(Apres de ?jours on continue à vacciner) = "<<periode<<endl;
    cout<<"nombre de fois de vaccination="<<kMax<<endl;
    cout<<"temps à commencer la vaccination="<<tstart<<endl;
    cout<<"temps à finir la vaccination="<<tfinal<<endl;
    cout<<"tauxVaccToutesVilles = "<<tauxVaccToutesVilles<<endl;


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
            sprintf(nameFile,"SimVaccinationVil%dEps%.3fTop%d.csv",i,epsilon,topology);
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

    //Vaccination
     //Les parametre pour vacciner
     int tpNombreFrequence=0;
     double pointVacciner=tstart + tpNombreFrequence*periode;
     int nombreSusceptibleVaccinees=0;
	 nombreTotalVaccine=0;
	 nombreTotalInfectes=0;

    //Calculer le temps où la programme marche selon le temps de CPU
    clock_t start, end;
    start = clock();

    //Initialiser les nombres de personnes de chaque groupe
    //MODIF YANN
	//srand(graine);
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

        //Vaccination
        //Vacciner pour toutes les villes selon chaque la taux de vaccination par ville
        
		if((t>pointVacciner)&&(tpNombreFrequence<= kMax) && nombreTotalVaccine < nombreMaxVaccine){
            cout<<"t="<<t<< " pointVacciner="<<pointVacciner<<" kMax="<<kMax<<endl;
            nombreSusceptibleVaccinees=vaccinerToutesVilles_V1(tauxVaccToutesVilles,nbVill,y,policy);
			nombreTotalVaccine += nombreSusceptibleVaccinees;
            pointVacciner +=periode;
            tpNombreFrequence++;
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
                //beta = creerBetaAleatoire(beta0,sigma2);
                CalculerLamda(nbVill,topology,beta0,epsilon,y,lamda);

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
                        table[nextville][sizeTable][6] = nombreSusceptibleVaccinees;
                        nombreSusceptibleVaccinees = 0;
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
        y[nextville][iI] = rnitN*0.01;
        y[nextville][iN] = rnitN;
        y[nextville][iS] = y[nextville][iN]*0.2;//-y[nextville][iI];
        y[nextville][iE] = 0;
        y[nextville][iR] = y[nextville][iN]-y[nextville][iS]-y[nextville][iI];
    }
    return;
}

/******************************************************************************************************/
void initialerVarSEIR(int nbVill,int nvar,unsigned long **&y,
                                 unsigned long S0, unsigned long E0, unsigned long I0, unsigned long R0){
    y = new unsigned long*[nbVill];
    for (int i = 0; i <  nbVill; i++)
    {
          y[i] = new unsigned long[nvar];
          y[i][iI] = I0;
          y[i][iN] = S0+E0+I0+R0;
          y[i][iS] = S0;
          y[i][iE] = E0;
          y[i][iR] = R0;
     }

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
    for(int v=0; v<nbVill; v++){//en forme d'un graphe complet
            sumIj=0.0;//
            for(int l=0; l<nbVill;l++){
                if(l !=v){
                    tp = (double)y[l][iI];
                    sumIj= sumIj+tp;
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
       // f[nextville][0]= nu;//*y[nextville][iS];
        //Infection
        f[nextville][3] = lamda[nextville]*y[nextville][iS];
        //Susceptible death
        f[nextville][2] = rmu*y[nextville][iS];
        //Birth
        f[nextville][1] = rmu*y[nextville][iN];
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
		nombreTotalInfectes++;
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
    //int ville = ivill+1;
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
	unsigned long S, E, I, R, N;
	//unsigned long event;
    int t;
    unsigned long sizeTable;
    char line[200];
    FILE * pFile;
    int ville;
    unsigned long SIZE;
    char nameFile[100];
    sprintf(nameFile,"%s%d%s%ld%s%.3f%s%.1f%s","SimuNbVille",nbVill,"Ninit",rinitN,"Epsilon",epsilon,"Tmax",tmax,"_SUMMARY.csv");
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

