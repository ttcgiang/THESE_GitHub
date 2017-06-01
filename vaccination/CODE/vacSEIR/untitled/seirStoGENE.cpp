#include "seirStoGENE.h"
#include "initialFUNC.h"

#include <float.h>
  /***************************************************************************************************************/
// Vaccination
// but est de vacciner p% de la population S->R
// Avec tous les T jours,
// La même chose pour toutes les villes
int vaccinerToutesVilles_V1(double tauxVaccToutesVilles,int nbVilles, unsigned long **y)//,double_vector *politique)
{
    unsigned long nbPersonneVaccine=0;
    int i=0;

    //cout<<"Avant S="<<y[0][iS]<<" E="<<y[0][iE]<<" I="<<y[0][iI]<<" R="<<y[0][iR]<<endl;

    //if (politique == NULL)
	for(int i=0;i<nbVilles;i++){
		nbPersonneVaccine=(y[i][iS]*tauxVaccToutesVilles)/100;
		y[i][iS] = y[i][iS]-nbPersonneVaccine;
		y[i][iR] = y[i][iR]+nbPersonneVaccine;
	}
  /*  else
         for(int i=0;i<nbVilles;i++)
	if (y[i][iS]*(*politique)[0]+y[i][iE]*(*politique)[1]+y[i][iI]*(*politique)[2]+y[i][iR]*(*politique)[3]+100.0*(*politique)[4] > 0) {
//	if (y[i][iS]*(*politique)[0]+100.0*(*politique)[1] > 0) {
	nbPersonneVaccine=(y[i][iS]*tauxVaccToutesVilles)/100;
	y[i][iS] = y[i][iS]-nbPersonneVaccine;
	y[i][iR] = y[i][iR]+nbPersonneVaccine;
}

    //cout<<"Apres S="<<y[0][iS]<<" E="<<y[0][iE]<<" I="<<y[0][iI]<<" R="<<y[0][iR]<<endl;
*/

	return nbPersonneVaccine;
}




/*
Objectif:
	Fonction est utilisée dans le cas, les valeurs des parameters, et des variables sont différents pour chaque ville.
	Simuler le modèle SEIR avec plusieurs populations et beta sinusoidal.        

Entrée: 
	nbVilles: nombre de sous-populations.
        gamma: vecteur de taux de guérison (/jour).
        sigma: vecteur de inverse de la durée moyenne d'infection (/jour).
        beta0: vecteur de valeur moyenne du taux de contact (/ind/jour).
        beta1: vecteur de phase du forçage
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
  vector<vector<vector<unsigned long> > > simulerSEIRGENE(int nbVilles,double sigma,double gamma,double mu,
                                               double beta0,double beta1,double *phi,
                                               double** arr_rho,
                                               double seed,double unitTIME, double tmax, double typeRNG,
                                               double T, unsigned long **valeursSEIR0,
		//vaccination
		double tstart,double tstop,double totalVac,double periodVac,double percentVac){
      /*
    cout<<"nbVilles="<<nbVilles<<"  sigma="<<sigma<<"  gamma="<<gamma<<endl;
    cout<<" mu = "<<mu<<cout<<" beta0="<<beta0<<endl;
    cout<<"beta1 = "<<beta1<< "  seed = "<<seed<< " unitTIME="<<unitTIME<<endl;
    cout<<"tmax= "<<tmax<<endl;

    for(int i=0; i<nbVilles; i++){
          cout<<"    "<<phi[i]<<endl;
      }

    for(int i=0; i<nbVilles; i++){
        for (int j=0; j<6; j++)
            cout<<"    " << valeursSEIR0[i][j];
        cout<<endl;
    }
    */
      // Parameters for pour SIMULATION
    // Number of varibales = 5 (S, E, I, R, N);
    // Number of events = 9
      int nevent = 9; //nbVarSEIRN = 6;
    //unsigned long S0 = S, E0 =E, I0 =I, R0 = R;
    double p1 = 0.0, p2 = 0.0, p3 = 0.0;
    double t = 0.0, tstep = 0.0;
    double fsumVilles = 0.0;
    int event = 0, nextville = 0;
    double pointTemps = unitTIME;

    //Initializing necessaire vectors
    double *lamda = new double[nbVilles];
    int *matrCumulI = new int[nbVilles];
    for(int i=0; i<nbVilles; i++) matrCumulI[i]=0.0;
    double *beta_sinusoidal;
    double** arr_xi;create2D(nbVilles,nbVilles,0.0,arr_xi);
    
    // Saving the values of all variables and all propensity functions at moment t
    unsigned long **y;
    double **f;
    f = new double*[nbVilles];
    for (int nextville = 0; nextville <  nbVilles; nextville++)
    {
        f[nextville] = new double[nevent];
        for(int i = 0; i< nevent ; i++) f[nextville][i] = 0.0;
    }

    //Initialising the values of variables
    y = valeursSEIR0;

    //Creating a table that saves the values of variables when simulation runs
    vector< unsigned long> tpmSEIRNP;
    vector<vector<vector<unsigned long> > > tableResVilles(nbVilles);
    for(int i=0; i<nbVilles; i++){
        tpmSEIRNP.push_back(0.0);tpmSEIRNP.push_back(y[i][iS]); tpmSEIRNP.push_back(y[i][iE]); tpmSEIRNP.push_back(y[i][iI]);tpmSEIRNP.push_back(y[i][iR]);tpmSEIRNP.push_back(y[i][iN]);tpmSEIRNP.push_back(0);
        tableResVilles[i].push_back(tpmSEIRNP);
        tpmSEIRNP.clear();
    }
    double xpro = 0.0, fville = 0.0, fprop =0.0;

     //Vaccination
     //Les parametre pour vacciner
     int tpNombreFrequence=0;
     double pointVacciner=tstart + tpNombreFrequence*periodVac;
     int nombreSusceptibleVaccinees=0;
int	 nombreTotalVaccine=0;
int	 nombreTotalInfectes=0;
    int tempsPolitiqueVaccination  = tstop-tstart;
    int  kMax = (int)(tempsPolitiqueVaccination /periodVac) +1;



    //PROGRAM starts.......................................................

    //Calculating the time where the program runs according to the time of CPU
    clock_t start, end;
    start = clock();
    static int seed0 = seed;
    int  *idum = &seed0;

    //number 'seed'
    srand(seed);


    //simulation
    while(t<tmax){

        //Calculating Value of \beta at time t
        beta_sinusoidal=calculerBeta(nbVilles,beta0,beta1,T,t,phi);


        //calculating Value of Prob "see"
         calculerProbVoir(nbVilles, y, arr_rho,arr_xi);


        //Calculating Value of \lamda at time t
        calculerLamda(nbVilles,arr_xi,beta_sinusoidal,arr_rho,y,lamda);

        //Calculating Values of the propensity functions
        calculerF(nbVilles,mu,sigma,gamma,lamda,y,f);

        //Generating the uniform random numbers
        if(typeRNG==FAST) {// using the method of Yann
            p1 = ran1(idum,GENERATE, (char*)"a");//probability for time
            p2 = ran1(idum,GENERATE, (char*)"a");//probability for time
            p3 = ran1(idum,GENERATE, (char*)"a");//probability to choose event
        }
        else
            if(typeRNG==GOOD){// using the method of C++
                p1 = rand()/(double)RAND_MAX;//probability for time
                p2 = rand()/(double)RAND_MAX;//probability for time
                p3 = rand()/(double)RAND_MAX;//probability to choose event
            }

        //Calculating the sum of all propensity probabilities of all cities
        fsumVilles = 0.0;
        for(int i = 0; i<nbVilles; i++){
            for(int j = 0; j<nevent; j++){
                fsumVilles += f[i][j];
            }
        }

        //Calculating the time interval \tstep and updating time
        if(fsumVilles > 0.0) {
            tstep = (double)(-log(p1)/fsumVilles);          
        }
        else {
            cout<<"no event"<<endl;
            cout<<"parce que: fsumVilles = "<<fsumVilles<<endl;
            break;
        }
       
        //Vaccination
        //Vacciner pour toutes les villes selon chaque la taux de vaccination par ville
        
         if((t>pointVacciner)&&(tpNombreFrequence<= kMax) && nombreTotalVaccine < totalVac){
            cout<<"t="<<t<< " pointVacciner="<<pointVacciner<<" kMax="<<kMax<<endl;
            nombreSusceptibleVaccinees=vaccinerToutesVilles_V1(percentVac,nbVilles,y);//,policy);
			nombreTotalVaccine += nombreSusceptibleVaccinees;
            pointVacciner +=periodVac;
            tpNombreFrequence++;
        } 






        //Choosing randomly a city where event can occur
        nextville = choisirVille(nbVilles,nevent,p2,f);

        // Choosing randomly one event in the city chosen
        if(nextville != -1){
            event = choisirEvenement(nevent,nextville,p3,f);
        }
        else{
            printf ("(nextville == -1)fsumVilles = %.20f \n", fsumVilles);
            printf ("p2 = %.5f \n", p2);
            printf ("xpro = %.20f \n", xpro);
            printf ("fville = %.20f \n", fville);
             for(int i=0; i<nbVilles; i++){
                 cout<<"ville = "<<i<<endl;
                 for(int j=0; j<nevent; j++)
                     cout<<"f"<<j<<" = "<<f[i][j]<<"  ";
                 cout<<endl;
             }
        }

        //Doing the event \m
        if(event != -1){
            fairEvenementM(event,sigma,nextville,y,matrCumulI);

            //Saving the values of variables at time t in a 3D table.
            if(t> pointTemps){
                /*
                cout<<"beta_sinusoidal"<<endl;
                    for(int i=0; i<nbVilles; i++){
                          cout<<"    "<<beta_sinusoidal[i]<<endl;
                    }
                    cout<<"RHOOOOOOO"<<endl;
                    for(int i=0; i<nbVilles; i++){
                        for (int j=0; j<nbVilles; j++)
                            cout<<"    " << arr_rho[i][j];
                        cout<<endl;
                    }
                    */
               for(int i= 0; i< nbVilles; i++){
                    tpmSEIRNP.clear();
                    tpmSEIRNP.push_back(pointTemps);tpmSEIRNP.push_back(y[i][iS]); tpmSEIRNP.push_back(y[i][iE]); tpmSEIRNP.push_back(y[i][iI]);tpmSEIRNP.push_back(y[i][iR]);tpmSEIRNP.push_back(y[i][iN]);tpmSEIRNP.push_back(matrCumulI[i]);
                    tableResVilles[i].push_back(tpmSEIRNP);
                    matrCumulI[i] = 0.0;
                }
                pointTemps +=unitTIME;
            }
        }
        else{
            printf ("(event == -1)fville = %.20f \n", fsumVilles);
            printf ("p2 = %.5f \n", p3);
            printf ("xpro = %.20f \n", xpro);
            printf ("fprop = %.20f \n", fprop);
             for(int i=0; i<nbVilles; i++){
                 cout<<"ville = "<<i<<endl;
                 for(int j=0; j<nevent; j++)
                     cout<<"f"<<j<<" = "<<f[i][j]<<"  ";
                 cout<<endl;
             }
        }
        // update time t
          t = t + tstep;
       //continue
    }
    end = clock();

    // CPU time is used:
    double totalDure = (double)(end - start)/CLOCKS_PER_SEC;

    //Freeing the memory of the pointers
    delete []lamda;// 1D pointer

    // Freeing the 2D pointers
    for (int nextville = 0; nextville <  nbVilles; nextville++)
    {
        delete []y[nextville];
        delete []f[nextville];
        delete []arr_xi[nextville];
    }
    delete []y; delete[]f;delete[]arr_xi;

    // result
    return(tableResVilles);  
}

  /***** The end of the main program******/

  /*****************************************************************************/
//____________________________________________________________________________________________________
//------------rate functions( calculating transition rates for the algorithm 'adaptive-tau'--------------------------------------------
/*
Goal:
    Calculating all rates for all propensity functions.
    Be created to help stochastic simulation that uses the approximate algorithm 'adaptive-tau'
    Transfer from R to C++, and from C++ to R
    R data types (SEXP) are matched to C++ objects in a class hierarchy.
*/

   double* seirratefunc(double* x,int nbVilles,double beta0,double beta1,
                         double mu, double sigma, double gamma,double *phi,
                         double** arr_rho, double T, double m_T){

        int nevent = 9, nbVarSEIRN = 6;

        //Initialisation
        double **arr_xi;
        double *lamda = new double[nbVilles];
        double *beta_sinusoidal; create1D(nbVilles,0.0,beta_sinusoidal);
        double **f;
        f = new double*[nbVilles];
        arr_xi = new double*[nbVilles];
        for (int i = 0; i <  nbVilles; i++)
        {
            f[i] = new double[nevent];
            arr_xi[i] = new double[nevent];
            for(int j = 0; j< nevent ; j++){
                f[i][j] = 0.0;
                arr_xi[i][j] = 0.0;
            }
        }

        //Tranform a 1D vector to 2D vector **y
        unsigned long **y;
        y = new unsigned long  *[nbVilles];
        for (int i = 0; i <  nbVilles; i++)
        {
             y[i] = new unsigned long [nbVarSEIRN];
             y[i][iS] = (unsigned long)x[0*nbVilles+i];
             y[i][iE] = (unsigned long)x[1*nbVilles+i];
             y[i][iI] = (unsigned long)x[2*nbVilles+i];
             y[i][iR] = (unsigned long)x[3*nbVilles+i];
             y[i][iN] = y[i][iS]+y[i][iE]+y[i][iI]+y[i][iR];
         }

        //Caculating the values of \beta (contact rate) for all cities at time t
        beta_sinusoidal=calculerBeta(nbVilles,beta0,beta1,T,m_T,phi);

        //calculating Value of Prob "see"
        calculerProbVoir(nbVilles, y, arr_rho, arr_xi);

        //Calculating the values of \lamda (infection rate) for all cities at time t
        calculerLamda(nbVilles,arr_xi,beta_sinusoidal,arr_rho,y,lamda);

        //Calculating the values of all propensity fonctions
        calculerF(nbVilles,mu,sigma,gamma,lamda,y,f);

        //Gather the result
        vector<double> vecRes;
        for(int j = 1; j< nevent ; j++)
        {
               for (int i = 0; i <  nbVilles; i++)
               vecRes.push_back(f[i][j]);
        }

        int nbvecRes=vecRes.size();
        double *resRates = new double[nbvecRes];

        for(int i=0; i<vecRes.size();i++){
            resRates[i] = vecRes[i];
        }

        // Free the memory of all pointers
        delete []lamda;
        for (int i = 0; i <  nbVilles; i++)
        {
            delete []y[i];
            delete []f[i];
            delete []arr_xi[i];
        }
        delete []y; delete[]f; delete[]arr_xi;

      // result
      return(resRates);
    }

/************LA FIN **************************/
