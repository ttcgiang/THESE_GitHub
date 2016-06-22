#include "seir_stoch.h"


//generic varibale


/******************************************************************************************************/
//Initialiser les valeurs initiales pour les parametres de SIM
vector<double> initialervalParSIM_00(int argc, char** argv){
    vector<double> res;
    // temps maximum de simulation
    double tmax = 10.0; res.push_back(tmax);

    // nombre de sous-populations dans une metapopulation
    double nbVilles = 1.0;  res.push_back(nbVilles);

    // taux d'exposition E -> I
    double sigma = (double)1/7; res.push_back(sigma);

    // taux d'exposition I -> R
    double gamma = (double)1/7; res.push_back(gamma);

    // taux de mortalite et natalite
    double rmu = (double)1/(70*365); res.push_back(rmu);
    // taux d'infection de l'extrieur
    double epsilon = 0.0;  res.push_back(epsilon);

    // type de modle de struction spatiale
    double topology = 0; res.push_back(topology);

    // taux de contact entre deux population
    double rho = 0.000; res.push_back(rho);

    // pas de temps auquel les valeurs des variables sont enregistres (indexTEMPS)
    double unitTIME = 1.0; res.push_back(unitTIME);

    // graine du gnrateur de nombres alatoires.
    double seed = 10; res.push_back(seed);

    // valeur initiale de susceptibles dans une sous-population.
    double S0 = 741559;   double E0 = 2794; double I0 = 1675; double R0 = 9253972; double N0 = S0+E0+I0+R0;
    res.push_back(S0); res.push_back(E0); res.push_back(I0); res.push_back(R0); res.push_back(N0);

    // valeur moyenne du taux de contact (/ind/indexTEMPS).
    double beta0 = (double)1000/365; res.push_back(beta0);

    // beta1: phase du forage
    double beta1 = 0.1; res.push_back(beta1);

    // phiMIN: phase minimum dans la mtapopulation.
    double phiMIN = 0.0; res.push_back(phiMIN);

    // phiMAX: phase maximun dans la mtapopulation.
    double phiMAX = 0.0; res.push_back(phiMAX);

    // nombre de fois de simulation
    double nbSimu = 5; res.push_back(nbSimu);

    // periode dans la fionction de BETA
    double periodeBETA = 365; res.push_back(periodeBETA);

    // periode dans la fionction de BETA
    double typeRNG = 0; res.push_back(typeRNG);

   //Parameters pour conserver les entrees
    char str[255] = "";
   for (int i=0; i<argc; i++)
       {
           // Recuperer le premier point de la fonction de transformation
             strcpy(str,argv[i]);
             if ((strcmp(str,"-sigma") == 0))
             {
                 i++;
                 res[isigma] =(double) atof(argv[i]);
             }
             else if((strcmp(str,"-gamma")) == 0)
             {
                 i++;
                 res[igamma] = atof(argv[i]);
             }

             else if((strcmp(str,"-mu"))==0)
             {
                 i++;
                 res[imu] = atof(argv[i]);
             }

             else if((strcmp(str, "-nbVilles"))==0)
             {
                 i++;
                 res[inbVilles]= atof(argv[i]);
             }
             else if((strcmp(str,"-topology"))==0)
             {
                 i++;
                 res[itopology] = atof(argv[i]);
             }
             else if((strcmp(str,"-rho"))==0){
                 i++;
                 res[irho] = atof(argv[i]);
             }
             else if((strcmp(str,"-epsilon"))==0){
                 i++;
                 res[iepsilon] = atof(argv[i]);
             }
             else if((strcmp(str,"-tmax"))==0){
                 i++;
                 res[itmax] = atof(argv[i]);
             }
             else if((strcmp(str,"-unitTIME"))==0){
                 i++;
                 res[iunitTIME] = atof(argv[i]);

             }
             else if((strcmp(str,"-seed"))==0){
                 i++;
                 res[iseed] = atof(argv[i]);
             }
             else if((strcmp(str,"-S0"))==0){
                 i++;
                 res[iS0] = atof(argv[i]);
             }
             else if((strcmp(str,"-E0"))==0){
                 i++;
                 res[iE0] = atof(argv[i]);
             }
             else if((strcmp(str,"-I0"))==0){
                 i++;
                 res[iI0] = atof(argv[i]);
             }
             else if((strcmp(str,"-R0"))==0){
                 i++;
                 res[iR0] = atof(argv[i]);
             }
             else if((strcmp(str,"-N0"))==0){
                 i++;
                 res[iN0] = atof(argv[i]);
             }
             else if((strcmp(str,"-beta0"))==0){
                 i++;
                 res[ibeta0] = atof(argv[i]);
             }
             else if((strcmp(str,"-beta1"))==0){
                 i++;
                 res[ibeta1] = atof(argv[i]);
             }
             else if((strcmp(str,"-phiMIN"))==0){
                 i++;
                 res[iphiMIN] = atof(argv[i]);
             }
             else if((strcmp(str,"-phiMAX"))==0){
                 i++;
                 res[iphiMAX] = atof(argv[i]);
             }
             else if((strcmp(str,"-nbSimu"))==0){
                 i++;
                 res[iNbSimu] = atof(argv[i]);
             }
             else if((strcmp(str,"-periodeBETA"))==0){
                 i++;
                 res[iPeriodeBETA] = atof(argv[i]);
             }
             else if((strcmp(str,"-typeRNG"))==0){
                 i++;
                 res[itypeRNG] = atof(argv[i]);
             }

        }

return res;

}

/******************************************************************************************************/
//Le but de cette fonction est d'initialiser les nombres de prosonnes dans chaque groupe (S, E, I, R, N) pour chaque ville
void initialerY_00(int nbVilles,int nvar, unsigned long **&y, unsigned long S0, unsigned long E0, unsigned long I0,
                                                        unsigned long R0, unsigned long N0){
    y = new unsigned long*[nbVilles];
    for (int nextville = 0; nextville <  nbVilles; nextville++)
    {
        y[nextville] = new unsigned long[nvar];
        y[nextville][iN] = N0;
        y[nextville][iS] = S0;
        y[nextville][iE] = E0;
        y[nextville][iI] = I0;
        y[nextville][iR] = R0;
    }
    return;
}

/******************************************************************************************************/
//Le but de cette fonction est de faire l'vnement m
void fairEvenementM_00(int m, vector<double> sigma, int nextville, unsigned long **y, double *matrCumulI){

    switch ( m ) {
     case 1: // vnement : une personne est ne
        y[nextville][iS]=y[nextville][iS]+1;
        y[nextville][iN]=y[nextville][iN]+1;
       break;

     case 2:  // vnement : une personne est morte
        y[nextville][iS]=y[nextville][iS]-1;
        y[nextville][iN]=y[nextville][iN]-1;
       break;

     case 3: // vnement : une personne expose est morte
        if(sigma[nextville]==INFINITY)   cout<<"invalid with m=3, when sigma==INFINITY "<<endl;
        else{
            y[nextville][iE]=y[nextville][iE]-1;
            y[nextville][iN]=y[nextville][iN]-1;
        }
       break;

     case 4: // vnement : une personne infecte est morte
        y[nextville][iI]=y[nextville][iI]-1;
        y[nextville][iN]=y[nextville][iN]-1;
       break;

     case 5:// vnement : une personne gurie est morte
        y[nextville][iR]=y[nextville][iR]-1;
        y[nextville][iN]=y[nextville][iN]-1;
      break;

     case 6:// vnement : une personne est expose S->E
        if(sigma[nextville]==INFINITY){
            y[nextville][iS]=y[nextville][iS]-1;
            y[nextville][iI]=y[nextville][iI]+1;
            matrCumulI[nextville]++;
        }
        else{
            y[nextville][iS]=y[nextville][iS]-1;
            y[nextville][iE]=y[nextville][iE]+1;
        }

      break;

     case 7:// vnement : une personne expose est infect E->I
        if(sigma[nextville]==INFINITY)  cout<<"invalid with m=7, when sigma==INFINITY "<<endl;
        else{
             y[nextville][iE]=y[nextville][iE]-1;
             y[nextville][iI]=y[nextville][iI]+1;
             matrCumulI[nextville]++;
        }
        break;

     case 8:   // vnement : une personne infect est gurie
        y[nextville][iI]=y[nextville][iI]-1;
        y[nextville][iR]=y[nextville][iR]+1;
        break;

     default:
        cout<<"prob dansfairEvenementM_00"<<endl;
        cout<<" m = "<<m<<endl;
       break;
     }
}

/******************************************************************************************************/
 /*
   Fonction est d'initialiser l'adresse pour un pointeur  1 dimention.
   */

 void creerAdresse1Point_00(int nbVilles,double *&Matrice){
     Matrice = new double[nbVilles];
     for(int i = 0; i< nbVilles ; i++) {
         Matrice[i] = 0.0;
     }

 }

 /******************************************************************************************************/
 /*
   Fonction est d'initialiser l'adresse pour un pointeur  2 dimention.
   */
void creerAdresse2Point_00(int nbVilles, int nevent, double **&Matrice){
    Matrice = new double*[nbVilles];
    for (int nextville = 0; nextville <  nbVilles; nextville++)
    {
        Matrice[nextville] = new double[nevent];
        for(int i = 0; i< nevent ; i++) Matrice[nextville][i] = 0.0;
    }

}
/******************************************************************************************************/
 /*
   Fonction cre un chiffre entre 0 et 1 au hasard.
   */
double rngFast(int *idum){
        return(ran1(idum,GENERATE, (char*)"a"));
    }

double rngGood(){
       return rand()/(double)RAND_MAX;
    }

/******************************************************************************************************/
/* Le but de cette fonction est de sauvegarder les resultats dans les fichier de text,
 un fichier par ville
 */
 void writeInFile_00(vector<double> valParSIM,vector<vector<vector<double> > > tabVille){

     //Creer les FILES pour sauvegarder les resultats
     // le nombre de villes utilis
     int nbVilles =(int) valParSIM[inbVilles];

     double S0, E0, I0, R0, N0;
     double t, S, E,P,R,I, N;

     char line0[255] ="time\t\tS\t\tE\t\tI\t\tR\t\tN\t\tP \n";
     char line[255] = "";
     int sizeTable = 0;

     for(int i=0; i<nbVilles; i++){
         char nameFilePAM[255]="";
         S0 = tabVille[i][0][iS];
         E0 = tabVille[i][0][iE];
         I0 = tabVille[i][0][iI];
         R0 = tabVille[i][0][iR];
         N0 = tabVille[i][0][iN];
         sprintf(nameFilePAM,"Output/Vil%d_S%.0fE%.0fI%.0fR%.0fN%.0f_NbVil%.0fTmax%.1fBeta0_%.3fBeta1_%.3fSigma%.3fGamma%.3fmu%.5fRho%.3fEpsilon%.5fTypeRNG%.0f_SIM.csv",
                                    i,S0,E0,I,R0,N0,valParSIM[inbVilles],valParSIM[itmax],valParSIM[ibeta0],valParSIM[ibeta1],valParSIM[isigma],valParSIM[igamma],valParSIM[imu],valParSIM[irho],valParSIM[iepsilon],valParSIM[itypeRNG]);
             FILE *pFile;
             pFile = fopen(nameFilePAM,"w+");
             fputs (line0,pFile);
             sizeTable = tabVille[i].size();
             for(int j=0; j<sizeTable;j++){
                 t = tabVille[i][j][it];
                 S = tabVille[i][j][iS];
                 E = tabVille[i][j][iE];
                 I = tabVille[i][j][iI];
                 R = tabVille[i][j][iR];                 
                 N = tabVille[i][j][iN];
                 P = tabVille[i][j][iP];
                 sprintf(line, "%.1f\t\t%.1f\t\t%.1f\t\t%.1f\t\t%.1f\t\t%.1f\t\t%.1f \n",t ,S, E, I,R, N,P);
                 fputs (line,pFile);
             }
             fclose (pFile);
     }
    return;
}

 /******************************************************************************************************/
 void writeInFileToTal_00(string nameFile,vector<vector<double> >  tableauToTal){
    char line[255];
    FILE * pFile;
    pFile = fopen (nameFile.c_str(),"w+");
    int size = tableauToTal.size();
    for(int i=0; i<size; i++){
        sprintf(line, "%.2f \t\t %.0f \t\t %.0f \t\t\t %.0f \t\t\t %.0f \t\t\t %.0f \n",tableauToTal[i][it], tableauToTal[i][iS], tableauToTal[i][iE],tableauToTal[i][iI],tableauToTal[i][iR],tableauToTal[i][iN],tableauToTal[i][iP]);
        fputs (line,pFile);
    }
    fclose (pFile);
    return;
}

 vector<vector<double> >  tableauToTal_00(vector<vector<vector<double> > > tabVille){
    vector<vector<double> > res;
    vector<double> tmpTotal;
    double temps=0.0,totalS=0, totalI=0, totalE=0, totalR=0, totalN=0, totalP=0;
    int nbVilles = tabVille.size();
    int size = tabVille[0].size();
    for(int j=0; j<size; j++){
        totalS=0, totalI=0, totalE=0, totalR=0, totalN=0;
        temps = tabVille[0][j][it];
        for(int i=0; i<nbVilles; i++){
            totalS += tabVille[i][j][iS];
            totalE += tabVille[i][j][iE];
            totalI += tabVille[i][j][iI];
            totalR += tabVille[i][j][iR];
            totalN += tabVille[i][j][iN];
            totalP += tabVille[i][j][iP];
        }
        tmpTotal.push_back(temps);tmpTotal.push_back(totalS); tmpTotal.push_back(totalE); tmpTotal.push_back(totalI); tmpTotal.push_back(totalR); tmpTotal.push_back(totalN); tmpTotal.push_back(totalP);
        res.push_back(tmpTotal);
        tmpTotal.clear();
    }
    return res;
}

/******************************************************************************************************/
 // Fonction  calculer le taux de contact beta_i qui subit un forcage saisonier, pour la population i.
 // beta_i(t) = beta0 *(1+ beta1*cos(2*pi*t/T + phi_i);
  vector<double>  calculerBeta_00(int nbVilles, vector<double> beta0, vector<double> beta1,
                       double periodeBETA,double t, vector<double> phi){
      vector<double> beta=vector<double>(nbVilles,0.0);
      for(int i=0; i<nbVilles; i++){
          beta[i] = beta0[i]*(1 + beta1[i] *cos(2*Pi*t/periodeBETA + phi[i]));
      }
      return(beta);
  }
  /******************************************************************************************************/
  // Fonction  calculer les phases pour chaque la population
  vector<double> calculerPhases_00(int nbVilles, double phiMAX, double phiMIN){
      vector<double> res = vector<double>(nbVilles,0);
      double unit_phi = 0.0;
      if(nbVilles==1) unit_phi = phiMIN;
      else
       unit_phi = (double)(phiMAX-phiMIN)/(nbVilles-1);
      double phi_tp=0.0;
      for(int i=0; i<nbVilles; i++){
          res[i] = phiMIN + phi_tp;
          phi_tp = phi_tp + unit_phi;
         // cout<<"phi"<<i<<" = "<<res[i]<<endl;
      }

      return res;
  }
  /******************************************************************************************************/
  // Fonction  calculer les valeurs de couplage pour entre deux population i et j
  double** calculerTauxCouplage_00(int nbVilles, double rho0){
      double **arr_rho;
      arr_rho = new double*[nbVilles];
      for (int i = 0; i<nbVilles; i++)
      {
          arr_rho[i] = new double[nbVilles];
          for(int j = 0; j< nbVilles ; j++) {
              if(j==i) arr_rho[i][j] = 1.0;
              else
              arr_rho[i][j] = rho0;
                //  cout<<"arr_rho"<<i<<j<<"   "<< arr_rho[i][j];
          }
          //cout<<endl;
      }
      return arr_rho;
  }

  /******************************************************************************************************/
  // Fonction  calculer les valeurs de epsilon (taux d'infection de l'extérieur)pour entre deux population i et j
  double **calculerTauxExterieur_00(int nbVilles, double epsilon0){
      double **epsilon;
      epsilon = new double*[nbVilles];
      for (int i = 0; i<nbVilles; i++)
      {
          epsilon[i] = new double[nbVilles];
          for(int j = 0; j< nbVilles ; j++) {
              if(j==i) epsilon[i][j] = 0.0;
              else
              epsilon[i][j] = epsilon0;
                //  cout<<"epsilon"<<i<<j<<"   "<< epsilon[i][j];
          }
         // cout<<endl;
      }
      return epsilon;
  }

  /******************************************************************************************************/
  //Le but de cette fonction est de calculer les valeurs de Lamda (c'est le taux de transmission S->E) pour chaque ville
  void calculerLamda_00(int nbVilles, double **epsilon, vector<double>beta, double **rho,
                     unsigned long **y, double *lamda){
      double total_ij = 0.0;
      double totalRHO = 0.0;

      for(int i=0; i<nbVilles; i++){
          total_ij = 0.0;
          totalRHO = 0.0;
          for(int j=0; j<nbVilles;j++){
              if(j!=i){
                  total_ij = total_ij + (double)(rho[i][j]*(((1-epsilon[i][j])*beta[i]*y[j][iN] + epsilon[i][j]*beta[j]*y[i][iN]))*y[j][iI])/(y[i][iN]*y[j][iN]);
                  totalRHO = totalRHO + rho[i][j];
              }
          }
          if(totalRHO < 1.0){
              lamda[i] = (double)(((1 - totalRHO)*(beta[i]*y[i][iI]/y[i][iN])) + (total_ij));
          }
          else
          {
              lamda[i] = -123456789.0;
              cout<<" totalRHO is more than 1."<<endl;
              cout<<" Review the function calculerLamda_00. "<<endl;
		cout<<"matrix of  epsilon "<<endl;
				for(int i=0; i<nbVilles; i++){
					for(int j=0; j<nbVilles; j++){
						cout<<"  "<< epsilon[i][j];
					}
				cout<<endl;
				 }
		cout<<"matrix of  coupling  by contact!"<<endl;
				for(int i=0; i<nbVilles; i++){
					for(int j=0; j<nbVilles; j++){
						cout<<"  "<< rho[i][j];
					}
				cout<<endl;
				 }
		cout<<"value of y"<<endl;
			for(int i=0; i<nbVilles; i++) {
				cout<<"S="<<y[i][iS]<<"  E="<<y[i][iE]<<"  I="<<y[i][iI]<<"  R="<<y[i][iR]<<"N="<<y[i][iS]<<endl;
				cout<<"beta"<<i<<" ="<<beta[i]<<endl;
			}

          }
      }
  }
  /******************************************************************************************************/
  //Le but de cette fonction est de calculer les valeurs de la fonction de propensit
  void calculerF_00(int nbVilles,vector<double> rmu, vector<double> sigma, vector<double> gamma,
                    double *lamda, unsigned long **y,double **f){
      for(int nextville=0; nextville< nbVilles;nextville++){
          if(sigma[nextville]==INFINITY){
              //Exposed death
              f[nextville][3] = 0.0;
              //Becoming infectious E -> I
              f[nextville][7] = 0.0;
          }
          else{
              //Exposed death
              f[nextville][3] = rmu[nextville]*y[nextville][iE];
              //Becoming infectious E -> I
              f[nextville][7] = sigma[nextville]*y[nextville][iE];
          }
          //Birth
          f[nextville][1] = rmu[nextville]*y[nextville][iN];

          //Susceptible death
          f[nextville][2] = rmu[nextville]*y[nextville][iS];

          //Exposed death
         // f[nextville][3] = rmu*y[nextville][iE];

           //Infected death
          f[nextville][4] = rmu[nextville]*y[nextville][iI];

          //Recovered death
          f[nextville][5] = rmu[nextville]*y[nextville][iR];
	
          //Infection S -> E
          f[nextville][6] = lamda[nextville]*y[nextville][iS];

          //Becoming infectious E -> I
         // f[nextville][7] = sigma*y[nextville][iE];

          //recovery I -> R
          f[nextville][8] = gamma[nextville]*y[nextville][iI];
      }
  }
/******************************************************************************************************/
  // initialerY_01 : initialise le nombre d'individus dans chaque groupe S, E, I, R, N pour chaque sous-population.
  // Chaque sous-population est différente.
  // nbVilles: nombre de sous-populations
  // valeursSEIR0 : c'est un vecteur qui a des noeuds.
  // Le nombre de noeuds est égal au nombre de sous-populations.
  // Chaque noeud contient les valeurs initiales des variables S, E, I, R, N pour chaque sous-population.

  void initialerY_01(int nbVilles,int nvar, unsigned long **&y, vector<vector<double> > villeSEIR0, 
			unsigned long S0, unsigned long E0, unsigned long I0, unsigned long R0){
      y = new unsigned long*[nbVilles];
      int szVilleSEIR0 = villeSEIR0.size();
    //  cout<<"szVilleSEIR0  = "<<szVilleSEIR0<<endl;
      for (int i = 0; i <  nbVilles; i++)
      {
           y[i] = new unsigned long[nvar];
           if(szVilleSEIR0 == 0){
               y[i][iI] = I0;
               y[i][iN] = S0+E0+I0+R0;
               y[i][iS] = S0;
               y[i][iE] = E0;
               y[i][iR] = R0;

           }
           else{
               bool vilExist = false;
               int ixVil=-1;
               for(int j=0; j<szVilleSEIR0; j++){
                   if(i==villeSEIR0[j][0]){
                       vilExist = true;
                       ixVil = j;
                       break;
                   }
               }
               if(vilExist){
                   y[i][iS] = villeSEIR0[ixVil][iS];
                   y[i][iE] = villeSEIR0[ixVil][iE];
                   y[i][iI] = villeSEIR0[ixVil][iI];
                   y[i][iR] = villeSEIR0[ixVil][iR];
                   y[i][iN] = y[i][iS]+y[i][iE]+y[i][iI]+y[i][iR];
               }
               else{
                   y[i][iI] = I0;
                   y[i][iN] = S0+E0+I0+R0;
                   y[i][iS] = S0;
                   y[i][iE] = E0;
                   y[i][iR] = R0;
               }

          }
      }
      return;
  }
  //
  double sumProb(int nbVilles, int nevent, double **f){
      double fsumVilles = 0.0;
      for(int i = 0; i<nbVilles; i++){
          for(int j = 0; j<nevent; j++){
              fsumVilles += f[i][j];
          }
      }
      return(fsumVilles);
  }

  //
  int choisirVille(int nbVilles, int nevent, double p2, double **f)
    {
        //Caculer le total fsum de toutes les villes
        double fsumVilles = sumProb(nbVilles,nevent,f);
        double fville = 0.0, xpro=0.0;
        //Selectionne aleatoirement une ville ou se produira l'evenement
        if(p2==1.0) return(nbVilles-1);
        else
          if(p2==0.0)return(0);
        else
              if((0.0<p2)&&(p2<1.0)){
              for(int ville=0; ville<nbVilles; ville++){
                  for(int i=0; i< nevent;i++){
                      fville+=f[ville][i];
                  }
                  xpro = fville/fsumVilles;
                 // cout<<"ville="<<ville<<"   p2="<<p2<<"  xpro="<<xpro<<endl;
                  if(p2 <= xpro) return(ville);
                  }
         }
        else return(-1);
    }

  //
  int choisirEvenement(int nevent, int nextville, double p3, double **f){
      double fville=0.0, xpro=0.0, fprop=0.0;
      for(int i=0; i< nevent;i++){
          fville+=f[nextville][i];
      }

      if(p3==1.0) return(nevent-1);
      else
        if(p3==0.0)return(1);
      else
            if((0.0<p3)&&(p3<1.0)){
                for(int j=0; j< nevent; j++){
                    fprop += f[nextville][j];
                    xpro = (double)(fprop/fville);
                   // cout<<"p3="<<p3<<"  fprop="<<fprop<<"  xpro="<<xpro<<endl;
                    if(p3 <= xpro) {
                        return(j);
                    }
                }
       }
      else return(-1);

  }
//
//resize un vector with length given
vector<double> resizeVector1D(vector<double> V,int n)
{
    if(V.empty()) return(V);
    else{
        int legV=V.size();
        int size=ceil(n/legV);
        if(n > legV){
            vector<double> res;
            if((n%2)==0){
                for(int i=0; i<legV; i++)
                    for(int j=0; j<size; j++)
                           res.push_back(V[i]);
                res.resize(n);
            }
            else{
                for(int i=0; i<n; i++)
                    for(int j=0; j<legV; j++)
                        res.push_back(V[j]);
                res.resize(n);
            }
            return(res);
        }
        else if(n<legV){
            V.resize(n);
            return(V);
        }
        else return(V);
    }
}

  /***************************************************************************************************************/
// C'est la fonction principale simulerSEIR avec beta constant  la fois beta sinusoidal
//
  vector<vector<vector<double> > > simulerSEIR(int nbVilles,vector<double>sigma,vector<double>gamma,vector<double>mu,
                                               vector<double>beta0,vector<double>beta1,vector<double>phi,
                                               double** arr_rho, double** epsilon,
                                               double seed,double unitTIME, double tmax, double typeRNG,
                                               double T, unsigned long S,  unsigned long E,  unsigned long I,  unsigned long R,
                                               vector<vector<double> >valeursSEIR0)
{
    // parametres pour SIMULATIONvector<double>rho, vector<double>epsilon,
    // Initialer le nombre de cas = 5 (S, E, I, R, N); le nombre d'evenements = 9
    int nevent = 9, nbVarSEIRN = 6;
    
    unsigned long S0 = S, E0 =E, I0 =I, R0 = R;
    
    double p1 = 0.0, p2 = 0.0, p3 = 0.0;
    double t = 0.0, tstep = 0.0;
    double fsumVilles = 0.0;
    int event = 0, nextville = 0;
    //Creer une variable de temps
    double pointTemps = unitTIME;

    // La programme marche :
    //Initialisation
    double *lamda = new double[nbVilles];
    double *matrCumulI = new double[nbVilles];
    for(int i=0; i<nbVilles; i++) matrCumulI[i]=0.0;
    vector<double> beta_sinusoidal = vector<double>(nbVilles,0);
    
    //Sauvegarder les valeurs dans chaque cas et dans chaque evenement
    unsigned long **y;
    double **f;
    f = new double*[nbVilles];
    for (int nextville = 0; nextville <  nbVilles; nextville++)
    {
        f[nextville] = new double[nevent];
        for(int i = 0; i< nevent ; i++) f[nextville][i] = 0.0;
    }

    //Initialiser les valeurs des variables
    initialerY_01(nbVilles,nbVarSEIRN,y,valeursSEIR0,S0,E0,I0,R0);
    /*
    for(int i=0; i<nbVilles; i++){
        cout<<"ville = "<<i<<" S0="<<y[i][iS]<<" E0="<<y[i][iE]<<" I0="<<y[i][iI]<<" R0="<<y[i][iR]<<" N0="<<y[i][iN]<<endl;
        cout<<"beta0="<<beta0[i]<<endl;
         cout<<"beta1="<<beta1[i]<<endl;
          cout<<"sigma="<<sigma[i]<<endl;
           cout<<"gamma="<<gamma[i]<<endl;
           cout<<"phi="<<phi[i]<<endl;
    }
    cout<<"seed="<<seed<<endl;
    */


    //Creer les tables qui sauvegardent des resultats quand le programme marche
    vector<double> tpmSEIRNP;
    vector<vector<vector<double> > > tableResVilles(nbVilles);
    for(int i=0; i<nbVilles; i++){
        tpmSEIRNP.push_back(0.0);tpmSEIRNP.push_back(y[i][iS]); tpmSEIRNP.push_back(y[i][iE]); tpmSEIRNP.push_back(y[i][iI]);tpmSEIRNP.push_back(y[i][iR]);tpmSEIRNP.push_back(y[i][iN]);tpmSEIRNP.push_back(0);
        tableResVilles[i].push_back(tpmSEIRNP);
        tpmSEIRNP.clear();
    }
    double xpro = 0.0, fville = 0.0, fprop =0.0;

    //PROGRAMME COMMENCE.......................................................
    //Initialiser les nombres de personnes de chaque groupe

    //Calculer le temps o la programme marche selon le temps de CPU
    clock_t start, end;
    start = clock();
    static int seed0 = seed;
    int  *idum = &seed0;
    srand(seed);
    //Caculer le temps o la raction se produit
    while(t<tmax){

        //Initialer les valeurs de beta pour chaque ville
        beta_sinusoidal=calculerBeta_00(nbVilles,beta0,beta1,T,t,phi);
        //Calculer les valeurs de Lamda
        calculerLamda_00(nbVilles,epsilon,beta_sinusoidal,arr_rho,y,lamda);
        //Calculer les valleurs de la foction de propensit
        calculerF_00(nbVilles,mu,sigma,gamma,lamda,y,f);

//        cout<<"t"<<"  ="<<t<<endl;
//        for (int i = 0; i <  nbVilles; i++)
//        {
//            cout<<"ville="<<i<<endl;
//            cout<<" S="<<y[i][iS]<<" E="<<y[i][iE]<<" I="<<y[i][iI]<<" R="<<y[i][iR]<<" N"<<y[i][iN]<<endl;
//            cout<<" beta="<<beta_sinusoidal[i]<<endl;
//            cout<<" lamda="<<lamda[i]<<endl;
//        }
//        for (int i = 0; i <  nbVilles; i++)
//        {
//            for(int j = 0; j< nevent ; j++) cout<<"     f"<<i<<j<<"  ="<<f[i][j]<<endl;
//        }
/*
	 for(int i=0; i<nbVilles; i++)                
                 for(int j=0; j<nevent; j++)
                     cout<<"f"<<i<<j<<"= "<<f[i][j]<<endl;
*/

        //Generate uniform random numbers
        if(typeRNG==FAST) {
            p1 = ran1(idum,GENERATE, (char*)"a");
            p2 = ran1(idum,GENERATE, (char*)"a");
            p3 = ran1(idum,GENERATE, (char*)"a");
        }
        else
            if(typeRNG==GOOD){
                p1 = rand()/(double)RAND_MAX;
                p2 = rand()/(double)RAND_MAX;
                p3 = rand()/(double)RAND_MAX;
            }

      //  cout<<"p1 = "<<p1<<"  p2 = "<<p2 <<"    p3 ="<<p3<<endl;

        //Caculer le total fsum de toutes les villes
        fsumVilles = 0.0;
        for(int i = 0; i<nbVilles; i++){
            for(int j = 0; j<nevent; j++){
                fsumVilles += f[i][j];
            }
        }
       //Determine time interval and update time
        if(fsumVilles > 0.0) {
            tstep = (double)(-log(p1)/fsumVilles);
            t = t + tstep;
        }
        else {
            cout<<"no event"<<endl;
            cout<<"parce que: fsumVilles = "<<fsumVilles<<endl;
            break;
        }
        //Selectionne aleatoirement une ville ou se produira l'evenement
        nextville = choisirVille(nbVilles,nevent,p2,f);
        if(nextville != -1){
            // Selectionne aleatoirement la nature du prochain evenement dans nextville
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
        if(event != -1){
            //Faire l'evenement m
            fairEvenementM_00(event,sigma,nextville,y,matrCumulI);

            //Sauvegarder la valeur de t, le nombre de personnes dans le groupe I et l'vnement qui se produit
            if(t> pointTemps){
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
       //Le programme fini!
    }
    end = clock();
    // Le temps CPU utilis:
    double totalDure = (double)(end - start)/CLOCKS_PER_SEC;
   // printf("Le temps CPU utilisé est %.20f secondes avec typeRNG == %d.\n",totalDure,typeRNG);

    //Librer la mémoire
    // Librer la memoire des pointeurs à 1 dimension
    delete []lamda; //delete []beta_sinusoidal; delete []phi;
    // Librer la memoire des pointeurs à deux dimension
    for (int nextville = 0; nextville <  nbVilles; nextville++)
    {
        delete []y[nextville];
        delete []f[nextville];
    }
    delete []y; delete[]f;
   // cout<<"Le programme SimulerSEIR est avec succes!"<<endl;

    return(tableResVilles);
    //The end the programme principal
}
  /////////////////////////////////////////////////////////////////

  vector<vector<vector<double> > > resSimSEIR(vector<double> valParSIM,vector<vector<double> >valSEIR0)
  {

      // le nombre de villes utilis
      int nbVilles=(int) valParSIM[inbVilles];

      //Taux de E->I
      vector<double> tpsigma; tpsigma.push_back(valParSIM[isigma]);
      vector<double> sigma=resizeVector1D(tpsigma,nbVilles);

      //Taux de I->R
       vector<double> tpgamma; tpgamma.push_back(valParSIM[igamma]);
       vector<double> gamma=resizeVector1D(tpgamma,nbVilles);

      //Taux de natalite
     // double mu = valParSIM[imu];
       vector<double> tpmu; tpmu.push_back(valParSIM[imu]);
       vector<double> mu=resizeVector1D(tpmu,nbVilles);

       //  double beta0 = valParSIM[ibeta0];
       vector<double> tpbeta0; tpbeta0.push_back(valParSIM[ibeta0]);
       vector<double> beta0=resizeVector1D(tpbeta0,nbVilles);

       vector<double> tpbeta1; tpbeta1.push_back(valParSIM[ibeta1]);
       vector<double> beta1=resizeVector1D(tpbeta1,nbVilles);


       double phiMIN = valParSIM[iphiMIN];
       double phiMAX = valParSIM[iphiMAX];
       vector<double> phi = calculerPhases_00(nbVilles,phiMAX,phiMIN);

      //Taux de couplage entre deux villes
      double rho = valParSIM[irho];
      double **arr_rho = calculerTauxCouplage_00(nbVilles,rho);

      //Proportion des infections de l'extrieur
      double epsilon0 = valParSIM[iepsilon];
      double **epsilon = calculerTauxExterieur_00(nbVilles,epsilon0);

      //int nbSimu = valParSIM[iNbSimu];

      double unitTIME = valParSIM[iunitTIME];
      double seed = valParSIM[iseed];
      double tmax = valParSIM[itmax];
      double T = valParSIM[iPeriodeBETA];
      int typeRNG = (int) valParSIM[itypeRNG];
      unsigned long S0 = valParSIM[iS0], E0 =valParSIM[iE0], I0 =valParSIM[iI0], R0 = valParSIM[iR0];

      vector<vector<vector<double> > >  res = simulerSEIR(nbVilles,sigma,gamma,mu,beta0,beta1,phi,
                                             arr_rho,epsilon,seed,unitTIME,tmax,typeRNG,T,S0,E0,I0,R0,valSEIR0);

      for (int nextville = 0; nextville <  nbVilles; nextville++)
      {
          delete []arr_rho[nextville];
          delete []epsilon[nextville];
      }

      return(res);


  }


//____________________________________________________________________________________________________
//------------rate functions(tinh toan ti so de transitions)--------------------------------------------
extern "C" {
   // double* Cratefunc(double* x,double* params, double m_T){
    double* seirratefunc(double* x,int nbVilles,vector<double>beta0,vector<double> beta1,
                         vector<double> mu, vector<double> sigma, vector<double> gamma,vector<double>phi,
                         double** arr_rho, double** epsilon, double T, double m_T){

        int nevent = 9, nbVarSEIRN = 6;
        //Initialisation
        double *lamda = new double[nbVilles];
        vector<double> beta_sinusoidal = vector<double>(nbVilles,0.0);
        double **f;
        f = new double*[nbVilles];	
        for (int i = 0; i <  nbVilles; i++)
        {
            f[i] = new double[nevent];	    
            for(int j = 0; j< nevent ; j++) f[i][j] = 0.0;
        }
        //chuyen doi vector x 1 chieu thanh 2 chieu **y
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
        //Initialer les valeurs de beta pour chaque ville
        beta_sinusoidal=calculerBeta_00(nbVilles,beta0,beta1,T,m_T,phi);
        //Calculer les valeurs de Lamda
        calculerLamda_00(nbVilles,epsilon,beta_sinusoidal,arr_rho,y,lamda);
        //Calculer les valleurs de la foction de propensit
        calculerF_00(nbVilles,mu,sigma,gamma,lamda,y,f);
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
       //Librer la mémoire
       // Librer la memoire des pointeurs à 1 dimension
        delete []lamda;
        // Librer la memoire des pointeurs à deux dimension
        for (int i = 0; i <  nbVilles; i++)
        {
            delete []y[i];
            delete []f[i];	 
        }
      delete []y; delete[]f; 
      return(resRates);
    }
}
}
/******************************************************************************************************/
