
#include "batchVac_01.h"

#include <float.h>
//generic varibale
int nbSuscVac = 0; //somme de personnes susceptibles vaccinees
int nbPersInf = 0; //somme de personnes susceptibles infectees
int nbParSIM = 21;


/******************************************************************************************************/
//Initialiser les valeurs initiales pour les parametres de SIM
double *initialervalParSIM(int argc, char** argv){
    double *res = new double[nbParSIM];
    // temps maximum de simulation
    res[itmax] = 5*365;

    // nombre de sous-populations dans une metapopulation
     res[inbVilles] = 1.0;

    // taux d'exposition E -> I
    res[isigma] = (double)1/8;

    // taux d'exposition I -> R
    res[igamma] = (double)1/5;

    // taux de mortalite et natalite
    res[imu] = (double)1/(70*365);
    // taux d'infection de l'extrieur
    res[iepsilon] = 0.0;

    // taux de contact entre deux population
    res[irho] = 0.01;

    // pas de temps auquel les valeurs des variables sont enregistres (indexTEMPS)
    res[iunitTIME] = 1.0;

    // graine du gnrateur de nombres alatoires.
    res[iseed] = 10;

    // valeur initiale de susceptibles dans une sous-population.
    res[iS0] = 7415;   res[iE0] = 27; res[iI0] = 16; res[iR0] = 92539;
    res[iN0] =  res[iS0]+ res[iE0]+ res[iI0]+ res[iR0];

    // valeur moyenne du taux de contact (/ind/indexTEMPS).
     res[ibeta0] = (double)1000/365;

    // beta1: phase du forage
     res[ibeta1] = 0.0;

    // phiMIN: phase minimum dans la mtapopulation.
    res[iphiMIN] = 0.0;

    // phiMAX: phase maximun dans la mtapopulation.
     res[iphiMAX] = 0.0;

    // nombre de fois de simulation
     res[inbSimu] = 5;

    // periode dans la fionction de BETA


    // periode dans la fionction de BETA
     res[itypeRNG] = 0;

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
                 res[inbSimu] = atof(argv[i]);
             }
             else if((strcmp(str,"-PeriBETA"))==0){
                 i++;
                 res[iPeriBETA] = atof(argv[i]);
             }
             else if((strcmp(str,"-typeRNG"))==0){
                 i++;
                 res[itypeRNG] = atof(argv[i]);
             }

        }

return res;

}

/******************************************************************************************************/
/*
Faire l'événement m dans la sous-population "nextville".
        INPUT:
                m: événement m
                nextville: sous-population où l'événement m se produit.
                y: table qui enregistre les nombres d'individus dans chaque groupe.
        sigma: c'est la condition pour le modèle SEIR ou le SIR
        matrCumulI: c'est la matrice qui contient le nombre d'infectés pendant une unité de temps
        OUTPUT:
                table qui enregistre les nombres d'individus dans chaque groupe après que l'événement m a fini.
  */
void fairEvenementM(int m, double sigma, int nextville, unsigned long **y, int *matrCumulI){

    switch ( m ) {
     case 1: // événement : une personne est néé
        y[nextville][iS]=y[nextville][iS]+1;
        y[nextville][iN]=y[nextville][iN]+1;
       break;

     case 2:  // événement : une personne est morte
        y[nextville][iS]=y[nextville][iS]-1;
        y[nextville][iN]=y[nextville][iN]-1;
       break;

     case 3: // événement : une personne exposée est morte
        if(sigma==INFINITY)   cout<<"invalid with m=3, when sigma==INFINITY "<<endl;// modèle SIR
        else{//modèle SEIR
            y[nextville][iE]=y[nextville][iE]-1;
            y[nextville][iN]=y[nextville][iN]-1;
        }
       break;

     case 4: // événement : une personne infectée est morte
        y[nextville][iI]=y[nextville][iI]-1;
        y[nextville][iN]=y[nextville][iN]-1;
       break;

     case 5:// événement : une personne guriée est morte
        y[nextville][iR]=y[nextville][iR]-1;
        y[nextville][iN]=y[nextville][iN]-1;
      break;

     case 6:// événement : une personne est exposée S->E
        if(sigma==INFINITY){//modèle SIR
            y[nextville][iS]=y[nextville][iS]-1;
            y[nextville][iI]=y[nextville][iI]+1;
            matrCumulI[nextville]++;
        }
        else{//modèle SEIR
            y[nextville][iS]=y[nextville][iS]-1;
            y[nextville][iE]=y[nextville][iE]+1;
            nbPersInf+=1;
        }

      break;

     case 7:// événement : une personne exposée est infectée E->I
        if(sigma==INFINITY)  cout<<"invalid with m=7, when sigma==INFINITY "<<endl;//modèle SIR
        else{//modèle SEIR
             y[nextville][iE]=y[nextville][iE]-1;
             y[nextville][iI]=y[nextville][iI]+1;
             matrCumulI[nextville]++;
        }
        break;

     case 8: // événement : une personne infecée est guriée
        y[nextville][iI]=y[nextville][iI]-1;
        y[nextville][iR]=y[nextville][iR]+1;
        break;

     default:
        cout<<"prob dans fairEvenementM_00"<<endl;
        cout<<" m = "<<m<<endl;
       break;
     }
}

/******************************************************************************************************/
unsigned long **initialerVarSEIR(int nbVilles,int nvar,
                                 unsigned long S0, unsigned long E0, unsigned long I0, unsigned long R0){
     unsigned long ** y;
     y = new unsigned long*[nbVilles];

     for (int i = 0; i <  nbVilles; i++)
     {
          y[i] = new unsigned long[nvar];
          y[i][iI] = I0;
          y[i][iN] = S0+E0+I0+R0;
          y[i][iS] = S0;
          y[i][iE] = E0;
          y[i][iR] = R0;
     }
     return y;
 }
/******************************************************************************************************/
 /*
   Fonction est d'initialiser l'adresse pour un pointeur  1 dimention.
   */

 void create1D(int nbVilles,double *&vect){
     vect = new double[nbVilles];
     for(int i = 0; i< nbVilles ; i++) {
         vect[i] = 0.0;
     }

 }

 /******************************************************************************************************/
 /*
   Fonction est d'initialiser l'adresse pour un pointeur  2 dimention.
   */
void create2D(int nbVilles, int nevent, double **&vect){
    vect = new double*[nbVilles];
    for (int nextville = 0; nextville <  nbVilles; nextville++)
    {
        vect[nextville] = new double[nevent];
        for(int i = 0; i< nevent ; i++) vect[nextville][i] = 0.0;
    }

}
/******************************************************************************************************/
 /*
   Fonction creer un chiffre entre 0 et 1 au hasard.
   */

double rngGood(){
       return rand()/(double)RAND_MAX;
    }

/******************************************************************************************************/
/*
Objectif:
    Stocker les resultats dans les fichier de text, un fichier par sous-population.
Entrée:
    valParSIM: vecteur qui contient les valeurs de parameters.
    tabVille: table qui stocke les résultats des villes après avoir fait la simulation.
Sortie:
    Fichiers qui retiennent les résulats des villes. Leurs noms viennent des noms des paramèters.
*/
 void writeInFile(double*valParSIM,vector<vector<vector< unsigned long> > > tabVille){

     //Creer les FILES pour sauvegarder les resultats
     // le nombre de villes utilis
     int nbVilles =(int) valParSIM[inbVilles];

     unsigned long S0, E0, I0, R0, N0;
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
         sprintf(nameFilePAM,"Output/Vil%d_S%luE%luI%luR%luN%lu_NbVil%.0fTmax%.1fBeta0_%.3fBeta1_%.3fSigma%.3fGamma%.3fmu%.5fRho%.3fEpsilon%.5fTypeRNG%.0f_SIM.csv",
                                    i,S0,E0,I0,R0,N0,valParSIM[inbVilles],valParSIM[itmax],valParSIM[ibeta0],valParSIM[ibeta1],valParSIM[isigma],valParSIM[igamma],valParSIM[imu],valParSIM[irho],valParSIM[iepsilon],valParSIM[itypeRNG]);
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
                 //P = tabVille[i][j][iP];
                 sprintf(line, "%.1f\t\t%.1f\t\t%.1f\t\t%.1f\t\t%.1f\t\t%.1f\t\t%.1f \n",t ,S, E, I,R, N);
                 fputs (line,pFile);
             }
             fclose (pFile);
     }
    return;
}

/******************************************************************************************************/
/*
Objectif:
    Créer une table à deux dimensions qui stocke le résultat total des toutes villes selon unité de temps.
    Après:
        enregistrer le résultat total à un fichier en text.
*/
 void writeInFileToTal_00(string nameFile,vector<vector<double> >  tableauToTal){
    char line[255];
    FILE * pFile;
    pFile = fopen (nameFile.c_str(),"w+");
    int size = tableauToTal.size();
    for(int i=0; i<size; i++){
        sprintf(line, "%.2f \t\t %.0f \t\t %.0f \t\t\t %.0f \t\t\t %.0f \t\t\t %.0f \t\t\t %.0f\n",
                tableauToTal[i][it], tableauToTal[i][iS], tableauToTal[i][iE],tableauToTal[i][iI],tableauToTal[i][iR],tableauToTal[i][iN]);
        fputs (line,pFile);
    }
    fclose (pFile);
    return;
}

 /************************/
 /*
 Objectif:
     Fichier total.
 */
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
            //totalP += tabVille[i][j][iP];
        }
        tmpTotal.push_back(temps);tmpTotal.push_back(totalS); tmpTotal.push_back(totalE); tmpTotal.push_back(totalI); tmpTotal.push_back(totalR); tmpTotal.push_back(totalN); tmpTotal.push_back(totalP);
        res.push_back(tmpTotal);
        tmpTotal.clear();
    }
    return res;
}

/******************************************************************************************************/
/*
Calcule le taux de contact beta qui suit un forçage saisonier, pour toutes les sous-populations
 avec beta_i(t) = beta0 *(1+ beta1*cos(2*pi*t + phi_i).
 INPUT :
    nbVilles: nombre de subpopulations.
        beta0 : valeur moyenne du taux de contact (/ind/jours).
        beta1 : l'amplitude du taux de contact (/ind/jours).
        phi_i : le phase du forçage pour la sous-population i, [0, 2*pi].
    t: au moment t de simulation.
 OUTPUT :
        un vecteur en 1 dimension contient les noeuds. Le nombre de noeuds est égal au nombre du sous-populations. La valeur de chaque noeud est la valeur de beta pour chaque sous-population.
 */
  double *calculerBeta(int nbVilles, double beta0, double beta1,
                       double periBETA,double t, double*phi){
      double *beta; create1D(nbVilles,beta);
      double coef=2*Pi*t/periBETA;
      for(int i=0; i<nbVilles; i++){
          beta[i] = beta0*(1 + beta1 *cos(coef + phi[i]));
      }
/*
      char line[255];
      FILE * pFile;
      pFile = fopen ("beta.txt","a+");

      sprintf(line, "%.5f \t\t %.5f \t\t %.5f \n",t,beta[0],beta[1]);
      fputs (line,pFile);

      fclose (pFile);
*/

       return(beta);
  }
  /******************************************************************************************************/
 /*
 Objectif:
    Calculer le phase pour chaque la population avec phi_i = i*(phiMAX - phiMIN)/(nbVille-1), avec i de 0 à
(nbVilles-1)

 Entrée :
        nbVille : nombre de sous-populations, chiffre entier.
        phiMAX : phase maximun,[0,2*pi].
        phiMIN : phase minimum, [0, 2*pi].

 Sortie :
       Un vecteur en 1 dimension contient les noeuds.
       Le nombre de noeuds est gal au nombre du sous-populations. La valeur de chaque noeud est le phase pour 				chaque sous-population.
  */
  double *calculerPhases(int nbVilles, double phiMAX, double phiMIN){
      double *res; create1D(nbVilles,res);
      double unit_phi = 0.0;
      if(nbVilles==1) unit_phi = phiMIN;
      else
       unit_phi = (double)(phiMAX-phiMIN)/(nbVilles-1);
      double phi_tp=0.0;
      for(int i=0; i<nbVilles; i++){
          res[i] = phiMIN + phi_tp;
          phi_tp = phi_tp + unit_phi;
      }

      return res;
  }
  /******************************************************************************************************/
 /*
Object:
    Fonction crée un table en deux dimensions qui contient les taux de couplage entre ville i et j,
    qui montre la probabilite pour que un ind x de la ville i visite la ville j.

Entrée:
    rho0: taux de couplage initiale

*/
  double** calculerProbVisiter(int nbVilles, double rho0){
      double **arr_rho;
      arr_rho = new double*[nbVilles];
      for (int i = 0; i<nbVilles; i++)
      {
          arr_rho[i] = new double[nbVilles];
          for(int j = 0; j< nbVilles ; j++) {
              double totaltp= (double)((nbVilles-1)*rho0);
              //\rho_ii
              if(j==i) arr_rho[i][j] = 1.0 -totaltp;
              else
              arr_rho[i][j] = rho0;
               //   cout<<"arr_rho"<<i<<j<<" =   "<< arr_rho[i][j];
          }
       //  cout<<endl;
      }
      return arr_rho;
  }

  /******************************************************************************************************/
  /* Fonction  calculer les valeurs de \xi
    \xi qui montre la probabilite
    pour que un indi y rencontre un indi x dans ville j qui vient de ville i
  */
  void calculerProbVoir(int nbVilles, unsigned long **y, double **arr_rho, double**arr_xi){
      //double **arr_xi; create2D(nbVilles,nbVilles,arr_xi);
      for (int i = 0; i<nbVilles; i++)
      {
         //calculer le demominateur
          double denom=0.0;
          for(int k=0; k<nbVilles; k++)
              denom = denom + y[k][iN]*arr_rho[k][i];
          //
          for(int j = 0; j< nbVilles ; j++) {
              arr_xi[i][j]= y[j][iN]*arr_rho[j][i]/denom;
          }
      }
  }

  /******************************************************************************************************/
/*Calculer les valeurs de Lamda (c'est le taux d'infection S -> E) pour toutes les sous-populations avec beta sinusoidal, le nouveau modle SEIR au moment t.
  INPUT:
          nbVilles: nombre de sous-populations, chiffre entier positif.
          beta: matrix en 1 dimension qui contient les valeurs de contact des toutes les sous-populations, (/ind/jours).
          arr_rho: matrice en 2 dimension, elle est la probabilite pour que un individu x visite une autre city j, [0,1]
          y: table servant à sauvegarder le nombre de groupe S, E, I, R au moment t.
          lamda: table pour sauvegarder les valeurs lamda au moment t, (/ind/jours).
          arr_xi: prob pour que un individu x de city i rencontre un indivudu y de city j.

  OUTPUT:
          valeus de taux d'infection des toutes les sous-populations.
  */
  void calculerLamda(int nbVilles, double **arr_xi, double *beta, double **arr_rho,
                     unsigned long **y, double *lamda){
      double part1 = 0.0, part2=0.0;

      for(int i=0; i<nbVilles; i++){
           double subpart1 = 0.0;
           part2 = 0.0;
          for(int j=0; j<nbVilles; j++){
              double rateIN = (double)y[j][iI]/y[j][iN];
                        //calculer part1
              subpart1 = subpart1 + arr_xi[i][j]*beta[j]*rateIN;
              //
                        //calculer part2
              if(j!=i)
                  part2 = part2 + arr_rho[i][j]*beta[j]*rateIN;
          }
          part1=arr_rho[i][i]*subpart1;

         lamda[i] = (double)(part1 + part2);
      }
  }
  /******************************************************************************************************/
/*
calculer les valeurs de la fonction de propensité, dans ce cas, beta est sinusoidal.
        INPUT:
                nbVill: nombre de sous-populations, chiffre entier.
                rmu: taux de naissance et taux de mortalité. (/jours).
                sigma: taux de transmission E -> I (/jours).
                gamma: taux de transmission I -> R (/jours).
                lamda: table pour sauvegarder les valeur de fonction de propensité (/ind/jours).
                y:  table pour  sauvegarder les nombres de personnes dans chaque groupe.
                f: table pour sauvegarder les valeurs de la fonction de propensité.

        OUTPUT:
                valeurs de fonctions de propensités des toutes les sous-populations.
  */
  void calculerF(int nbVilles,double rmu, double sigma, double gamma,
                    double *lamda, unsigned long **y,double **f){
      for(int nextville=0; nextville< nbVilles;nextville++){
     //cas:
          if(sigma==INFINITY){//modèle SIR
              //Exposed death
              f[nextville][3] = 0.0;
              //Becoming infectious E -> I
              f[nextville][7] = 0.0;
          }
          else{//modèle SEIR
              //Exposed death
              f[nextville][3] = rmu*y[nextville][iE];
              //Becoming infectious E -> I
              f[nextville][7] = sigma*y[nextville][iE];
          }

          //Birth
          f[nextville][1] = rmu*y[nextville][iN];

          //Susceptible death
          f[nextville][2] = rmu*y[nextville][iS];

           //Infected death
          f[nextville][4] = rmu*y[nextville][iI];

          //Recovered death
          f[nextville][5] = rmu*y[nextville][iR];

          //Infection S -> E
          f[nextville][6] = lamda[nextville]*y[nextville][iS];

          //recovery I -> R
          f[nextville][8] = gamma*y[nextville][iI];
      }
  }
/******************************************************************************************************/
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
/******************/
  /*
Objectif:
    Fonction choisit une ville où un événement va se produire.
Entrée:
    nbVille : nombre de sous-populations
        nevent : nombre d'événement qui peut se produire.
        p2 : probabilite de propensité
        f : fonctions de propensité
Sortie:
       ville choisie.
   */
  int choisirVille(int nbVilles, int nevent, double p2, double **f){
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
/******************/
   /*
Objectif:
    Fonction choisit l'événement qui se produit dans la sous-population choisie "nextville"
Entrée:
       nevent : nombre d'événement qui peut se produire dans cette sous-population.
       nextville: sous-population choisie
       p3 : probabilité s'établit pour qu'un événement peut se produire.
       f : fonctions de pensilite de la souspopulation choisie.
Sortie:
       événement choisi.
   */
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
  /******************/
/*
Objectif:
    Redimensionner un vecteur à une dimension avec la longueur donnée.
Entrée:
    V: vecteur redimensionnée
    n: longueur à la quelle le vecteur V doit arriver
Sortie:
    Vecteur à une dimension

*/
vector<double> resizeVector1D(vector<double> V,int n){
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
/******************************************************************************************************/
// Modifier a partir de la version vaccinerToutesVilles_V_1
int vaccinerToutesVilles_V_Nhan(double tauxglobVAC, int nbVilles, unsigned long **y, double_vector *politique)
{
    unsigned long nbPersonneVaccine=0;
        int i;
    //cout<<"In vaccinerToutesVilles_V1"<<endl;
    if (politique == NULL){
        cout<<"politique == NULL"<<endl;
                for(int i=0;i<nbVilles;i++){
                         // nbPersonneVaccine = (y[i][iS]*tauxglobVAC)/100;
                        nbPersonneVaccine=(y[i][iS]*tauxglobVAC)/100; // Nhan : changer tauxglobVAC est entre 0,1
                        y[i][iS] = y[i][iS]-nbPersonneVaccine;
                        //y[i][iR] = y[i][iR]+nbPersonneVaccine;
                        y[i][iV] = y[i][iV]+nbPersonneVaccine; // Nhan : (w0% de S) va dans un compartiment qui s'appelle V (pour vacciné
                                                               // et pas dans R, car il faut que l'on garde dans R que ceux qui
                                                               // viennent de I)
                }
    }
    else
    {
         double tauxglobVACtemp = (*politique)[0];
        //cout<<"politique != NULL"<<endl;
                for(i=0;i<nbVilles;i++){
                 cout << "S=" << y[i][iS] << "E=" << y[i][iE] << "I=" << y[i][iI] << "R=" << y[i][iR] << "V=" << y[i][iV] << endl;

                    //  if (y[i][iS]*(*politique)[0]+y[i][iE]*(*politique)[1]+y[i][iI]*(*politique)[2]+y[i][iR]*(*politique)[3] > 0) {
                 //   if (y[i][iS]*(*politique)[0] + y[i][iI]*(*politique)[1] + y[i][iR]*(*politique)[2] > 100000) {
                if (y[i][iE]*((*politique)[1]) + y[i][iI]*((*politique)[2]) > 500*((*politique)[3])) {
                    //if (y[i][iE]*((*politique)[1]) > y[i][iI]*((*politique)[2])) { //> 200.0*((*politique)[3])
                        // nbPersonneVaccine = (y[i][iS]*tauxglobVAC)/100;
                        nbPersonneVaccine = (y[i][iS]*tauxglobVACtemp); // Nhan : changer tauxglobVAC est entre 0,1
                        //y[i][iS] = std::max((y[i][iS]-nbPersonneVaccine), 0);
                        //cout << "nbPersonneVaccin=" << nbPersonneVaccine << "=" << y[i][iS] << endl;
                        y[i][iS] = (y[i][iS]-nbPersonneVaccine);
                        //y[i][iR] = y[i][iR]+nbPersonneVaccine;
                        y[i][iV] = y[i][iV]+nbPersonneVaccine; // Nhan : (w0% de S) va dans un compartiment qui s'appelle V (pour vacciné
                                                               // et pas dans R, car il faut que l'on garde dans R que ceux qui
                                                               // viennent de I)
                    }
                }
    }

        return nbPersonneVaccine;
}

/******************************************************************************************************/

  /***************************************************************************************************************/
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
  vector<vector<vector<unsigned long> > > simulerSEIR(int nbVilles,double sigma,double gamma,double mu,
                                               double beta0,double beta1,double *phi,
                                               double** arr_rho,
                                               double seed,double unitTIME, double tmax, double typeRNG,
                                               double T, unsigned long **valeursSEIR0,
                                               bool valParVAC, double tauxVaccToutesVilles,
                                               vector<double>*politiqueEVOL){
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
    double** arr_xi;create2D(nbVilles,nbVilles,arr_xi);

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

    vector<double> vectVAC;
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
        do{
            p1 = rand()/(double)RAND_MAX;//probability for time
            p2 = rand()/(double)RAND_MAX;//probability for time
            p3 = rand()/(double)RAND_MAX;//probability to choose event
        }while(p1 == 0.0 || p2 == 0.0 || p3 == 0.0);


        if(valParVAC){
            int nombreSusceptibleVaccinees = vaccinerToutesVilles_V_Nhan(tauxVaccToutesVilles,nbVilles,y,politiqueEVOL);
            nbSuscVac += nombreSusceptibleVaccinees;
             if(nombreSusceptibleVaccinees != 0) vectVAC.push_back(t);
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
  }

  /***** The end of the main program******/
//____________________________________________________________________________________________________

/************LA FIN **************************/
