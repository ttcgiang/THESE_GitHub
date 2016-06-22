#include "seir_stoch.h"
#include "seir_det.h"
#include "equilibrium.h"
#include <R.h>
#include <Rdefines.h>
using namespace std;

//arr_rho is matrix or a number
//epsilon is a matrix or a number

extern "C" {
//lay gia tri tat  ca cac thanh pho....
/*  vector<vector<vector<double> > > simulerSEIR(int nbVilles,vector<double>sigma,vector<double>gamma,vector<double>mu,
                                               vector<double>beta0,vector<double>beta1,vector<double>phi,
                                               double** arr_rho, double** epsilon,
                                               double seed,double unitTIME, double tmax, double typeRNG,
                                               double T, unsigned long S,  unsigned long E,  unsigned long I,  unsigned long R,
                                               vector<vector<double> >valeursSEIR0)

*/

//SEXP getValeurPOPS(SEXP argc0, SEXP argv0, SEXP initvarib) {
SEXP getValeurPOPS(SEXP nbVilles,SEXP sigma,SEXP gamma,SEXP mu,
                   SEXP beta0, SEXP beta1,SEXP phi,
                   SEXP arr_rho, SEXP epsilon,
                   SEXP seed,SEXP unitTIME, SEXP tmax, SEXP typeRNG,
                   SEXP T, SEXP initvarib) {
    try{
	//nbVilles
	if(length(nbVilles)!=1)
		error("invalid number of cities.");
	int cnbVilles = INTEGER_VALUE(nbVilles);

	//sigma
	if (!isVector(sigma)  ||  !isReal(sigma)) {
            error("invalid vector of sigma");
        }
	double *msigma = REAL(sigma);
	vector<double> csigma=vector<double>(cnbVilles,0);
	for(int i=0; i<cnbVilles; i++) csigma[i] = msigma[i];
	//gamma
	if (!isVector(gamma)  ||  !isReal(gamma)) {
            error("invalid vector of gamma");
        }
	double *mgamma = REAL(gamma);
	vector<double> cgamma=vector<double>(cnbVilles,0);
	for(int i=0; i<cnbVilles; i++) cgamma[i] = mgamma[i];
	//mu
	if (!isVector(mu)  ||  !isReal(mu)) {
            error("invalid vector of mu");
        }
	double *mmu = REAL(mu);
	vector<double> cmu=vector<double>(cnbVilles,0);
	for(int i=0; i<cnbVilles; i++) cmu[i] = mmu[i];
	//beta0
	if (!isVector(beta0)  ||  !isReal(beta0)) {
            error("invalid vector of beta0");
        }
	double *mbeta0 = REAL(beta0);
	vector<double> cbeta0=vector<double>(cnbVilles,0);
	for(int i=0; i<cnbVilles; i++) cbeta0[i] = mbeta0[i];
	//beta1
	if (!isVector(beta1)  ||  !isReal(beta1)) {
            error("invalid vector of beta1");
        }
	double *mbeta1 = REAL(beta1);
	vector<double> cbeta1=vector<double>(cnbVilles,0);
	for(int i=0; i<cnbVilles; i++) cbeta1[i] = mbeta1[i];
	//phi
	if (!isVector(phi)  ||  !isReal(phi)) {
            error("invalid vector of phi");
        }
	double *mphi = REAL(phi);
	vector<double> cphi=vector<double>(cnbVilles,0);
	for(int i=0; i<cnbVilles; i++) cphi[i] = mphi[i];
	//arr_rho
	if (!isMatrix(arr_rho)  ||  !isReal(arr_rho)) {
            error("invalid value of coupling rate, it should be a symmetric matrix.");
        }
	double **carr_rho;	
	if(isMatrix(arr_rho)){
		if(INTEGER(getAttrib(arr_rho, R_DimSymbol))[0]>1){
			if ((cnbVilles != INTEGER(getAttrib(arr_rho, R_DimSymbol))[0])||(INTEGER(getAttrib(arr_rho, R_DimSymbol))[0]!=INTEGER(getAttrib(arr_rho, R_DimSymbol))[1])) {
          			error("mismatch between number of cities (%i) and number "
            		      		"of rows in coupling matrix (%i)", cnbVilles,INTEGER(getAttrib(arr_rho, R_DimSymbol))[0]);
	        	}
			else{
				double *marr_rho = REAL(arr_rho);
				for(int i=0; i<cnbVilles; i++){
					carr_rho[i] = new double[cnbVilles];
					for(int j=0; j<cnbVilles; j++){
						carr_rho[i][j] = marr_rho[i*cnbVilles+j];
					}
				}
					
			}
		}
		else
		{
			  double * rho = REAL(arr_rho);
     			  carr_rho = calculerTauxCouplage_00(cnbVilles,rho[0]);
			  /*
			  cout<<"matrix of coupling!"<<endl;
				for(int i=0; i<cnbVilles; i++){
					for(int j=0; j<cnbVilles; j++){
						cout<<"  "<< carr_rho[i][j];
					}
				cout<<endl;
				 }
			*/
		}		
	}
	else
		  error("invalid value of coupling rate, it should be a symmetric matrix.");

	//epsilon
	if (!isMatrix(epsilon)  ||  !isReal(epsilon)) {
            error("invalid value of proportion of infections by contact, it should be a symmetric matrix.");
        }
	double **carr_epsilon;	
	if(isMatrix(epsilon)){
		 if(INTEGER(getAttrib(epsilon, R_DimSymbol))[0]>1){
			if ((cnbVilles != INTEGER(getAttrib(epsilon, R_DimSymbol))[0])||(INTEGER(getAttrib(epsilon, R_DimSymbol))[0]!=INTEGER(getAttrib(epsilon, R_DimSymbol))[1])) {
          			  error("mismatch between number of cities (%i) and number "
            			      "of rows in ' proportion of infections by contact' matrix (%i)", cnbVilles,INTEGER(getAttrib(epsilon, R_DimSymbol))[0]);
	        	}
			else{
				double *marr_epsilon = REAL(epsilon);
				for(int i=0; i<cnbVilles; i++){
					carr_epsilon[i] = new double[cnbVilles];
					for(int j=0; j<cnbVilles; j++){
						carr_epsilon[i][j] = marr_epsilon[i*cnbVilles+j];
					}
				}
			}		
		}
		else{
			 double* epsilon0 = REAL(epsilon);
     			  carr_epsilon =  calculerTauxExterieur_00(cnbVilles,epsilon0[0]);
				/*
				cout<<"matrix of  proportion of infections by contact!"<<endl;
				for(int i=0; i<cnbVilles; i++){
					for(int j=0; j<cnbVilles; j++){
						cout<<"  "<< carr_epsilon[i][j];
					}
				cout<<endl;
				 }
				*/
		}
	}
	else
		error("invalid value of proportion of infections by contact, it should be a symmetric matrix.");
	//seed
	if (!isReal(seed)  ||  length(seed) != 1) {
            error("invalid value of seed,it should be a real number");
        }
	double cseed=NUMERIC_VALUE(seed);
	//unitTIME
	if (!isReal(unitTIME)  ||  length(unitTIME) != 1) {
            error("invalid value of unitTIME,it should be a real number");
        }
	double cunitTIME=NUMERIC_VALUE(unitTIME);
	//tmax
	if (!isReal(tmax)  ||  length(tmax) != 1) {
            error("invalid value of simulation time,it should be a real number");
        }
	double ctmax=NUMERIC_VALUE(tmax);
	//typeRNG
	if (!isReal(typeRNG)  ||  length(typeRNG) != 1) {
            error("invalid value of type generating a random,it should be 0 or 1");
        }
	double ctypeRNG=NUMERIC_VALUE(typeRNG);
	// SEXP T
	if (!isReal(T)  ||  length(T) != 1) {
            error("invalid value of period, it should be a real number");
        }
	double cT=NUMERIC_VALUE(T);
	//vector initvarib
        if (!isVector(initvarib)  ||  !isReal(initvarib)) {
            error("invalid vector of initial values");
        }

/*
    // get params
    int argc = INTEGER_VALUE(argc0);
    char* argv[255];
    for(int i=0; i< argc; i++){
        argv[i] = R_alloc(strlen(CHAR(STRING_ELT(argv0, i))), sizeof(char));
        strcpy( argv[i], CHAR(STRING_ELT(argv0, i)));
    }
*/
    //get the initial values of variable S,E,I,R for each city
    int sumnbvarib = length(initvarib);
    double *m_X = REAL(initvarib);
    vector<vector<double> > villeSEIR0;
    double*ville,*S0,*E0,*I0,*R0;
    int nbVilles = (int)sumnbvarib/4;
    vector<double> tmpVille;
    for(int i=0; i<sumnbvarib; i++){
	for(int j=0; j<nbVilles;j++){
		tmpVille.push_back(j);tmpVille.push_back(m_X[0*nbVilles+j]);
		tmpVille.push_back(m_X[1*nbVilles+j]);tmpVille.push_back(m_X[2*nbVilles+j]);tmpVille.push_back(m_X[3*nbVilles+j]);
		villeSEIR0.push_back(tmpVille);
		tmpVille.clear();
	}
    }
   // vector<double> valParSIM = initialervalParSIM_00(argc,argv);
    vector<vector<vector<double> > > res = simulerSEIR(cnbVilles,csigma,cgamma,cmu,cbeta0,cbeta1,cphi,carr_rho,carr_epsilon,
						cseed,cunitTIME,ctmax,ctypeRNG,cT,0,0,0,0,villeSEIR0);

    //int nbVille = valParSIM[inbVilles];
    int size2D = res[0].size();
    SEXP listVilles, list_names_villes;
    char names_villes[255];
    // objects in out list:
    PROTECT(list_names_villes = allocVector(STRSXP,cnbVilles));
    for(int i=0; i<cnbVilles;i++){
        sprintf(names_villes,"pop%d",i);
        SET_STRING_ELT(list_names_villes,i,mkChar(names_villes));
    }

    PROTECT(listVilles = allocVector(VECSXP, cnbVilles));
    SEXP myT,myS,myE,myI,myR,myN,myP, list, list_names;
    double *valeurT, *valeurS, *valeurE,*valeurI,*valeurR,*valeurN,*valeurP;
    const char *names[7] = {"time(day)", "S", "E", "I", "R", "N","P"};

     for(int idxVille=0; idxVille<cnbVilles;idxVille++){
         // creating an double vector:
         PROTECT(myT = NEW_NUMERIC(size2D));
         valeurT = NUMERIC_POINTER(myT);

         PROTECT(myS = NEW_NUMERIC(size2D));
         valeurS = NUMERIC_POINTER(myS);

         PROTECT(myE = NEW_NUMERIC(size2D));
         valeurE = NUMERIC_POINTER(myE);

         PROTECT(myI = NEW_NUMERIC(size2D));
         valeurI = NUMERIC_POINTER(myI);

         PROTECT(myR = NEW_NUMERIC(size2D));
         valeurR = NUMERIC_POINTER(myR);

         PROTECT(myN = NEW_NUMERIC(size2D));
         valeurN = NUMERIC_POINTER(myN);

	 PROTECT(myP = NEW_NUMERIC(size2D));
         valeurP = NUMERIC_POINTER(myP);

         // Creating a character string vector of the "names" attribute of the objects in out list:
         PROTECT(list_names = allocVector(STRSXP,7));
         for(int i = 0; i < 7; i++)
             SET_STRING_ELT(list_names,i,mkChar(names[i]));

         // Creating a list with 2 vector elements:
          PROTECT(list = allocVector(VECSXP, 7));

         for(int j=0; j<size2D; j++){
            valeurT[j]= res[idxVille][j][0];
            valeurS[j]= res[idxVille][j][1];
            valeurE[j]= res[idxVille][j][2];
            valeurI[j]= res[idxVille][j][3];
            valeurR[j]= res[idxVille][j][4];
            valeurN[j]= res[idxVille][j][5];
	    valeurP[j]= res[idxVille][j][6];
        }

        // attaching myint vector to list:
        SET_VECTOR_ELT(list, 0, myT);
        // attaching mydouble vector to list:
        SET_VECTOR_ELT(list, 1, myS);
        SET_VECTOR_ELT(list, 2, myE);
        SET_VECTOR_ELT(list, 3, myI);
        SET_VECTOR_ELT(list, 4, myR);
        SET_VECTOR_ELT(list, 5, myN);
        SET_VECTOR_ELT(list, 6, myP);
        // and attaching the vector names:
        setAttrib(list, R_NamesSymbol, list_names);
        SET_VECTOR_ELT(listVilles, idxVille, list);
        UNPROTECT(9);
     }

     setAttrib(listVilles,R_NamesSymbol,list_names_villes);
     UNPROTECT(2);
	
     for (int i = 0; i <  cnbVilles; i++)
    {
        delete []carr_rho[i];
        delete []carr_epsilon[i];
    }
    delete []carr_rho; delete[]carr_epsilon;

     return listVilles;
    }catch (exception &e) {
            error(e.what());
            return(0);
        }
    
  }
 }
}


