#include "seirStoGENE.h" //library of stochastic functions
#include "initialFUNC.h"
#include "seirStoYANN.h"


//Base library of R
#include <R.h>
#include <Rdefines.h>

using namespace std;

//arr_rho is matrix or a number
//epsilon is a matrix or a number

extern "C" {

// Gathering all values of parameters and variables from R to C++
SEXP getValeurPOPSGENE(SEXP nbVilles, SEXP sigma, SEXP gamma, SEXP mu,
                   SEXP beta0, SEXP beta1, SEXP phi,
                   SEXP rho0, //SEXP epsilon,
                   SEXP seed, SEXP unitTIME, SEXP tmax, SEXP typeRNG,
                   SEXP T, SEXP initvarib) {
    try{
        //Number of cities in the metapopulation
        if(!isReal(nbVilles) || length(nbVilles)!=1)
            error("invalid number of cities.");
        int cnbVilles = INTEGER_VALUE(nbVilles);

        //sigma
        if (!isReal(sigma) || length(sigma)!=1) {
            error("invalid parameter sigma");
        }
        double csigma = NUMERIC_VALUE(sigma);

        //gamma
        if (!isReal(gamma) || length(gamma)!=1) {
            error("invalid parameter gamma");
        }
        double cgamma = NUMERIC_VALUE(gamma);

        //mu
        if (!isReal(mu) || length(mu) !=1) {
            error("invalid parameter mu");
        }
        double cmu = NUMERIC_VALUE(mu);

        //beta0
        if (!isReal(beta0) || length(beta0)!=1) {
            error("invalid parameter beta0");
        }
        double cbeta0 = NUMERIC_VALUE(beta0);

        //beta1
        if (!isReal(beta1) || length(beta1)!=1) {
            error("invalid vector parameter beta1");
        }
        double cbeta1 = NUMERIC_VALUE(beta1);

        //phi vector phi
        if (!isVector(phi)  ||  !isReal(phi)) {
                error("invalid vector of phi");
            }
        double *mphi = REAL(phi);
        double* cphi; create1D(cnbVilles,0.0,cphi);
        for(int i=0; i<cnbVilles; i++) cphi[i] = mphi[i];


        //taux rho0, prob pour que un individu x visite une ville j
        if (!isReal(rho0) || length(rho0)!=1) {
            error("invalid vector parameter rho0");
        }
        double crho0 = NUMERIC_VALUE(rho0);

        //number 'seed'
        if (!isReal(seed)  ||  length(seed) != 1) {
            error("invalid value of seed,it should be a real number");
        }
        double cseed=NUMERIC_VALUE(seed);

        // time unit to gather the values of variables SEIR
        if (!isReal(unitTIME)  ||  length(unitTIME) != 1) {
            error("invalid value of unitTIME,it should be a real number");
        }
        double cunitTIME=NUMERIC_VALUE(unitTIME);

        //simulation time tmax
        if (!isReal(tmax)  ||  length(tmax) != 1) {
            error("invalid value of simulation time,it should be a real number");
        }
        double ctmax=NUMERIC_VALUE(tmax);

        //type of the random number generator 1: Yann, 0: C++
        if (!isReal(typeRNG)  ||  length(typeRNG) != 1) {
            error("invalid value of type generating a random,it should be 0 or 1");
        }
        double ctypeRNG=NUMERIC_VALUE(typeRNG);

        // Disease period
        if (!isReal(T)  ||  length(T) != 1) {
            error("invalid value of period, it should be a real number");

        }
        double cT=NUMERIC_VALUE(T);

        //vector of initial values of variables SEIR for all cities
        // Tranform the initial values of variables SEIR from R to C++
        //get the initial values of variable S,E,I,R for each city from R
        if (!isVector(initvarib)  ||  !isReal(initvarib)) {
            error("invalid vector of initial values");
        }
        int sumnbvarib = length(initvarib);
        double *m_X = REAL(initvarib);
        int nvar = 6;
        int rnbVilles = (int)sumnbvarib/4;
        unsigned long **villeSEIR0;  villeSEIR0 = new unsigned long*[rnbVilles];
        for (int i = 0; i < rnbVilles; i++)
        {
            villeSEIR0[i] = new unsigned long[nvar];
            villeSEIR0[i][0]=0;
            villeSEIR0[i][iS] = m_X[0*rnbVilles+i];
            villeSEIR0[i][iE] = m_X[1*rnbVilles+i];
            villeSEIR0[i][iI] = m_X[2*rnbVilles+i];
            villeSEIR0[i][iR] = m_X[3*rnbVilles+i];
            villeSEIR0[i][iN] = villeSEIR0[i][iS]+ villeSEIR0[i][iE]+ villeSEIR0[i][iI]+ villeSEIR0[i][iR];
        }
        //printing
        for(int i=0; i<rnbVilles; i++){
            for (int j=0; j<nvar; j++)
                cout<<"    " << villeSEIR0[i][j];
            cout<<endl;
        }

        //SIMULATION

        // initialiser les valeurs des phi /arr_rho
       // double *cphi = calculerPhases(cnbVilles,cphiMAX,cphiMIN);

        //rho: pro de visite
         double** carr_rho= calculerProbVISITER(cnbVilles, crho0);

        // Doing simulation by using the function 'simulerSEIR'
        vector<vector<vector<unsigned long> > > res = simulerSEIRGENE(cnbVilles,csigma,cgamma,cmu,cbeta0,cbeta1,cphi,carr_rho,
                                                           cseed,cunitTIME,ctmax,ctypeRNG,cT,villeSEIR0);


        //RESULTAT
        // Gathering the result, then send this from C++ to R
        int size2D = res[0].size();
        SEXP listVilles, list_names_villes;
        char names_villes[255];
        // objects in our list:
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


      // freeing the memory of pointer
        for (int i = 0; i <  cnbVilles; i++)
        {
            delete []carr_rho[i];            
        }
        delete []carr_rho;
        delete []cphi;
        return listVilles;

    }catch (exception &e) {
        error(e.what());
        return(0);
    }
}
}

/***************The end*********************/

