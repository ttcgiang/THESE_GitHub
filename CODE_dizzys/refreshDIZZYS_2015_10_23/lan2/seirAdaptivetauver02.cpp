/* $Id: adaptivetau.cpp 298 2014-11-19 21:13:54Z pjohnson $
    --------------------------------------------------------------------------
    C++ implementation of the "adaptive tau-leaping" algorithm described by
    Cao Y, Gillespie DT, Petzold LR. The Journal of Chemical Physics (2007).
    Author: Philip Johnson <plfjohnson@emory.edu>


    Copyright (C) 2010 Philip Johnson

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    If building library outside of R package (i.e. for debugging):
        R CMD SHLIB adaptivetau.cpp
    --------------------------------------------------------------------------
*/

#include <iostream>
#include <vector>
#include <limits>
#include <sstream>
#include <stdexcept>
#include "seirStoYANN.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>
#include "Rwrappers.h"


using namespace std;

enum EStepType {
    eExact = 0,
    eExplicit,
    eImplicit
};

const bool debug = false;

class CEarlyExit : public runtime_error {
public:
    CEarlyExit(const string &w) : runtime_error(w) {}
};

//use below rather than R's "error" directly (which will not free memory, etc.)
#ifdef throwError
#undef throwError
#endif
#define throwError(e) { ostringstream s; s << e; throw runtime_error(s.str()); }
#ifdef throwEarlyExit
#undef throwEarlyExit
#endif
#define throwEarlyExit(e) { ostringstream s; s << e << "; results returned only up until this point"; throw CEarlyExit(s.str());}
// CODE of TRAN THI CAM GIANG that are added in the c++ codes of the package "adaptive tau"
extern "C" { 
/*****/
/*
    Goal: calculating the propensity functions for all subpopulation.
        Used for SEIR model
    */
double* seirratefunc(double* x,int nbVilles,double** arr_probVISITER,double** arr_probVENIRk,double** arr_probINFECTER,
        double** arr_nbCONTACT0,double** arr_nbCONTACT1,
        double periDISE, double m_T,double *phi,double** arr_nbCONTACTFORCING,double mu,double sigma, double gamma){
    // 4 variables (S, I, R, N) and 6 events
    int nbVarSEIRN=5, nevent=8;
    // index for the positions of variables
    int S=0, E=1,I=2, R=3, N=4;
    //Initializing
    double *lamda = new double[nbVilles];

    //Converting 1D vector to 2D vector **y
    unsigned long **y;
    y = new unsigned long  *[nbVilles];
    for (int i = 0; i<nbVilles; i++)
    {
        y[i] = new unsigned long [nbVarSEIRN];
        y[i][0] = (unsigned long)x[0*nbVilles+i];//iS
        y[i][1] = (unsigned long)x[1*nbVilles+i];//iE
        y[i][2] = (unsigned long)x[2*nbVilles+i];//iI
        y[i][3] = (unsigned long)x[3*nbVilles+i];//iR
        y[i][4] = y[i][S]+y[i][E]+y[i][I]+y[i][R];//iN

    }

      // recalculating the matrix "COME from the city k
        calculerProbVENIRk(nbVilles,y,arr_probVISITER,arr_probVENIRk);

       //recalculing the nbCONATCT forcing
        calculateNbCONTACTForcing(nbVilles,arr_nbCONTACT0,arr_nbCONTACT1,periDISE,m_T,phi,arr_nbCONTACTFORCING);
      //Calculating Value of \lamda at time t
        calLamdaYANNGIANG(nbVilles,arr_probVISITER,y,arr_nbCONTACTFORCING,arr_probINFECTER,arr_probVENIRk,lamda);

        //Calculating the values of all propensity functions
        double **f;
        f = new double*[nbVilles];
        for (int i = 0; i <  nbVilles; i++)
        {
            f[i] = new double[nevent];
            for(int j = 0; j< nevent ; j++) f[i][j] = 0.0;
        }

        for(int i=0; i< nbVilles;i++){
            //S born
            f[i][0] = mu*y[i][N];
            //S die
            f[i][1] = mu*y[i][S];
             //E die
            f[i][2] = mu*y[i][E];
            //I die
            f[i][3] = mu*y[i][I];
            //R die
            f[i][4] = mu*y[i][R];
            //S-> E
            f[i][5] = lamda[i]*y[i][S];
            //E-> I
            f[i][6] = sigma*y[i][S];
            //I->R
            f[i][7] = gamma*y[i][I];
        }

        // vector of result
        vector<double> vecRes;
        for(int j = 0; j< nevent ; j++)
        {
               for (int i = 0; i <  nbVilles; i++)
               vecRes.push_back(f[i][j]);
        }
        int nbvecRes=vecRes.size();
        double *resRates = new double[nbvecRes];
        for(int i=0; i<vecRes.size();i++){
            resRates[i] = vecRes[i];
        }

        // Freeing the memory of all dynnamic pointers
        delete []lamda;
        for (int i = 0; i<nbVilles; i++)
        {
            delete []y[i];
            delete []f[i];
        }
        delete []y; delete[]f;
    return(resRates);
  }
/*****/
/*
   Goal: calculating the propensity functions for all subpopulation.
       Used for SIR model
   */
double* sirratefunc(double* x,int nbVilles,double** arr_probVISITER,double** arr_probVENIRk,double** arr_probINFECTER,
       double** arr_nbCONTACT0,double** arr_nbCONTACT1,
       double periDISE, double m_T,double *phi,double** arr_nbCONTACTFORCING, double mu,double gamma){
   // 4 variables (S, I, R, N) and 6 events
   int nbVarSIRN=4, nevent=6;
   // index for the positions of variables
   int S=0, I=1, R=2, N=3;
   //Initializing
   double *lamda = new double[nbVilles];

   //Converting 1D vector to 2D vector **y
   unsigned long **y;
   y = new unsigned long  *[nbVilles];
   for (int i = 0; i<nbVilles; i++)
   {
       y[i] = new unsigned long [nbVarSIRN];
       y[i][0] = (unsigned long)x[0*nbVilles+i];//iS
       y[i][1] = (unsigned long)x[1*nbVilles+i];//iI
       y[i][2] = (unsigned long)x[2*nbVilles+i];//iR
       y[i][3] = y[i][S]+y[i][I]+y[i][R];//iN

   }

     // recalculating the matrix "COME from the city k
       calculerProbVENIRk(nbVilles,y,arr_probVISITER,arr_probVENIRk);

      //recalculing the nbCONATCT forcing
       calculateNbCONTACTForcing(nbVilles,arr_nbCONTACT0,arr_nbCONTACT1,periDISE,m_T,phi,arr_nbCONTACTFORCING);
     //Calculating Value of \lamda at time t
       calLamdaYANNGIANG(nbVilles,arr_probVISITER,y,arr_nbCONTACTFORCING,arr_probINFECTER,arr_probVENIRk,lamda);

       //Calculating the values of all propensity functions
       double **f;
       f = new double*[nbVilles];
       for (int i = 0; i <  nbVilles; i++)
       {
               f[i] = new double[nevent];
     for(int j = 0; j< nevent ; j++) f[i][j] = 0.0;
       }

       for(int i=0; i< nbVilles;i++){
           //S born
           f[i][0] = mu*y[i][N];
           //S die
           f[i][1] = mu*y[i][S];
           //I die
           f[i][2] = mu*y[i][I];
           //R die
           f[i][3] = mu*y[i][R];
           //S-> I
           f[i][4] = lamda[i]*y[i][S];
           //I->R
           f[i][5] = gamma*y[i][I];
       }

       // vector of result
       vector<double> vecRes;
       for(int j = 0; j< nevent ; j++)
       {
              for (int i = 0; i <  nbVilles; i++)
              vecRes.push_back(f[i][j]);
       }
       int nbvecRes=vecRes.size();
       double *resRates = new double[nbvecRes];
       for(int i=0; i<vecRes.size();i++){
           resRates[i] = vecRes[i];
       }

       // Freeing the memory of all dynnamic pointers
       delete []lamda;
       for (int i = 0; i<nbVilles; i++)
       {
           delete []y[i];
           delete []f[i];
       }
     delete []y; delete[]f;
     return(resRates);
 }

/*****/
/*****/
}
/**************/
class CStochasticEqns {
public:
/*
    CStochasticEqns(SEXP initVal, SEXP nu,
                    SEXP rateFunc, SEXP rateJacobianFunc,
                    SEXP params, double* changeBound, SEXP maxTauFunc,
                    SEXP detTrans, SEXP haltTrans) {
*/
CStochasticEqns(SEXP initVal, int *nu, unsigned int numTrans,
                    SEXP rateFunc, SEXP rateJacobianFunc,
                    // for model SEIR
                    int nbVilles,double sigma,double gamma,double mu,
                    // for metapopulation
                    double** arr_probVISITER,double** arr_nbCONTACT0,double** arr_nbCONTACT1, double** arr_probINFECTER, double** arr_probVENIRk,
	     double *phiPHASE,double tmax,double periDISE,
                    //
                    double* changeBound, SEXP maxTauFunc, SEXP detTrans) {
        // copy initial values into new vector (keeping in SEXP vector
        // allows easy calling of R function to calculate rates)
        m_NumStates = length(initVal);
        SEXP x;
        PROTECT(x = allocVector(REALSXP, m_NumStates));
        copyVector(x, coerceVector(initVal,REALSXP));
        if (isNull(getAttrib(initVal, R_NamesSymbol))) {
            m_VarNames = NULL;
        } else {
            SEXP namesO = getAttrib(initVal, R_NamesSymbol);

            PROTECT(m_VarNames = allocVector(STRSXP, length(namesO)));

            copyVector(m_VarNames, namesO);
            setAttrib(x, R_NamesSymbol, m_VarNames);
        }
        m_X = REAL(x);

        // copy full-size Nu matrix into sparse matrix
        m_Nu.resize(numTrans);
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            for (unsigned int j = 0;  j < numTrans;  ++j) {
                if (nu[j*m_NumStates + i] != 0) {
                    m_Nu[j].push_back(SChange(i, nu[j*m_NumStates + i]));
                }
            }
        }
        m_TransCats.resize(numTrans, eNoncritical);

        // potentially flag some transitions as "deterministic"
        x_SetCat(detTrans, eDeterministic);
        if (m_TransByCat[eDeterministic].size() == m_TransCats.size()) {
            throwError("At least one transition must be stochastic (all "
                       "transitions are currently flagged as "
                       "deterministic).");
        }

        // potentially flag some transitions as "halting" (which are
        // always critical)
        x_SetCat(haltTrans, eHalting);
        m_TransByCat[eCritical] = m_TransByCat[eHalting];

        // needed for ITL
        x_IdentifyBalancedPairs();
        x_IdentifyRealValuedVariables();

        // prepare R function for evaluation by setting up arguments
        // (current X values, parameters, current time)
        SEXP s_time;
        PROTECT(s_time = allocVector(REALSXP, 1));
        m_T = REAL(s_time);
        PROTECT(m_RateFunc = lang4(rateFunc, x, params, s_time));
        if (!rateJacobianFunc  ||  isNull(rateJacobianFunc)) {
            m_RateJacobianFunc = NULL;
        } else {
            PROTECT(m_RateJacobianFunc = lang4(rateJacobianFunc, x,
                                               params, s_time));
        }
        m_Rates = NULL;

        //default parameters to adaptive tau leaping algorithm
        m_Epsilon = 0.05;
        m_Ncritical = 10;
        m_Nstiff = 100;
        m_ExactThreshold = 10;
        m_Delta = 0.05;
        m_NumExactSteps[eExact] = 100;
        m_NumExactSteps[eExplicit] = 100;
        m_NumExactSteps[eImplicit] = 10;
        m_ITLConvergenceTol = 0.01;
        m_MaxTau = numeric_limits<double>::infinity();
        m_MaxSteps = 0; // special case 0 == no limit

        //useful additional parameters
        m_ExtraChecks = true;
        m_VerboseTracing = 0;
        m_RateChangeBound = changeBound;
        if (!maxTauFunc  ||  isNull(maxTauFunc)) {
            m_MaxTauFunc = NULL;
        } else {
            PROTECT(m_MaxTauFunc = lang4(maxTauFunc, x, params, s_time));
        }

        //check initial conditions to make sure legit
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            if (m_X[i] < 0) {
                throwError("initial value for variable " << i+1 <<
                           " must be positive (currently " << m_X[i] << ")");
            }
            if (!m_RealValuedVariables[i]  &&  (m_X[i] - trunc(m_X[i]) > 1e-5)){
                if (m_VarNames != NULL) {
                    throwError("initial value for variable " << i+1 <<
                               " ('"<<CHAR(STRING_PTR(m_VarNames)[i])<<"') " <<
                               "must be an integer (currently " << m_X[i] << ")");
                } else {
                    throwError("initial value for variable " << i+1 <<
                               " must be an integer (currently " << m_X[i] << ")");
                }
            }
        }
        
        *m_T = 0;
        m_LastTransition = -1;
        m_PrevStepType = eExact;
        GetRNGstate();
    }
    ~CStochasticEqns(void) {
        int cnt = 3;
        if (m_RateJacobianFunc != NULL) {
            ++cnt;
        }
        if (m_Rates != NULL) {
            ++cnt;
        }
        if (m_MaxTauFunc != NULL) {
            ++cnt;
        }
        if (m_VarNames != NULL) {
            ++cnt;
        }
        UNPROTECT(cnt);
    }
    void SetTLParams(SEXP list) {
        SEXP names = getAttrib(list, R_NamesSymbol);

        for (int i = 0;  i < length(names);  ++i) {
            if (strcmp("epsilon", CHAR(STRING_PTR(names)[i])) == 0) {
                if (!isReal(VECTOR_ELT(list, i))  ||
                    length(VECTOR_ELT(list, i)) != 1) {
                    throwError("invalid value for parameter '" <<
                               CHAR(STRING_PTR(names)[i]) << "'");
                }
                m_Epsilon = REAL(VECTOR_ELT(list, i))[0];
            } else if (strcmp("delta", CHAR(STRING_PTR(names)[i])) == 0) {
                if (!isReal(VECTOR_ELT(list, i))  ||
                    length(VECTOR_ELT(list, i)) != 1) {
                    throwError("invalid value for parameter '" <<
                               CHAR(STRING_PTR(names)[i]) << "'");
                }
                m_Delta = REAL(VECTOR_ELT(list, i))[0];
            } else if (strcmp("maxtau", CHAR(STRING_PTR(names)[i])) == 0) {
                if (!isReal(VECTOR_ELT(list, i))  ||
                    length(VECTOR_ELT(list, i)) != 1) {
                    throwError("invalid value for parameter '" <<
                               CHAR(STRING_PTR(names)[i]) << "'");
                }
                m_MaxTau = REAL(VECTOR_ELT(list, i))[0];
            } else if (strcmp("extraChecks",
                              CHAR(STRING_PTR(names)[i])) == 0) {
                if (!isLogical(VECTOR_ELT(list, i))  ||
                    length(VECTOR_ELT(list, i)) != 1) {
                    throwError("invalid value for parameter '" <<
                               CHAR(STRING_PTR(names)[i]) << "'");
                }
                m_ExtraChecks = LOGICAL(VECTOR_ELT(list, i))[0];
            } else if (strcmp("verbose",
                              CHAR(STRING_PTR(names)[i])) == 0) {
                if (!isInteger(VECTOR_ELT(list, i))  ||
                    length(VECTOR_ELT(list, i)) != 1) {
                    throwError("invalid value for parameter '" <<
                               CHAR(STRING_PTR(names)[i]) << "'");
                }
                m_VerboseTracing = INTEGER(VECTOR_ELT(list, i))[0];
            } else {
                warning("ignoring unknown parameter '%s'",
                        CHAR(STRING_PTR(names)[i]));
            }
        }
    }

    // Functions below are a hack suggested by Simon Urbanek (although
    // he "would not recommend for general use") to check if the user has
    // asked to interrupt execution.  The problem with calling
    // R_CheckUserInterrupt directly is it longjmps and we don't have
    // a chance to free memory from the heap.
    // http://tolstoy.newcastle.edu.au/R/e13/devel/11/04/1049.html
    static void chkIntFn(void*) { R_CheckUserInterrupt(); }
    bool checkUserInterrupt(void) {
        return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
    }

    void EvaluateATLUntil(double tF) {
        unsigned int c = 0;
        //add initial conditions to time series
        m_TimeSeries.push_back(STimePoint(0, m_X, m_NumStates));
        //main loop
        while (*m_T < tF  &&  (m_MaxSteps == 0 || c < m_MaxSteps)  &&
               (m_LastTransition < 0  ||
                m_TransCats[m_LastTransition] != eHalting)) {
            x_UpdateRates();
            x_SingleStepATL(tF);
            if (++c % 10 == 0  &&  checkUserInterrupt()) {
                throwEarlyExit("simulation interrupted by user at time " << *m_T
                               << " after " << c << " time steps.");
            }
        }
        //save RNG state back to R (could also do in destructor, but
        //no harm in extra calls to PutRNGstate and avoids potential
        //problems with PROTECTing return value from GetTimeSeriesSEXP
        PutRNGstate();
    }
    void EvaluateExactUntil(double tF) {
        unsigned int c = 0;
        //add initial conditions to time series
        m_TimeSeries.push_back(STimePoint(0, m_X, m_NumStates));
        m_LastTransition = -1;
        //main loop
        while (*m_T < tF  &&  (m_MaxSteps == 0 || c < m_MaxSteps)  &&
               (m_LastTransition < 0  ||
                m_TransCats[m_LastTransition] != eHalting)) {
            x_UpdateRates();
            x_SingleStepExact(tF);
            if (++c % 10 == 0  &&  checkUserInterrupt()) {
                throwEarlyExit("simulation interrupted by user at time " << *m_T
                               << " after " << c << " time steps.");
            }
        }
        //save RNG state back to R (could also do in destructor, but
        //no harm in extra calls to PutRNGstate and avoids potential
        //problems with PROTECTing return value from GetTimeSeriesSEXP
        PutRNGstate();
    }

    SEXP GetResult(void) const {
        if (m_TransByCat[eHalting].size() == 0) {
            return GetTimeSeriesSEXP();
        } else {
            CRList res(2, true);
            res.SetSEXP(0, GetTimeSeriesSEXP(), "dynamics");
            CRVector<int> lastTrans(1, false);
            lastTrans[0] = (m_LastTransition < 0  ||
                            m_TransCats[m_LastTransition] != eHalting) ?
                NA_INTEGER : m_LastTransition+1;
            res.SetSEXP(1, lastTrans, "haltingTransition");
            return res;
        }
    }

    SEXP GetTimeSeriesSEXP(void) const {
        SEXP res;
        PROTECT(res = allocMatrix(REALSXP, m_TimeSeries.size(), m_NumStates+1));
        double *rvals = REAL(res);
        for (unsigned int t = 0;  t < m_TimeSeries.size();  ++t) {
            rvals[t] = m_TimeSeries[t].m_T;
            for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                rvals[(i+1) * m_TimeSeries.size() + t] = m_TimeSeries[t].m_X[i];
            }
        }

        SEXP dimnames, colnames;
        PROTECT(dimnames = allocVector(VECSXP, 2));
        PROTECT(colnames = allocVector(VECSXP, m_NumStates+1));
        SET_VECTOR_ELT(dimnames, 1, colnames);
        SET_VECTOR_ELT(colnames, 0, mkChar("time"));
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            if (m_VarNames  &&  (unsigned int)length(m_VarNames) > i) {
                SET_VECTOR_ELT(colnames, i+1,
                               STRING_PTR(m_VarNames)[i]);
            } else {
                char name[10];
                snprintf(name, 10, "x%i", i+1);
                SET_VECTOR_ELT(colnames, i+1, mkChar(name));
            }
        }
        setAttrib(res, R_DimNamesSymbol, dimnames);

        UNPROTECT(3);
        return res;
    }

protected:
    enum ETransCat {
        eNormal = 0,
        eCritical,
        eDeterministic,
        eHalting
    };
    typedef vector<ETransCat> TTransCats;
    typedef vector<int> TTransList;
    typedef vector<pair<unsigned int, unsigned int> > TBalancedPairs;
    typedef vector<bool> TBools;
    typedef double* TStates;
    typedef double* TRates;
    struct SChange {
        short int m_State;
        short int m_Mag;
    };
    typedef vector< vector<SChange> > TTransitions;
    struct STimePoint {
        STimePoint(double t, double *x, int n) {
            m_T = t;
            m_X = new double[n];
            memcpy(m_X, x, n*sizeof(double));
        }
        double m_T;
        double *m_X;
    };
    class CTimeSeries : public vector<STimePoint> {
    public:
        ~CTimeSeries(void) {
            for (iterator i = begin();  i != end();  ++i) {
                delete[] i->m_X; i->m_X = NULL;
            }
        }
    };

protected:
    void x_IdentifyBalancedPairs(void);
    void x_IdentifyRealValuedVariables(void);
    void x_SetCat(SEXP trans, ETransCat cat);

    void x_AdvanceDeterministic(double deltaT, bool clamp = false);
    void x_SingleStepExact(double tf);
    void x_SingleStepETL(double tau);
    void x_SingleStepITL(double tau);
    void x_SingleStepATL(double tf);

    void x_UpdateRates(void) {
        if (m_ExtraChecks) {
            for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                if (m_X[i] < 0) {
                    throwError("negative variable: " << i+1 << " is " <<
                               m_X[i] << " (check rate function "
                               "and/or transition matrix)");
                } else if (ISNAN(m_X[i])) {
                    throwError("NaN variable: " << i+1 << " is " <<
                               m_X[i] << " (check rate function "
                               "and/or transition matrix)");
                }
            }
        }

        //sigh.  Cover ourselves in case rateFunc code messes with
        //RNGs (which really should NOT be the case -- we rely on the
        //rate functions being deterministic!), but ran into this
        //problem when using a Rcpp rate function (which adds
        //arbitrary calls to Get/Put RNG).
        PutRNGstate(); 

        // make sure our rates are protected!
        if (m_Rates != NULL) { 
            UNPROTECT(1);
        }
        SEXP res = eval(m_RateFunc, R_EmptyEnv);
        PROTECT(res);
        m_Rates = REAL(res);

        if ((unsigned int) length(res) != m_Nu.size()) {
            throwError("invalid rate function -- returned number of rates ("
                       << length(res) << ") is not the same as specified by "
                       "the transition matrix (" << m_Nu.size() << ")!");
        }
        if (m_ExtraChecks) {
            for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
                if (ISNAN(m_Rates[j])) {
                    throwError("invalid rate function -- rate for transition "
                               << j+1 << " is not a number (NA/NaN)! (check "
                               "for divison by zero or similar)");
                }
                if (m_Rates[j] < 0) {
                    throwError("invalid rate function -- rate for transition "
                               << j+1 << " is negative!");
                }
            }
        }
    }
    double* x_CalcJacobian(void) {
        SEXP res = eval(m_RateJacobianFunc, R_EmptyEnv);
        if (!isMatrix(res)) {
            throwError("invalid Jacobian function -- should return a " <<
                       m_NumStates << " by " << m_Nu.size() << " matrix");
        }
        unsigned int nrow = INTEGER(getAttrib(res, R_DimSymbol))[0];
        unsigned int ncol = INTEGER(getAttrib(res, R_DimSymbol))[1];
        if (nrow != m_NumStates  ||  ncol != m_Nu.size()) {
            throwError ("invalid Jacobian function -- returned a " << nrow
                        << " by " << ncol << " matrix instead of the expected "
                        << m_NumStates << " by " << m_Nu.size() <<
                        " (variables by transitions)");
        }
        return REAL(res);
    }
    double x_CalcUserMaxTau(void) {
        if (!m_MaxTauFunc) { throwError("logic error at line " << __LINE__) }
        SEXP res = eval(m_MaxTauFunc, R_EmptyEnv);
        if (length(res) != 1  || !isReal(res)) {
            throwError("invalid return value from maxTau function (should be "
                       "a single real number)");
        }
        return REAL(res)[0];
    }

    unsigned int x_PickCritical(double prCrit) const;

    double x_TauEx(void) const {
        double tau = numeric_limits<double>::infinity();
        vector <double> mu(m_NumStates, 0);
        vector <double> sigma(m_NumStates, 0);

        for (TTransList::const_iterator j = m_TransByCat[eNormal].begin();
             j != m_TransByCat[eNormal].end();  ++j) {
            for (unsigned int i = 0;  i < m_Nu[*j].size();  ++i) {
                const SChange &c = m_Nu[*j][i];
                mu[c.m_State] += c.m_Mag * m_Rates[*j];
                sigma[c.m_State] += c.m_Mag * c.m_Mag * m_Rates[*j];
            }
        }
        //cerr << "-=| mu:";
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            //cerr << "\t" << mu[i];
            double val = max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                             1.)/fabs(mu[i]);
            //cerr << "/" << val;
            if (val < tau) {
                tau = val;
                if (tau < 0) {
                    throwError("tried to select tau < 0; most likely means "
                               "your rate function generated a negative rate");
                }
            }
            val = pow(max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                          1.),2) / sigma[i];
            if (val < tau) {
                tau = val;
                if (tau < 0) {
                    throwError("tried to select tau < 0; most likely means "
                               "your rate function generated a negative rate");
                }
            }
        }
        //cerr << endl;

        return tau;
    }

    double x_TauIm(void) const {
        if (!m_RateJacobianFunc) {
            return 0;
        }
        vector<bool> equil(m_TransCats.size(), false);
        for (TBalancedPairs::const_iterator i = m_BalancedPairs.begin();
             i != m_BalancedPairs.end();  ++i) {
            if (fabs(m_Rates[i->first] - m_Rates[i->second]) <=
                m_Delta * min(m_Rates[i->first], m_Rates[i->second])) {
                equil[i->first] = true;
                equil[i->second] = true;
            }
        }

        vector<double> mu(m_NumStates, 0);
        vector<double> sigma(m_NumStates, 0);
        for (TTransList::const_iterator j = m_TransByCat[eNormal].begin();
             j != m_TransByCat[eNormal].end();  ++j) {
            if (!equil[*j]) {
                for (unsigned int i = 0;  i < m_Nu[*j].size();  ++i) {
                    const SChange &c = m_Nu[*j][i];
                    mu[c.m_State] += c.m_Mag * m_Rates[*j];
                    sigma[c.m_State] += c.m_Mag * c.m_Mag * m_Rates[*j];
                }
            }
        }

        double tau = numeric_limits<double>::infinity();
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            double val = max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                             1.)/fabs(mu[i]);
            if (val < tau) {
                tau = val;
            }
            val = pow(max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                          1.),2) / sigma[i];
            if (val < tau) {
                tau = val;
            }
        }
        
        return tau;
    }

private:
    bool m_ExtraChecks; //turns on extra checks on rates returned by
                        //user-supplied rate function. Slower, but if
                        //the rate function does have a bug, this will
                        //give a more meaningful error message.
    int m_VerboseTracing; //trace algorithm verbosely

    // parameters to tau leaping algorithm
    unsigned int m_Ncritical;
    double m_Nstiff;
    double m_Epsilon;
    double m_ExactThreshold;
    double m_Delta;
    unsigned int m_NumExactSteps[3];
    double m_ITLConvergenceTol;
    double m_MaxTau;
    unsigned int m_MaxSteps;

    // time-dependent variables
    double *m_T;    // *current* time
    TStates m_X;    // *current* state variables
    TRates m_Rates; // *current* rates (must be updated if m_X changes!)
    EStepType m_PrevStepType; // type of last step
    int m_LastTransition; // id of transition taken if critical/exact; -1 o.w.

    // constant variables
    unsigned int m_NumStates; //total number of states
    SEXP m_VarNames;          //variable names (if any)
    TTransitions m_Nu;        //state changes caused by transitions
    TTransCats  m_TransCats;  //i.e. normal, deterministic, halting
    TTransList  m_TransByCat[4];//i.e. critical, normal, deterministic, halting
    TBalancedPairs m_BalancedPairs;
    TBools m_RealValuedVariables;
    SEXP m_RateFunc;         //R function to calculate rates as f(m_X)
    SEXP m_RateJacobianFunc; //R function to calculate Jacobian of rates as f(m_X) [optional!]
    double *m_RateChangeBound; //see Cao (2006) for details
    SEXP m_MaxTauFunc; //R function to calculate maximum leap given curr. state

    CTimeSeries m_TimeSeries;
};


/*---------------------------------------------------------------------------*/
// PRE : m_Nu initialized
// POST: all balanced pairs of transitions identified & saved
void CStochasticEqns::x_IdentifyBalancedPairs(void) {
    for (unsigned int j1 = 0;  j1 < m_Nu.size();  ++j1) {
        for (unsigned int j2 = j1 + 1;  j2 < m_Nu.size();  ++j2) {
            if (m_Nu[j1].size() != m_Nu[j2].size()) {
                continue;
            }
            unsigned int i;
            for (i = 0;  i < m_Nu[j1].size()  &&
                     m_Nu[j1][i].m_State == m_Nu[j2][i].m_State  &&
                     m_Nu[j1][i].m_Mag == -m_Nu[j2][i].m_Mag;  ++i);
            if (i == m_Nu[j1].size()) {
                m_BalancedPairs.push_back(TBalancedPairs::value_type(j1, j2));
                if (debug) {
                    cerr << "balanced pair " << j1 << " and " << j2 << endl;
                }
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
// PRE : boolean vector flagging transitions as category "cat"
// POST: appropriate transCats set
void CStochasticEqns::x_SetCat(SEXP trans, ETransCat cat) {
    if (!trans  || isNull(trans)) { return; } //NULL may be passed as a flag
    if (isLogical(trans)) {
        CRVector<bool> logic(trans);
        if (logic.size() > m_TransCats.size()) {
            throwError("length of logical vector specifying deterministic or "
                       "halting transitions is greater than the total number "
                       "of transitions!");
        }
        for (unsigned int i = 0;  i < logic.size();  ++i) {
            if (logic[i]) {
                m_TransCats[i] = cat;
                m_TransByCat[cat].push_back(i);
            }
        }
    } else {
        CRVector<int> w(coerceVector(trans, INTSXP), true);
        for (unsigned int i = 0;  i < w.size();  ++i) {
            if (w[i] > m_TransCats.size()) {
                throwError("one of your list(s) of transitions references a "
                           "transition that doesn't exist (" << w[i] << ") "
                           "when last transition is " << m_TransCats.size() <<
                           ")")
            }
            m_TransCats[w[i]-1] = cat;
            m_TransByCat[cat].push_back(w[i]-1);
        }
    }
}

/*---------------------------------------------------------------------------*/
// PRE : m_Nu initialized, deterministic transition set (if any)
// POST: all variables identified will take real values
// (i.e. either non-integer nu or modified by a deterministic transition)
void CStochasticEqns::x_IdentifyRealValuedVariables(void) {
    m_RealValuedVariables.clear();
    m_RealValuedVariables.resize(m_NumStates, false);

    for (TTransList::const_iterator j = m_TransByCat[eDeterministic].begin();
         j != m_TransByCat[eDeterministic].end();  ++j) {
        for (unsigned int i = 0;  i < m_Nu[*j].size();  ++i) {
            m_RealValuedVariables[m_Nu[*j][i].m_State] = true;
        }
    }
                    
}

/*---------------------------------------------------------------------------*/
// PRE : list of critical transitions & their total rate
// POST: one picked according to probability
unsigned int CStochasticEqns::x_PickCritical(double critRate) const {
    double r = runif(0,1);
    double d = 0;
    TTransList::const_iterator j = m_TransByCat[eCritical].begin();
    while (j != m_TransByCat[eCritical].end()) {
        d += m_Rates[*j]/critRate;
        if (d > r) {
            break;
        }
        ++j;
    }
    if (!(d >= r)) { throwError("logic error at line " << __LINE__) }
    return *j;
}

/*---------------------------------------------------------------------------*/
// PRE : time period to step; whether to clamp variables at 0
// POST: all determinisitic transitions updated by the expected amount
// (i.e. Euler method); if clamping, then negative variables set to 0.
void CStochasticEqns::x_AdvanceDeterministic(double deltaT, bool clamp) {
    for (TTransList::const_iterator j = m_TransByCat[eDeterministic].begin();
         j != m_TransByCat[eDeterministic].end();  ++j) {
        for (unsigned int i = 0;  i < m_Nu[*j].size();  ++i) {
            m_X[m_Nu[*j][i].m_State] += m_Nu[*j][i].m_Mag * m_Rates[*j] *
                deltaT;
            //clamp at zero if specified
            if (clamp  &&  m_X[m_Nu[*j][i].m_State] < 0) {
                m_X[m_Nu[*j][i].m_State] = 0;
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
// PRE : simulation end time; **transition rates already updated**
// POST: id of transition taken (if none, then -1) & time series updated.
void CStochasticEqns::x_SingleStepExact(double tf) {
    m_LastTransition = -1;
    double stochRate = 0;
    double detRate = 0;
    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] != eDeterministic) {
            stochRate += m_Rates[j];
        } else {
            detRate += m_Rates[j];
        }
    }

    double tau = stochRate > 0 ? rexp(1./stochRate) :
        detRate > 0 ? 1./detRate : tf - *m_T;
    if (stochRate == 0  ||  tau > tf - *m_T) {
        tau = tf - *m_T; // step is off end so just advance time
    } else {
        double r = runif(0,1);
        double d = 0;
        unsigned int j = 0;
        for (;  j < m_Nu.size()  &&  d < r;  ++j) {
            if (m_TransCats[j] != eDeterministic) {
                d += m_Rates[j]/stochRate;
            }
        }
        if (!(d >= r)) { throwError("logic error at line " << __LINE__) }
        --j;

        //take transition "j"
        if (m_VerboseTracing >= 1) {
            REprintf("%f: taking transition #%i\n", *m_T, j+1);
        }
        for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
            m_X[m_Nu[j][i].m_State] += m_Nu[j][i].m_Mag;
        }
        m_LastTransition = j;
    }

    //clamp deterministic at 0, assuming that it is unreasonable to
    //take a smaller step then exact.
    x_AdvanceDeterministic(tau, true);
    *m_T += tau;
    m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
}

/*---------------------------------------------------------------------------*/
// PRE : tau value to use for step, list of "critical" transitions
// POST: IMPLICIT tau step taken (m_X updated if so) (or overflow
// error thrown if tau was too big)
// NOTE: See equation (7) in Cao et al. (2007)
void CStochasticEqns::x_SingleStepITL(double tau) {
    if (m_VerboseTracing >= 1) {
        REprintf("%f: taking implicit step of tau = %f\n", *m_T, tau);
    }
    if (!m_RateJacobianFunc) { throwError("logic error at line " << __LINE__) }
    double *origX = new double[m_NumStates];
    double *origRates = new double[m_NumStates];
    memcpy(origX, m_X, sizeof(double)*m_NumStates);
    memcpy(origRates, m_Rates, sizeof(double)*m_NumStates);

    if (debug) {
        cerr << " origX: ";
        for (unsigned int i =0; i < m_NumStates;  ++i) {
            cerr << origX[i] << "\t";
        }
        cerr << endl;
    }

    // draw (stochastic) number of times each transition will occur
    vector<int> numTransitions(m_Nu.size(), 0);
    
    for (TTransList::const_iterator j = m_TransByCat[eNormal].begin();
         j != m_TransByCat[eNormal].end();  ++j) {
        if (m_Rates[*j]*tau > 1e8) {
            //for high rate, use normal to approx poisson.
            //should basically never yield negative, but just to
            //be sure, cap at 0
            numTransitions[*j] = max(0.,floor(rnorm(m_Rates[*j]*tau,
                                                    sqrt(m_Rates[*j]*tau))));
        } else {
            numTransitions[*j] = rpois(m_Rates[*j]*tau);
            //cerr << "nt[" << *j << "] " << numTransitions[*j] << endl;
        }
    }

    // Calculate equation (7) terms not involving x[t+tau] and call this alpha:
    //   alpha = x + nu.(P - tau/2 R(x))
    // Also initialize iterative search for x[t+tau] at expectation (reset m_X)
    double* alpha = new double[m_NumStates];
    memcpy(alpha, m_X, sizeof(double)*m_NumStates);
    for (TTransList::const_iterator j = m_TransByCat[eNormal].begin();
         j != m_TransByCat[eNormal].end();  ++j) {
        for (unsigned int k = 0;  k < m_Nu[*j].size();  ++k) {
            alpha[m_Nu[*j][k].m_State] += m_Nu[*j][k].m_Mag * 
                (numTransitions[*j] - (tau/2)*m_Rates[*j]);
            //reset m_X to expectation as our initial guess
            m_X[m_Nu[*j][k].m_State] += m_Nu[*j][k].m_Mag *
                (tau/2)*m_Rates[*j];
        }
    }
    //expectations may send states negative; clamp!
    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            m_X[i] = 0;
        }
    }

    if (debug) {
        cerr << " alpha:";
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            cerr << " " << alpha[i];
        }
        cerr << endl;
        cerr << "it " << 0 << " newX: ";
        for (unsigned int i =0; i < m_NumStates;  ++i) {
            cerr << m_X[i] << "\t";
        }
        cerr << endl;
    }

    //a few variables needed by LAPACK
    int N = m_NumStates;
    int nrhs = 1;
    int *ipiv = new int[m_NumStates];
    int info;
    double *matrixA = new double[m_NumStates*m_NumStates];
    double *matrixB = new double[m_NumStates];

    
    //Use Newton's method to solve implicit equation:
    //  Let Y = x[t+tau]
    //  F(Y) = Y - alpha - nu.((tau/2)*R(Y))
    //Solve Jacobian(F(Y0)) Y1 = -F(Y0) for Y1 to iteratively approach solution
    //This eqn expands to (I - nu.((tau/2)Jacobian(R(Y0)))) Y1 = -F(Y0) where
    //the Jacobian of rates is supplied by the user.  The term to the
    //left of Y1 is called matrix A by LAPACK and the term on right is
    //matrix B.
    //
    //Perhaps should adjust max # of iterations..
    bool converged = false;
    unsigned int c = 0;
    while (++c <= 20  &&  !converged) {
        // Check to make sure we haven't taken too big a step --
        // i.e. no state variables should go negative
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            if (m_X[i] < 0) {
                delete[] origRates;
                delete[] alpha;
                delete[] ipiv;
                delete[] matrixA;
                delete[] matrixB;
                memcpy(m_X, origX, sizeof(double)*m_NumStates);
                delete[] origX;
                throw overflow_error("tau too big");
            }
        }

        // define matrix A
        double* rateJacobian = x_CalcJacobian();
        memset(matrixA, 0, m_NumStates*m_NumStates*sizeof(double));
        for (TTransList::const_iterator j =
                 m_TransByCat[eNormal].begin();
             j != m_TransByCat[eNormal].end();  ++j) {
            for (unsigned int k = 0;  k < m_Nu[*j].size();  ++k) {
                for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                    //R stores matrices column-wise
                    //LAPACK stores matrices row-wise
                    matrixA[i*m_NumStates + m_Nu[*j][k].m_State] +=
                        m_Nu[*j][k].m_Mag * rateJacobian[(*j)*m_NumStates + i];
                }
            }
        }
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            for (unsigned int i2 = 0;  i2 < m_NumStates;  ++i2) {
                matrixA[i*m_NumStates + i2] *= -tau/2;
            }
            matrixA[i*m_NumStates + i] += 1;
        }

        // define matrix B
        // m_X is now our proposed x[t+tau].  Note that m_X has changed
        // even in our first iteration (initialized to expected value).
        x_UpdateRates();
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            matrixB[i] = alpha[i] - m_X[i];
        }
        for (TTransList::const_iterator j =
                 m_TransByCat[eNormal].begin();
             j != m_TransByCat[eNormal].end();  ++j) {
            for (unsigned int k = 0;  k < m_Nu[*j].size();  ++k) {
                matrixB[m_Nu[*j][k].m_State] += 
                    m_Nu[*j][k].m_Mag * (tau/2) * m_Rates[*j];
            }
        }


    if (debug) {
        cerr << "A:" << endl;
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            for (unsigned int i2 = 0;  i2 < m_NumStates;  ++i2) {
                cerr << matrixA[i2*m_NumStates + i] << "\t";
            }
            cerr << endl;
        }

        cerr << "B:" << endl;
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            cerr << matrixB[i] << "\t";
        }
        cerr << endl;

        cerr << "a:" << endl;
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            cerr << m_Rates[j] << "\t";
        }
        cerr << endl;
    }

        //solve linear eqn
        F77_NAME(dgesv)(&N, &nrhs, matrixA, &N, ipiv, matrixB, &N, &info); 
        if (info != 0) {
            warning("warning: lapack ran into trouble solving implicit equation");
            break;
        }
        //matrixB now contains solution (change in X)
        double normDelta = 0, normX = 0;
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            m_X[i] += matrixB[i];
            normDelta += matrixB[i]*matrixB[i];
            normX += m_X[i] * m_X[i];
        }
        //cerr << "\tNorms: " << normDelta << "\t" << normX << endl;
        converged = (normDelta < normX * m_ITLConvergenceTol);
        if (debug) {
            /*
            cerr << "Delta: ";
            for (unsigned int i =0; i < m_NumStates;  ++i) {
                cerr << matrixB[i] << "\t";
            }
            cerr << endl;
            */
            cerr << "it " << c << " newX: ";
            for (unsigned int i =0; i < m_NumStates;  ++i) {
                cerr << m_X[i] << "\t";
            }
            cerr << endl;
            /*
            x_UpdateRates();
            double t[m_NumStates];
            memcpy(t, alpha, sizeof(double)*m_NumStates);
            for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
                if (m_TransCats[j] == eNoncritical) {
                    for (unsigned int k = 0;  k < m_Nu[j].size();  ++k) {
                        t[m_Nu[j][k].m_State] +=
                            m_Nu[j][k].m_Mag * (tau/2) * m_Rates[j];
                    }
                }
            }
            cerr << "     newX: ";
            for (unsigned int i =0; i < m_NumStates;  ++i) {
                cerr << t[i] << "\t";
            }
            cerr << endl;
            */
        }

    } // end of iterating for Newton's method
    if (!converged) {
        warning("ITL solution did not converge!");
    }

    //restore original rates to execute deterministic transitions
    memcpy(m_Rates, origRates, sizeof(double)*m_NumStates);
    x_AdvanceDeterministic(tau);

    delete[] origRates;
    delete[] alpha;
    delete[] ipiv;
    delete[] matrixA;
    delete[] matrixB;

    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            memcpy(m_X, origX, sizeof(double)*m_NumStates);
            delete[] origX;
            throw overflow_error("tau too big");
        }
        if (!m_RealValuedVariables[i]) {
            m_X[i] = round(m_X[i]);
        }
    }
    delete[] origX;
    *m_T += tau;
}

/*---------------------------------------------------------------------------*/
// PRE : tau value to use for step, list of "critical" transitions
// POST: EXPLICIT tau step taken (m_X updated if so) (or overflow
// error thrown if tau was too big)
void CStochasticEqns::x_SingleStepETL(double tau) {
    if (m_VerboseTracing >= 1) {
        REprintf("%f: taking explicit step of tau = %f\n", *m_T, tau);
    }
    if (m_VerboseTracing >= 2) {
        REprintf("%f:    ", *m_T);
    }
    double *origX = new double[m_NumStates];
    memcpy(origX, m_X, sizeof(double)*m_NumStates);
    for (TTransList::const_iterator j = m_TransByCat[eNormal].begin();
         j != m_TransByCat[eNormal].end();  ++j) {
        double k;
        if (m_Rates[*j]*tau > 1e8) {
            //for high rate, use normal to approx poisson.
            //should basically never yield negative, but just to
            //be sure, bound at 0
            k = max(0.,floor(rnorm(m_Rates[*j]*tau, sqrt(m_Rates[*j]*tau))));
        } else {
            k = rpois(m_Rates[*j]*tau);
        }
        if (k > 0) {
            if (m_VerboseTracing >= 2) {
                REprintf("%fx#%i ", k, *j);
            }
            for (unsigned int i = 0;  i < m_Nu[*j].size();  ++i) {
                m_X[m_Nu[*j][i].m_State] +=  k * m_Nu[*j][i].m_Mag;
            }
        }
    }
    if (m_VerboseTracing >= 2) {
        REprintf("\n");
    }
    x_AdvanceDeterministic(tau);

    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            memcpy(m_X, origX, sizeof(double)*m_NumStates);
            delete[] origX;
            throw overflow_error("tau too big");
        }
    }

    *m_T += tau;
    delete[] origX;
}

/*---------------------------------------------------------------------------*/
// PRE : time at which to end simulation; **transition rates already updated**
// POST: single adaptive tau leaping step taken & time series updated.
// Implemented from Cao Y, Gillespie DT, Petzold LR. The Journal of Chemical
// Physics (2007).
void CStochasticEqns::x_SingleStepATL(double tf) {
    m_LastTransition = -1;
    EStepType stepType;

    //identify "critical" transitions
    double criticalRate = 0;
    double noncritRate = 0;
    {
        for (TTransList::const_iterator j =
                 m_TransByCat[eDeterministic].begin();
             j != m_TransByCat[eDeterministic].end();  ++j) {
            noncritRate += m_Rates[*j];
        }
        for (TTransList::const_iterator j =
                 m_TransByCat[eHalting].begin();
             j != m_TransByCat[eHalting].end();  ++j) {
            criticalRate += m_Rates[*j];
        }

        //reset (lop off) all non-halting criticals
        m_TransByCat[eCritical].resize(m_TransByCat[eHalting].size());
        m_TransByCat[eNormal].clear();
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] != eNormal) {
                continue;
            }
            unsigned int minTimes = numeric_limits<unsigned int>::max();
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                if (m_Nu[j][i].m_Mag < 0  &&
                    m_X[m_Nu[j][i].m_State]/abs(m_Nu[j][i].m_Mag) < minTimes) {
                    minTimes = m_X[m_Nu[j][i].m_State]/abs(m_Nu[j][i].m_Mag);
                }
            }
            if (minTimes < m_Ncritical) {
                criticalRate += m_Rates[j];
                m_TransByCat[eCritical].push_back(j);
            } else {
                noncritRate += m_Rates[j];
                m_TransByCat[eNormal].push_back(j);
            }
        }
    }

    if (debug) {
        cerr << "critical rate: " << criticalRate << "\t" << "noncrit rate: " << noncritRate << endl;
    }
    if (criticalRate + noncritRate == 0) {
        *m_T = tf;//numeric_limits<double>::infinity();
        m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
        return;
    }
    if (!R_finite(criticalRate + noncritRate)) {
        throwEarlyExit("Infinite transition rate at time " << *m_T);
    }

    // calc explicit & implicit taus
    double tau1, tau2;
    double tauEx = x_TauEx();
    double tauIm = x_TauIm();
    if (debug) {
        cerr << "tauEx: " << tauEx << "  tauIm:" << tauIm << endl;
    }
    if (tauEx*m_Nstiff < tauIm) {
        stepType = eImplicit;
        tau1 = tauIm;
    } else {
        stepType = eExplicit;
        tau1 = tauEx;
    }
    if (tau1 > tf - *m_T) { //cap at the final simulation time
        tau1 = tf - *m_T;
    }
    if (tau1 > m_MaxTau) {
        tau1 = m_MaxTauFunc ? min(tau1, x_CalcUserMaxTau()) : m_MaxTau;
        if (debug) {
            cerr << "maxtau: " << tau1 << " (" <<
                (m_MaxTauFunc ? x_CalcUserMaxTau() : m_MaxTau) << ")" << endl;
        }
    }

    bool tauTooBig;
    do {
        tauTooBig = false;
        if (!(tau1 > 0)) { throwError("logic error at line " << __LINE__) }
        if (tau1 < m_ExactThreshold / (criticalRate + noncritRate)) {
            if (debug) {
                cerr << "Taking exact steps.. (tau1 = " << tau1 << ")" << endl;
            }
            stepType = eExact;
            for (unsigned int i = 0;
                 i < m_NumExactSteps[m_PrevStepType]  &&  *m_T < tf;  ++i) {
                if (i > 0) {
                    x_UpdateRates();
                }
                x_SingleStepExact(tf);
                if (m_VerboseTracing >= 2) {
                    REprintf("%f -- ", *m_T);
                    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                        REprintf("%f ", m_X[i]);
                    }
                    REprintf("\n");
                }
                if (m_LastTransition >= 0  &&
                    m_TransCats[m_LastTransition] == eHalting) {
                    return;
                }
            }
        } else {
            try { //catch exception if tauTooBig
                tau2 = (criticalRate == 0) ? numeric_limits<double>::infinity():
                    rexp(1./criticalRate);
                if (stepType == eExplicit  ||
                    (tau1 > tau2  &&  stepType == eImplicit && tau2 <= tauEx)) {
                    if (debug) {
                        cerr << "going explicit w/ tau = " << min(tau1, tau2)
                             << endl;
                    }
                    x_SingleStepETL(min(tau1, tau2));
                } else {
                    if (debug) {
                        cerr << "going implicit w/ tau = " << tau1 << endl;
                    }
                    x_SingleStepITL(tau1);
                }
                if (tau1 > tau2) { //pick one critical transition
                    unsigned int j = x_PickCritical(criticalRate);
                    m_LastTransition = j;
                    if (debug) {
                        cerr << "hittin' the critical (" << j << ")" << endl;
                    }
                    if (m_VerboseTracing >= 1) {
                        REprintf("%f:    executing critical transition #%i\n",
                                 *m_T, j+1);
                    }
                    for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                        m_X[m_Nu[j][i].m_State] +=  m_Nu[j][i].m_Mag;
                        if (m_X[m_Nu[j][i].m_State] < 0) {
                            throwError("variable " << m_Nu[j][i].m_State+1 <<
                                       " went negative after executing "
                                       "transition " << j+1 << ".  Most likely "
                                       "either your rate calculation or "
                                       "transition matrix is flawed.");
                        }
                    }
                }

                m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
                if (m_VerboseTracing >= 2) {
                    REprintf("%f -- ", *m_T);
                    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                        REprintf("%f ", m_X[i]);
                    }
                    REprintf("\n");
                }
            } catch (overflow_error&) { //i.e. tauTooBig exception
                tauTooBig = true;
                if (m_VerboseTracing >= 1) {
                    REprintf("%f:    tau too big; cutting in half\n", *m_T);
                }
                tau1 /= 2;
            }
        }
    } while (tauTooBig);

    m_PrevStepType = stepType;
}

/*---------------------------------------------------------------------------*/
// Exported C entrypoints for calling from R
/*
extern "C" {
    SEXP simAdaptiveTau(SEXP s_x0, SEXP s_nu, SEXP s_f, SEXP s_fJacob,
                        SEXP s_params, SEXP s_tf,
                        SEXP s_deterministic, SEXP s_halting,
                        SEXP s_changebound,
                        SEXP s_tlparams, SEXP s_fMaxtau) {
        try{
        if (!isVector(s_x0)  ||  !(isReal(s_x0)  ||  isInteger(s_x0))) {
            error("invalid vector of initial values");
        }
        if (!isVectorList(s_nu)  &&
            (!isMatrix(s_nu)  ||
             INTEGER(getAttrib(s_nu, R_DimSymbol))[0] != length(s_x0))) {
            error("invalid transition specification");
        }
        if (!isFunction(s_f)) {
            error("invalid rate function");
        }
        if (!isNull(s_fJacob)  &&  !isFunction(s_fJacob)) {
            error("invalid Jacobian function");
        }
        if (!(isReal(s_tf)  ||  isInteger(s_tf))  ||  length(s_tf) != 1) {
            error("invalid final time");
        }
        if (!isVector(s_changebound)  ||  !isReal(s_changebound)  ||
            length(s_changebound) != length(s_x0)) {
            error("invalid relratechange");
        }
        if (!isNull(s_tlparams)  &&  !isVector(s_tlparams)) {
            error("tl.params must be a list");
        }
        if (!isNull(s_fMaxtau)  &&  !isFunction(s_fMaxtau)) {
            error("invalid maxTau function");
        }

        CStochasticEqns eqns(s_x0, s_nu,
                             s_f, s_fJacob, s_params, REAL(s_changebound),
                             s_fMaxtau, s_deterministic, s_halting);
        if (!isNull(s_tlparams)) {
            eqns.SetTLParams(s_tlparams);
        }
        try {
            eqns.EvaluateATLUntil(REAL(coerceVector(s_tf, REALSXP))[0]);
        } catch (CEarlyExit &e) {
            warning(e.what());
        }
        return eqns.GetResult();
        } catch (exception &e) {
            error(e.what());
            return R_NilValue;
        }
    }
    
    //-----------------------------------------------------------------------

    SEXP simExact(SEXP s_x0, SEXP s_nu, SEXP s_f, SEXP s_params, SEXP s_tf) {
        try {
        if (!isVector(s_x0)  ||  !(isReal(s_x0)  ||  isInteger(s_x0))) {
            error("invalid vector of initial values");
        }
        if (!isVectorList(s_nu)  &&
            (!isMatrix(s_nu)  ||
             INTEGER(getAttrib(s_nu, R_DimSymbol))[0] != length(s_x0))) {
            error("invalid transition specification");
        }
        if (!isFunction(s_f)) {
            error("invalid rate function");
        }
        if (!(isReal(s_tf)  ||  isInteger(s_tf))  ||  length(s_tf) != 1) {
            error("invalid final time");
        }

        CStochasticEqns eqns(s_x0, s_nu,
                             s_f, NULL, s_params, NULL, NULL,
                             R_NilValue, R_NilValue);
        try {
            eqns.EvaluateExactUntil(REAL(coerceVector(s_tf, REALSXP))[0]);
        } catch (CEarlyExit &e) {
            warning(e.what());
        }
        return eqns.GetResult();
        } catch (exception &e) {
            error(e.what());
            return R_NilValue;
        }
    }
}
*/
/********/
// Exported C entrypoints for calling from R for SEIR model and SIR model.
extern "C" {
   // Doing the event that occurs
   // for the SEIR model
    int* seirtransitions(int nbVilles){
        //Initializing
        //4 for number of states, 8 for number of events
        int *transitions = new int[nbVilles*4*nbVilles*8];
        for(int i=0; i<nbVilles*4*nbVilles*8;i++) transitions[i] =0;
        for (unsigned int i = 0;  i < nbVilles;  ++i) {
            //susceptible birth; S<- S+1
            transitions[0*nbVilles*4*nbVilles + i*4*nbVilles+0*nbVilles+i] = 1;
            //susceptible death : S <-S-1
            transitions[1*nbVilles*4*nbVilles + i*4*nbVilles+0*nbVilles+i] = -1;
            //exposed death : E <- E-1
            transitions[2*nbVilles*4*nbVilles + i*4*nbVilles+1*nbVilles+i] = -1;
            //infectious death : I<- I-1
            transitions[3*nbVilles*4*nbVilles + i*4*nbVilles+2*nbVilles+i] = -1;
            //recovered death : R<-R-1
            transitions[4*nbVilles*4*nbVilles + i*4*nbVilles+3*nbVilles+i] = -1;
            // infection : S <-S-1; E<- E+1
            transitions[5*nbVilles*4*nbVilles + i*4*nbVilles+0*nbVilles+i] = -1;
            transitions[5*nbVilles*4*nbVilles + i*4*nbVilles+1*nbVilles+i] = 1;
            // becoming infectious: E<-E-1; I<-I+1
            transitions[6*nbVilles*4*nbVilles + i*4*nbVilles+1*nbVilles+i] = -1;
            transitions[6*nbVilles*4*nbVilles + i*4*nbVilles+2*nbVilles+i] = 1;
            // recovery : I<- I-1; R <-R+1
            transitions[7*nbVilles*4*nbVilles + i*4*nbVilles+2*nbVilles+i] = -1;
            transitions[7*nbVilles*4*nbVilles + i*4*nbVilles+3*nbVilles+i] = 1;
        }
        // Printing the result on the screen to check the correction of the results
         /*
	        for(int i=0; i<nbVilles*4*nbVilles*8;i++)
	            cout<<"i= "<<i<<"  trans="<<transitions[i]<<"   ";
	        cout<<endl;
        */
        return(transitions);
}

/***/
    // Transitions for the SIR model
    int* sirtransitions(int nbVilles){
        //Initializing
        // 3 for number of states, 6 for number of events
        int *transitions = new int[nbVilles*3*nbVilles*6];
        for(int i=0; i<nbVilles*3*nbVilles*6;i++) transitions[i] =0;

        // event that occurs
        for (unsigned int i = 0;  i < nbVilles;  ++i) {
            //susceptible birth; S<- S+1
            transitions[0*nbVilles*3*nbVilles + i*3*nbVilles+0*nbVilles+i] = 1;
            //susceptible death : S <-S-1
            transitions[1*nbVilles*3*nbVilles + i*3*nbVilles+0*nbVilles+i] = -1;
            //infectious death : I<- I-1
            transitions[2*nbVilles*3*nbVilles + i*3*nbVilles+1*nbVilles+i] = -1;
            //recovered death : R<-R-1
            transitions[3*nbVilles*3*nbVilles + i*3*nbVilles+2*nbVilles+i] = -1;
            // infection : S <-S-1; I<- I+1
            transitions[4*nbVilles*3*nbVilles + i*3*nbVilles+0*nbVilles+i] = -1;
            transitions[4*nbVilles*3*nbVilles + i*3*nbVilles+1*nbVilles+i] = 1;
            // recovery : I<- I-1; R <-R+1
            transitions[5*nbVilles*3*nbVilles + i*3*nbVilles+1*nbVilles+i] = -1;
            transitions[5*nbVilles*3*nbVilles + i*3*nbVilles+2*nbVilles+i] = 1;
        }
        // Printing the result on the screen to check the correction of the results
        /*
         *             for(int i=0; i<nbVilles*3*nbVilles*6;i++)
         *                 cout<<"i= "<<i<<"  trans="<<transitions[i]<<"   ";
         *             cout<<endl;
         */
        return(transitions);
  }
    /****/
    /*
        Doing stochastic simulation for the SEIR and SIR models
        by using the "adaptive tau-leaping" algorithm described by Cao Y, Gillespie DT, Petzold LR. The Journal of Chemical Physics (2007)
        and by integrating the R and the C++.

        In order to the values of parameters and of variables from R.
          Step 1: we get these values from R that are in the R type.
          Step 2: we convert the values of the R type to these of the C++ type
        */
      SEXP ssesAdaptiveTauRETRY(SEXP s_x0, // initial values of variables
                           SEXP s_f, SEXP s_fJacob,
                           // values of parameters
                           SEXP sigma, SEXP gamma,SEXP mu,SEXP seed,
                           SEXP nbCONTACT0, SEXP nbCONTACT1, SEXP nbMulCONTACT,
                           SEXP phiPHASE, SEXP probVISITER, SEXP probINFECTER, SEXP nbVilles,
                           SEXP periDISE,SEXP tmax,// time of simulation
                           SEXP s_deterministic, SEXP s_changebound,
                           SEXP s_tlparams, SEXP s_fMaxtau) {
          try{
              if (!isVector(s_x0)  ||  !isReal(s_x0)) {
               error("invalid vector of initial values");
                }
               if (!isNull(s_fJacob)  &&  !isFunction(s_fJacob)) {
                   error("invalid Jacobian function");
               }
               if (length(tmax) != 1) {
                   error("invalid final time");
               }

               if (!isVector(s_deterministic)  ||  !isLogical(s_deterministic)) {
                  error("invalid deterministic parameter -- must be logical vector");
               }

               if (!isVector(s_changebound)  ||  !isReal(s_changebound)  || length(s_changebound) != length(s_x0)) {
                   error("invalid relratechange");
               }

               if (!isNull(s_tlparams)  &&  !isVector(s_tlparams)) {
                   error("tl.params must be a list");
               }

               if (!isNull(s_fMaxtau)  &&  !isFunction(s_fMaxtau)) {
                   error("invalid maxTau function");
               }

               //Getting the values of parameters and variables
               //number of subpopulations
      //Number of cities in the metapopulation
          if(!isReal(nbVilles) || length(nbVilles)!=1)
              error("invalid number of cities.");
          int cnbVilles = INTEGER_VALUE(nbVilles);

          //-gamma 0.2  //gamma
          if (!isReal(gamma) || length(gamma)!=1) {
              error("invalid parameter gamma");
          }
          double cgamma = NUMERIC_VALUE(gamma);

          //-mu 0.000039   //mu
          if (!isReal(mu) || length(mu) !=1) {
              error("invalid parameter mu");
          }
          double cmu = NUMERIC_VALUE(mu);

          //-seed 10    //number 'seed'
          if (!isReal(seed)  ||  length(seed) != 1) {
              error("invalid value of seed,it should be a real number");
          }
          double cseed=NUMERIC_VALUE(seed);

          //-nbCONTACT0 300
          if (!isReal(nbCONTACT0)  ||  length(nbCONTACT0) != 1) {
              error("invalid value of nbCONTACT0,it should be a real number");
          }
          double cnbCONTACT0=NUMERIC_VALUE(nbCONTACT0);

          //-nbCONTACT1 0.1
          if (!isReal(nbCONTACT1)  ||  length(nbCONTACT1) != 1) {
              error("invalid value of nbCONTACT1,it should be a real number");
          }
          double cnbCONTACT1=NUMERIC_VALUE(nbCONTACT1);

          //-nbMulCONTACT 1
          if (!isReal(nbMulCONTACT)  ||  length(nbMulCONTACT) != 1) {
              error("invalid value of nbMulCONTACT,it should be a real number");
          }
          double cnbMulCONTACT=NUMERIC_VALUE(nbMulCONTACT);

          //-phiPHASE  is a vector of n elements
           if (!isVector(phiPHASE)  ||  !isReal(phiPHASE)) {
              error("invalid vector of phase (phiPHASE)");
          }
          double *rphiPHASE = REAL(phiPHASE);
          double* cphiPHASE; create1D(cnbVilles,0.0,cphiPHASE);
          for(int i=0; i<cnbVilles; i++) cphiPHASE[i] = rphiPHASE[i];

          //-probVISITER 0.01
          if (!isReal(probVISITER)  ||  length(probVISITER) != 1) {
              error("invalid value of probVISITER,it should be a real number");
          }
          double cprobVISITER=NUMERIC_VALUE(probVISITER);

          //-probINFECTER 0.01
          if (!isReal(probINFECTER)  ||  length(probINFECTER) != 1) {
              error("invalid value of probINFECTER,it should be a real number");
          }
          double cprobINFECTER=NUMERIC_VALUE(probINFECTER);


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

          // Disease period
          if (!isReal(periDISE)  ||  length(periDISE) != 1) {
              error("invalid value of period, it should be a real number");
          }
          double cperiDISE=NUMERIC_VALUE(periDISE);

          // for metapopulation
          double ** arr_probVISITER = calculerProbVISITER(cnbVilles,cprobVISITER);

          // xay dung ma tran INFECTER
          double **arr_probINFECTER; create2D(cnbVilles,cnbVilles,cprobINFECTER,arr_probINFECTER);

          double **arrNbCONTACT0 = calculerNbCONTACT(cnbVilles,cnbCONTACT0,1);
          double **arrNbCONTACT1 = calculerNbCONTACT(cnbVilles,cnbCONTACT1,cnbMulCONTACT);
         // xay dung ma tran VENIR de la ville k
        // initialiser les valeurs des variables
           //separate initial values of variables
          unsigned int nbTotalVar= length(s_x0); //total number of states
          double *rs_x0 = REAL(s_x0);
          double* cs_x0; create1D(nbTotalVar,0.0,cs_x0);
          for(int i=0; i<nbTotalVar; i++) { cs_x0[i] = rs_x0[i]; cout <<"iii =" <<cs_x0[i]<<endl;}
          int nvar= 6;
          unsigned long **valSEIR0;
          valSEIR0 = new unsigned long*[cnbVilles];
          for (int i = 0; i <  cnbVilles; i++)
          {
              valSEIR0[i] = new unsigned long[nvar];
              valSEIR0[i][0]= 0;
              valSEIR0[i][iS]= (unsigned long)cs_x0[0*cnbVilles+i];//iS
              valSEIR0[i][iE]= (unsigned long)cs_x0[1*cnbVilles+i];//iE
              valSEIR0[i][iI]= (unsigned long)cs_x0[2*cnbVilles+i];//iI
              valSEIR0[i][iR]= (unsigned long)cs_x0[3*cnbVilles+i];//iR
              valSEIR0[i][iN]= valSEIR0[i][iS]+ valSEIR0[i][iE]+ valSEIR0[i][iI]+valSEIR0[i][iR];
          }

          cout<<"valSEIR0"<<endl;
          for(int i=0; i<cnbVilles; i++){
              for (int j=0; j<nvar; j++)
                  cout<<"    " << valSEIR0[i][j];
              cout<<endl;
          }
  cout<<"valSEIRAZZZZZZZ"<<endl;

          double **arr_probVENIRk;  create2D(cnbVilles,cnbVilles,0, arr_probVENIRk);
          calculerProbVENIRk(cnbVilles,valSEIR0,arr_probVISITER,arr_probVENIRk);
  cout<<"valSEIRAZZZZZZZ"<<endl;

          // vector of transitions (event occurs)
          int* transitions;
          unsigned int nbrate = 0;
          // check here is SEIR model or SIR model
          //-sigma 0.125   //sigma
          if (!isReal(sigma) || length(sigma)!=1) {
              error("invalid parameter sigma");
          }
          double csigma = NUMERIC_VALUE(sigma);
  cout<<"valSEIRAZZZZZZZ"<<endl;
        bool flagSEIR=TRUE;
        if(csigma==INFINITY) flagSEIR=FALSE;
  cout<<"valSEIRAZZZZZZZ"<<endl;
         // Doing simulation for SEIR/SIR models
        if(flagSEIR){//SEIR model
            transitions = seirtransitions(cnbVilles);
            nbrate=cnbVilles*8; //number of events for all subpopulations
        }
        else{//SIR model
            transitions = sirtransitions(cnbVilles);
            nbrate=cnbVilles*6; //number of events for all subpopulations
        }
   cout<<"valSEIRAZZZZZZZ"<<endl;
        cout<<"aaaaa "<<cnbVilles<<"  "<<csigma<<"  "<<cgamma<<"  "<<cmu<<"   "<<cphiPHASE<<"  "<<ctmax<<"  "<<cperiDISE<<endl;
        cout<<"valSEIRAAAAAA"<<endl;
          for(int i=0; i<cnbVilles; i++){
              for (int j=0; j<cnbVilles; j++){
                  cout<<"    " << arr_probVISITER[i][j];
       cout<<"    " << arrNbCONTACT0[i][j];
       cout<<"    " << arrNbCONTACT1[i][j];
   cout<<"    " << arr_probINFECTER[i][j];
   cout<<"    " << arr_probVENIRk[i][j];
              }
              cout<<endl;
          }
  cout<<"valSEIRAZZZZZZZ"<<endl;
         // stochastic simulation
       CStochasticEqns eqns(s_x0,transitions,nbrate,
                            s_f, s_fJacob,
                            //parameters
                            cnbVilles,csigma,cgamma,cmu,arr_probVISITER,arrNbCONTACT0,arrNbCONTACT1,
                            arr_probINFECTER,arr_probVENIRk,cphiPHASE,ctmax,cperiDISE,
                            //
                            REAL(s_changebound),
                            s_fMaxtau, s_deterministic);
       delete []transitions;
       for (int i = 0; i <  cnbVilles; i++){
           delete []arr_probVISITER[i];
           delete []arr_probINFECTER[i];
           delete []arrNbCONTACT0[i];
           delete []arrNbCONTACT1[i];
       }
       delete []arr_probVISITER; delete[]arr_probINFECTER;delete []arrNbCONTACT0; delete[]arrNbCONTACT1;

       if (!isNull(s_tlparams)) {
           eqns.SetTLParams(s_tlparams);
       }
       eqns.EvaluateATLUntil(REAL(tmax)[0]);
       return eqns.GetTimeSeriesSEXP();
          } catch (exception &e) {
              error(e.what());
              return R_NilValue;
          }
      }
/*****/
}