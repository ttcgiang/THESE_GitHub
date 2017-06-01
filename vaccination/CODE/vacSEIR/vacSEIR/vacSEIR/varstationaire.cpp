#include "initialFUNC.h"
#include "varstationaire.h"
/**********************************/
double** createUnifNumb(int nbVilles, double seuil, double min, double max ){
    double**res;create2D(nbVilles,nbVilles,0.0,res);
    std::random_device rd;
    std::mt19937 gen(rd());
    double tp=1000000;
    double sumtp=0.0;

    std::uniform_real_distribution<> dis(min, max);
    for (int i = 0; i < nbVilles; i++) {
          sumtp=0.0;
        for(int j=0; j<nbVilles;j++){
               // cout<<"i="<<i<<"  j="<<j<<endl;
//                 cout<<"seuil="<<seuil<<endl;
            if(j != i){            //voir des voisins
                do{
                    tp=dis(gen);
                }
                while(tp>seuil);
                sumtp +=tp;
                res[i][j]=tp;
                //  std::cout << tp << ' ';
            }
        }
         res[i][i]=1.0 - sumtp;

    }

    return (res);
}
/******/
// Tich hai ma tran
double** exponentMatrix(double** matrix,int nrow, int power){
    double**x; double **y; double**c;
    x = new double*[nrow]; y = new double*[nrow]; c = new double*[nrow];
    for (int i=0; i<nrow; i++)
    {
        x[i] = new double[nrow];
        y[i] = new double[nrow];
        c[i] = new double[nrow];
        for(int j=0; j<nrow; j++){
             x[i][j] = matrix[i][j];
             y[i][j] = matrix[i][j];
             c[i][j] = 0.0;
         }
    }


    //multiplier des matrices
    for(int tp=0;tp<power;tp++){
        for(int i=0;i<nrow;i++)
        {
            for(int j=0;j<nrow;j++)
            {
                c[i][j]=0;
                for(int k=0;k<nrow;k++)
                {
                    c[i][j]=c[i][j]+y[i][k]*x[k][j];
                }
            }
        }
        //copier la matrice
  //      cout<<"tp ="<<tp<<endl;
     //   cout<<"matrix y"<<endl;
        for(int i=0;i<nrow;i++)
        {
            for(int j=0;j<nrow;j++)
            {
                y[i][j]=c[i][j];
         //        cout<<"\t"<<y[i][j];
            }
         //   cout<<"\n\n";
        }
    }

    cout<<"\n-----------------------------------------------------------\n"<<endl;
    cout<<"\n\n Matrice stationaire :\n\n";

    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<nrow;j++)
        {
            cout<<"\t"<<y[i][j];
        }
        cout<<"\n\n";
    }

    for (int i = 0; i < nrow; i++)
    {
        delete []x[i];
        delete []c[i];
    }
    delete []x;delete[]c;
    return(y);
}
/******************************/
double* stationarydist(double** matrix,int nrow){
    double* res; create1D(nrow,0.0,res);
    /*
    double sumRow=0.0;
    for(int i=0; i<nrow; i++){
        sumRow +=matrix[0][i];
    }
    cout<<"sumRow===="<<sumRow<<endl;
    */
    for(int i=0; i<nrow; i++){
        res[i]=(double)matrix[0][i];
    }
    /*
    for(int j=0;j<nrow;j++)
        {
            cout<<"\t"<<res[j]<<"    ";
        }
        */
    return(res);

}
/**********************************/


