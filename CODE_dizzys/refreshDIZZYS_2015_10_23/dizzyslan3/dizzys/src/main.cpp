#include <limits>
#include "seir_stoch.h"
#include "rangenyan.h"
using namespace std;
extern "C" {
	int main(int argc, char *argv[])
	{

    		cout<<"\nBonjour \nCommencer  simuler le modele SEIR !"<<endl;

	    //Creer un repertoire pour sauvegarder tous les fichiers de resultats
	    char comRepertoire[255] = "";
	    char nomRepertoire[255] = "Output";
	    sprintf(comRepertoire, "%s%s", "mkdir ", nomRepertoire);
	    if(system(comRepertoire)) cout <<"Repertoire contenant les fichiers des resultats Output est en train de exister!"<<endl;
	    else
	    cout<<"Creation du repertoire Output est avec succes!"<<endl;

	    vector<double> a;
	    a.push_back(1);
	   // a.push_back(34);
	   // a.push_back(23);
	    a = resizeVector1D(a,5);

	     for(int i=0; i<a.size();i++) cout<<a[i]<<endl;

	    vector<double> valParSIM = initialervalParSIM_00(argc,argv);

	    vector<vector<double> > valSEIR0;
	    vector<double> ville1;
	    ville1.push_back(1);
	    ville1.push_back(20000);
	    ville1.push_back(0);
	    ville1.push_back(100);
	    ville1.push_back(10000);
	    ville1.push_back(30100);
	    valSEIR0.push_back(ville1);
	   vector<vector<vector<double> > >tabVille =  resSimSEIR(valParSIM,valSEIR0);
	   writeInFile_00(valParSIM,tabVille);
	    // 19/03/2013
	   double sigma = std::numeric_limits<double>::infinity();
	   if(sigma==INFINITY) cout<<"gianggiang"<<endl;
	   else
	       cout<<"khuongkhuong"<<endl;
	    return (0);
	 }
    }
}
