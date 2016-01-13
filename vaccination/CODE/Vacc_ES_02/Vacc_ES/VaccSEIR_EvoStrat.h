#ifndef VACCSEIR_EVOSTRAT_H
#define VACCSEIR_EVOSTRAT_H
// Situer les valeurs des situations
const int iS =1, iE = 2, iI = 3, iR=4, iN=5;
#define DEUX_PI ( 2.0 * 3.141592653589793238462643383279502884197169399375 ) // PI x 2
#define PI 3.141592653589
// Initialer une node pour chaque etat et chaque action
// Etat : est un état
// nbFoisVisite : est le nombre de fois de l'etat visité
// action : est une action choisie selon l'etat
// r : est une valeur de qualite de l'etat et l'action
struct NodeEtat{
    vector<int> Etat;
    int nbFoisVisite;
    vector<int> action;
    double Q;
};

int mymain(int,char**,double_vector*);

// Fonction d'initialisation la table qui contient les nombres de personnes dans chaque groupe
void InitialerY(int nbVill,int neq, unsigned long **&y, unsigned long rnitN, unsigned long rnitI);

void initialerVarSEIR(int nbVill,int nvar,unsigned long **&y,
                                 unsigned long S0, unsigned long E0, unsigned long I0, unsigned long R0);

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
#endif // VACCSEIR_EVOSTRAT_H
