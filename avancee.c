/* ***** LISTE DES FONCTIONS SEQUENTIELLES ***** */
void affichageMat(double**); 
void affichageMatSpy(double**); 
void affichageVect(double*); 
int PCG(double**, double*, double*, double**, struct info_t *); 
int CGR(double **, double *, double *, struct info_t *);
//void LUfact(double**, double**);
//int CG(double **, double *, double *, struct info_t *);
void matrixMarket(double **, char*);

/* ***** LISTE DES FONCTIONS PARALLELES ***** */
void poisson2D(double**);
void MPI_initialize(int, char **, struct info_t *);
double prodScal(double *, double *);
void prodVect(double **, double *, double *);
void vecteur_b(double*); 
void ilup(double**, int**, double**); 

/* ***** LISTE DES FONCTIONS PARALLELES VIA ***** */
double norme(double*); 
double produitT(double*, double**); 

/* ***** LISTE DES FONCTIONS QUI ONT ETE CORRIGEES ***** */
int PCG(double**, double*, double*, double**, struct info_t *); 
void LUfact(double**, double**);
int CGR(double **, double *, double *, struct info_t *);
