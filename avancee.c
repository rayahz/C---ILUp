/* ***** LISTE DES FONCTIONS SEQUENTIELLES ***** */
int PCG(double**, double*, double*, double**, struct info_t *); 
int CGR(double **, double *, double *, struct info_t *);
//void LUfact(double**, double**);
//int CG(double **, double *, double *, struct info_t *);
//double produitT(double*, double**);
void matrixMarket(double **, char*);
double get_timer();
double diff_time(double, double);

/* ***** LISTE DES FONCTIONS PARALLELES ***** */
void poisson2D(double**, struct info_t *);
void MPI_initialize(int, char **, struct info_t *);
double prodScal(double *, double *, int);
void prodMatVect(double **, double *, double *, struct info_t *);
void vecteur_b(double*, struct info_t *); 
void ilup(double**, int**, double**, struct info_t *); 
void affichageMat(double**, struct info_t *); 
void affichageMatSpy(double**, struct info_t *); 
void affichageVect(double*, struct info_t *); 
void print_time(struct info_t *, double );

/* ***** LISTE DES FONCTIONS PARALLELES VIA ***** */
double norme(double*); 

/* ***** LISTE DES FONCTIONS QUI ONT ETE CORRIGEES ***** */
int PCG(double**, double*, double*, double**, struct info_t *); 
void LUfact(double**, double**);
int CGR(double **, double *, double *, struct info_t *);

