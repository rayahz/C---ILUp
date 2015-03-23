/*
	AUTHORS : Walid HAMMACHE, Yanis KHORSI, Rayhana ZIARA
	Students in M1 IHPS 2014 / 2015
	Project regarding incomplet factorisation (ILU)
	LANGUAGE USED : C
*/
#include <projet.h>

int main(int argc, char **argv)
{
	struct timespec t1, t2;
	clock_gettime(CLOCK_REALTIME, &t1);
	
	/* *** DECLARATION DES VARIABLES *** */
	int i, **level;
	double **A, **LUc, *x, **LUi, *b;
	
	/* *** INITIALISATION DES VARIABLES *** */
	level = (int**) malloc(n * sizeof(int*));
	A = (double**) malloc(n * sizeof(double*));
	LUi = (double**) malloc(n * sizeof(double*));
	LUc = (double**) malloc(n * sizeof(double*));
	b = (double*) malloc(n * sizeof(double));
	x = (double*) calloc(n, sizeof(double));

	for(i = 0; i < n; i++)
	{
		level[i] = (int*) malloc(n *sizeof(int));
		A[i] = (double*) malloc(n * sizeof(double));
		LUi[i] = (double*) malloc(n * sizeof(double));
		LUc[i] = (double*) malloc(n * sizeof(double));
	}

	/* *** CORPS DU PROGRAMME PRINCIPAL *** */	
	struct info_t info;
	//MPI_initialize(argc, argv, &info);
	
	//printf("Matrice A\n");
	poisson2D(A);
	//matrixMarket(A, "e05r0200.mtx");
	//affichageMat(A);
	//affichageMatSpy(A);

	//printf("Matrice LUc\n");
	//LUfact(A, LUc);
	//affichageMat(LUc);
	//affichageMat(A);
	//affichageMatSpy(A);

	//printf("Matrice LUi\n");
	ilup(A, level, LUi);
	//affichageMat(LUi);
	//affichageMatSpy(LUi);

	//printf("Vecteur b\n");
	vecteur_b(b);
	//affichageVect(b);
	
	//printf("Ax = b (PCG)\n");
	PCG(A, x, b, LUi, &info);
	//affichageVect(x);

	//printf("Vecteur residu issu du PCG\n");
	//affichageVect(residu);
	
	/*printf("CGR\n");
	int iter = CGR(LUi, x, b);
	printf("x \n");
	affichageVect(x);
	printf("iter=%d\n",iter);*/
	
	/* *** LIBERATION DES RESSOURCES *** */
	free(A);
	free(level);
	free(LUc);
	free(b);
	free(x);
	free(LUi);
	
	//MPI_Finalize();
	
	clock_gettime(CLOCK_REALTIME, &t2);
	printf("%lg\n",1.*t2.tv_sec-t1.tv_sec+(t2.tv_nsec-t1.tv_nsec)/1e9);
	
	return 0;
}
