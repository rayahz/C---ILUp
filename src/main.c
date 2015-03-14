/*
	AUTHORS : Walid HAMMACHE, Yanis KHORSI, Rayhana ZIARA
	Students in M1 IHPS 2014 / 2015
	Project regarding incomplet factorisation (ILU)
	LANGUAGE USED : C
*/
#include <projet.h>

int main()
{
	/* *** DECLARATION DES VARIABLES *** */
	int i, **level;
	double **A, **LUc, *x, *residu, **LUi, *b;

	/* *** INITIALISATION DES VARIABLES *** */
	level = (int**) malloc(n * sizeof(int*));
	
	A = (double**) malloc(n * sizeof(double*));
	LUi = (double**) malloc(n * sizeof(double*));
	LUc = (double**) malloc(n * sizeof(double*));
	
	b = (double*) malloc(n * sizeof(double));
	residu = (double*) calloc(n, sizeof(double));
	x = (double*) calloc(n, sizeof(double));

	for(i = 0; i < n; i++)
	{
		level[i] = (int*) malloc(n *sizeof(int));
		A[i] = (double*) malloc(n * sizeof(double));
		LUi[i] = (double*) malloc(n * sizeof(double));
		LUc[i] = (double*) malloc(n * sizeof(double));
	}
	
	/* *** CORPS DU PROGRAMME PRINCIPAL *** */
	printf("Matrice A\n");
	poisson2D(A);
	//matrixMarket(A, "e05r0200.mtx");
	//affichageMat(A);
	affichageMatSpy(A);

	//printf("Matrice LUc\n");
	//LUfact(A, LUc);
	//affichageMat(LUc);
	//affichageMatSpy(A);

	printf("Matrice LUi\n");
	ilup(A, level, LUi);
	//affichageMat(LUi);
	affichageMatSpy(LUi);

	printf("Vecteur b\n");
	vecteur_b(b);
	affichageVect(b);
	
	printf("Ax = b (PCG)\n");
	PCG(A, x, b, LUi, residu);
	affichageVect(x);
	
	printf("Vecteur residu issu du PCG\n");
	//affichageVect(residu);
	
	/* *** LIBERATION DES RESSOURCES *** */
	free(A);
	free(level);
	free(LUc);
	free(b);
	free(x);
	free(residu);
	free(LUi);
	
	return 0;
}
