/*
	AUTHORS : Walid HAMMACHE, Yanis KHORSI, Rayhana ZIARA
	Students in M1 IHPS 2014 / 2015
	Project regarding incomplet factorisation (ILU)
	LANGUAGE USED : C
*/
#include <projet.h>

int main(int argc, char **argv)
{
	/* *** DECLARATION DES VARIABLES *** */
	int i;
	int **level;
	double **A, *b, **LUi, *x;
	double runtime, timer_start, timer_end;
	struct info_t info;
	
	MPI_initialize(argc, argv, &info);
	
	/* *** INITIALISATION DES VARIABLES *** */
	level = (int**) malloc(n * sizeof(int*));
	A = (double**) malloc(n * sizeof(double*));
	LUi = (double**) malloc(n * sizeof(double*));
	b = (double*) malloc(info.nloc * sizeof(double));
	x = (double*) calloc(info.nloc, sizeof(double));

	#pragma omp parallel for schedule(static)
	for(i = 0; i < n; i++)
	{
		level[i] = (int*) malloc(info.nloc * sizeof(int));
		A[i] = (double*) malloc(info.nloc * sizeof(double));
		LUi[i] = (double*) malloc(info.nloc * sizeof(double));
	}

	/* *** CORPS DU PROGRAMME PRINCIPAL *** */	
	timer_start = get_timer();
	
	//printf("Matrice A:\n");
	poisson2D(A, &info);
	//affichageMat(A, &info);
	//affichageMatSpy(A, &info);

	//printf("Vecteur b:\n");
	vecteur_b(b, &info);
	//affichageVect(b, &info);

	//prodMatVect(A, b, x, &info);
	//printf("vecteur x:\n");
	//affichageVect(x, &info);

	//printf("Matrice LUi\n");
	ilup(A, level, LUi, &info);
	//affichageMat(LUi, &info);
	affichageMatSpy(LUi,&info);

	timer_end = get_timer();

	runtime = diff_time(timer_end, timer_start);
	print_time(&info, runtime);

	/* *** LIBERATION DES RESSOURCES *** */
	free(A);
	free(b);
	free(level);
	free(LUi);
	free(x);
	
	MPI_Finalize();
	
	return 0;
}
