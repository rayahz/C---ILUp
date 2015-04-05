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
	int i, iterPCG = 0, iterCGR = 0;
	double **A, *b, **LUi, *x;
	double runtime, timer_start, timer_end;
	struct info_t info;
	
	MPI_initialize(argc, argv, &info);
	
	/* *** INITIALISATION DES VARIABLES *** */
	A = (double**) malloc(info.nloc * sizeof(double*));
	LUi = (double**) malloc(info.nloc * sizeof(double*));
	b = (double*) malloc(info.nloc * sizeof(double));
	x = (double*) calloc(info.nloc, sizeof(double));

	#pragma omp parallel for schedule(static)
	for(i = 0; i < info.nloc; i++)
	{
		A[i] = (double*) malloc(n * sizeof(double));
		LUi[i] = (double*) malloc(n * sizeof(double));
	}

	/* *** CORPS DU PROGRAMME PRINCIPAL *** */
	timer_start = get_timer();
	
	poisson2D(A, &info);
	//affichageMat(A, &info);
	//affichageMatSpy(A, &info);
	
	vecteur_b(b, &info);
	//affichageVect(b, &info);
	
	//ilup(A, LUi, &info);
	//affichageMat(LUi, &info);
	//affichageMatSpy(LUi,&info);
	
	//iterPCG = PCG(A, x, b, LUi, &info);
	//iterCGR = CGR(LUi, x, b, &info);
	
	prodMatVect(A, b, x, &info);
	//affichageVect(x, &info);
	
	timer_end = get_timer();
	runtime = diff_time(timer_end, timer_start);
	print_time(&info, runtime);
	
	ecrireFichier(&info, iterPCG, iterCGR, x);

	/* *** LIBERATION DES RESSOURCES *** */
	free(A);
	free(b);
	free(LUi);
	free(x);
	
	MPI_Finalize();
	
	return 0;
}
