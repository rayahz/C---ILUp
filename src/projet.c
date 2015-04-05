/*
	AUTHORS : Walid HAMMACHE, Yanis KHORSI, Rayhana ZIARA
	Students in M1 IHPS 2014 / 2015
	Project regarding incomplet factorisation (ILU)
	LANGUAGE USED : C
*/
#include <projet.h>

/*
  PROCEDURE : ecrireFichier
  DESCRIPTION : permet d'ecrire dans un fichier de sortie les données du programme
  IN : structure info_t, nb iteration PCG, nb iteration CGR et vecteur resultat
  OUT : /
*/
void ecrireFichier(struct info_t *info, int x, int y, double *V)
{
	FILE* fichier = NULL;
	char buffer1[256], buffer2[256];
	int nbthread = omp_get_num_threads(), i;
	double A[n];
	
	time_t timestamp = time(NULL);
	strftime(buffer1, sizeof(buffer1), "%A %d %B %Y - %X", localtime(&timestamp));
	strftime(buffer2, sizeof(buffer2), "%d%m%Y_%H%M%S.txt", localtime(&timestamp));
	
	fichier = fopen(buffer2, "w");

	if(fichier != NULL)
	{
		fprintf(fichier, "#######################################\n");
		fprintf(fichier, "### %s ###\n", buffer1);
		fprintf(fichier, "#######################################\n");
		fprintf(fichier, "nx = %d\n", nx);
		fprintf(fichier, "ny = %d\n", ny);
		fprintf(fichier, "n = nx * ny = %d \n", n);
		fprintf(fichier, "p = %d\n\n", p);
		
		#pragma omp parallel
		{
			nbthread++;
		}

		fprintf(fichier, "Nombre de processus MPI : %d\n", info->nproc);
		fprintf(fichier, "Nombre de thread : %d\n", nbthread);
		fprintf(fichier, "Temps d'execution totale : %e\n\n", info->temps);
		
		fprintf(fichier, "Nombre d'iteration PCG : %d\n", x);
		fprintf(fichier, "Nombre d'iteration CGR : %d\n\n", y);
		
		fprintf(fichier, "Vecteur resultat x\n");

		MPI_Gather(V, info->nloc, MPI_DOUBLE, A + info->rang * info->nloc, info->nloc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if(info->rang == 0)
		{
			for(i = 0; i < n; i++)
				fprintf(fichier, "%.2lf\t", A[i]);
			fprintf(fichier, "\n\n");
		}
		fclose(fichier);
	}
}

/*
  PROCEDURE : affichageMat
  DESCRIPTION : affichage de la matrice sur la sortie standard
  IN : matrice et structure info_t
  OUT : /
*/
void affichageMat(double **M, struct info_t *info)
{
	int i, j, k;

	for(k = 0; k < info->nproc; k++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if(info->rang == k)
		{
			for(i = 0; i < info->nloc; i++)
			{
				for(j = 0; j < n; j++)
					fprintf(stdout, "%.2lf   ", M[i][j]);
				fprintf(stdout, "\n");
			}
			
			if(info->rang == info->nproc - 1)
				fprintf(stdout, "\n");
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

/*
  PROCEDURE : affichageMatSpy
  DESCRIPTION : affichage d'une matrice sous forme spy
  IN : matrice et structure info_t
  OUT : /
*/
void affichageMatSpy(double **M, struct info_t *info)
{
	int i, j, k;
	
	for(k = 0; k < info->nproc; k++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if(info->rang == k)
		{
			for(i = 0; i < info->nloc; i++)
			{
				for(j = 0; j < n; j++)
					if(fabs(M[i][j]) > DBL_EPSILON)
						fprintf(stdout, "* ");
					else fprintf(stdout, "  ");
				fprintf(stdout, "\n");
			}
			
			if(info->rang == info->nproc - 1)
				fprintf(stdout, "\n");
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
}

/*
  PROCEDURE : affichageVect
  DESCRIPTION : affichage d'un vecteur sur la sortie standard
  IN : vecteur et structure info_t
  OUT : /
*/
void affichageVect(double *V, struct info_t *info)
{
	int i;
	double A[n];

	MPI_Gather(V, info->nloc, MPI_DOUBLE, A + info->rang * info->nloc, info->nloc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(info->rang == 0)
	{
		for(i = 0; i < n; i++)
			fprintf(stdout, "%.2lf   ", A[i]);
		fprintf(stdout, "\n\n");
	}
}

/*
  PROCEDURE : poisson2D
  DESCRIPTION : permet de remplir une matrice à l'aide de l'equation de poisson 2D
  IN : matrice à remplir et structure info_t
  OUT : matrice 2D A
*/
void poisson2D(double **A, struct info_t *info)
{
int i, v1 = 1, v2 = 0, v3 = nx, v4 = 0;

	if(info->rang == 0)
	{
		v1 = 0;
		v3 = 0;
	}

	if(info->rang == info->nproc - 1)
	{
		v2 = 1;
		v4 = nx;
	}

	// superdiagonale inférieure
	#pragma omp parallel for schedule(static)
	for(i = nx - v3; i < info->nloc; i++)
		A[i][i + info->ideb - nx] = -1.0;

	// diagonale inférieure
	#pragma omp parallel for schedule(static)
	for(i = 1 - v1; i < info->nloc; i++)
		A[i][i + info->ideb - 1] = -1.0;

	// diagonale principale
	#pragma omp parallel for schedule(static)
	for(i = 0; i < info->nloc; i++)
		A[i][i + info->ideb] = 4.0;

	// diagonale supérieure
	#pragma omp parallel for schedule(static)
	for(i = 0; i < info->nloc - v2; i++)
		A[i][i + info->ideb + 1] = -1.0;

	// superdiagonale supérieure
	#pragma omp parallel for schedule(static)
	for(i = 0; i < info->nloc - v4; i++)
		A[i][i + info->ideb + nx] = -1.0;
}

/*
  PROCEDURE : ilup
  DESCRIPTION : permet d'effectuer la factorisation LU incomplete au rang p
  IN : matrice initiale, matrice LUi et structure info_t
  OUT : matrice LUi sous forme LU incomplete 
*/
void ilup(double **A, double **LUi, struct info_t *info)
{
	int i, j, k;
	int **level;
	
	level = (int**) malloc(info->nloc * sizeof(int*));
	for(i = 0; i < n; i++)
		level[i] = (int*) malloc(n * sizeof(int));
	
	// initialisation de la matrice de remplissage 
	// si A(i,j) != 0
		// alors lev(i,j) = 0
		// sinon lev(i,j) = inf
	for(i = 0; i < info->nloc; i++)
	{
		for(j = 0; j < n; j++)
		{
			if((fabs(A[i][j]) > DBL_EPSILON) || (i == (j+info->rang*info->nloc)))
				level[i][j] = 0;
			else level[i][j] = inf;

			// copie de la matrice A dans la matrice LUi
			LUi[i][j] = A[i][j]; 
		}
	}
	
	int deb = info->ideb;
	if(info->rang==0)
		deb++;
		
	for(i = deb; i < info->ifin; i++)
	{
		for(k = 0; k <= i-1; k++)
		{
			if(level[i][k-info->nloc*info->rang] <= p)
			{
				LUi[i][k-info->nloc*info->rang] /= LUi[k][k-info->nloc*info->rang]; // calcul
			
				for(j = k +1; j < info->nloc; j++)	
				{
					LUi[i][j] -= LUi[i][k] * LUi[k][j]; // calcul
					
					// remplissage de la matrice de niveau de remplissage
						// factorisation symbolique : Le but de cette étape est de calculer la structure de remplissage des facteurs avant l’étape de factorisation numérique.
					if(fabs(LUi[i][j]) > DBL_EPSILON)
						level[i][j] = min(level[i][j], level[i][k] + level[k][j] + 1);
				}
			}
		}
		
		// tous les elements dont le niveau de remplissage est superieur à p sont remplaces par 0
		for(j = 0 ; j < info->nloc; j++)
		{
			if(level[i][j] > p)
				LUi[i][j] = 0.0;
		}
	}
	
	/*for(k = 0; k < info->nloc; k++)
	{
		for(i = k+1; i < n; i++)
		{
			if(level[i][k] <= p)
				LUi[i][k] = LUi[i][k] / LUi[k+info->nloc*info->rang][k];
		}

		for(j = k + 1; j < info->nloc; j++)
		{
			for(i = k +1 ; i < n; i++)
			{
				if(level[i][k] <= p)
				{
					LUi[i][j] -= LUi[i][k] * LUi[k][j]; 

					if(fabs(LUi[i][j]) > DBL_EPSILON)
						level[i][j] = min(level[i][j], level[i][k] + level[k][j] + 1);
				}
			}
		}
	}
	
	for(j = 0 ; j < info->nloc; j++)
	{
		for(i = 0 ; i < n; i++) 
		{
			if(level[i][j] > p)
				LUi[i][j] = 0.0;
		}
	}*/
	
	free(level);
}

/*
  PROCEDURE : LUfact
  DESCRIPTION : permet d'effectuer la factorisation LU complete 
  IN : matrice à factoriser et matrice factorisée 
  OUT : matrice LU --> A = LU
*/
void LUfact(double **A, double **LU)
{
	int i, j, k;
	
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
			LU[i][j] = A[i][j];
	}
	
	for(k = 0; k < n - 1; k++)
	{
		for(i = k + 1; i < n; i++)
		{
			LU[i][k] = LU[i][k] / LU[k][k];
			
			for(j = k + 1; j < n; j++)
				LU[i][j] = LU[i][j] - (LU[i][k] * LU[k][j]);
		}
	}
}

/*
  FONCTION : norme
  DESCRIPTION : permet de calculer la norme d'un vecteur
  IN : vecteur à calculer et structure info_t
  OUT : /
  RETOUR : scalaire contenant la norme du vecteur
*/
double norme(double *a, struct info_t *info)
{
	return sqrt(fabs(prodScal(a, a, info)));
}

/*
  FONCTION : produitT
  DESCRIPTION : permet d'effectuer le produit d'un vecteur transpose avec une matrice puis avec le vecteur (v' * a * v)
  IN : vecteur, matrice et structure info_t
  OUT : /
  RETOUR : scalaire contenant le resultat du produit
*/
double produitT(double *v, double **a, struct info_t *info)
{
	double *r1, r2 = 0.0;
	
	r1 = (double*) calloc(n, sizeof(double));
	prodMatVect(a, v, r1, info);
	r2 = prodScal(r1, v, info);
	free(r1);
	
	return r2;
}

/*
  FONCTION : PCG
  DESCRIPTION : permet d'effectuer le gradient conjugue preconditionne à partir d'une matrice conditionnée (LU incomplete)
  IN : matrice initiale, vecteur recherche, vecteur resultat, matrice A sous forme ILU et structure info_t
  OUT : vecteur resultat x
  RETOUR : scalaire contenant le nombre d'iteration
*/
int PCG(double **A, double *x, double *b, double **B, struct info_t *info)
{
	int iterG = 0, iterL = 0, i, j, maxiter = max(100.0, sqrt(n));
	double alpha, beta, M, Mold, bnrm2 = norme(b, info);
	double *Ap, *Br, *r, *v;

	Ap = (double*) malloc(info->nloc * sizeof(double));
	Br = (double*) malloc(info->nloc * sizeof(double));
	r = (double*) malloc(info->nloc * sizeof(double));
	v = (double*) calloc(info->nloc, sizeof(double));

	#pragma omp parallel for schedule(static)
	for(i = 0; i < info->nloc; i++)
		r[i] = b[i];

	prodMatVect(B, r, v, info);
	M = prodScal(r, v, info);

	if(bnrm2 == 0.0) 
		bnrm2 = 1.0;

	while(((norme(r, info) / bnrm2) > DBL_EPSILON) && (iterG < maxiter))
	{
		// A VIRER QUAND CA MARCHE, DE ICI
		iterL++;
		if(info->rang == 0)
			iterG = iterL;
		MPI_Bcast(&iterG, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if(info->rang != 0 && iterG != iterL)
			MPI_Abort(MPI_COMM_WORLD, 1);
		// A DE LA
		
		#pragma omp parallel for schedule(static)
		for(i = 0; i < info->nloc; i++)
			Ap[i] = 0.0;
			
		prodMatVect(A, v, Ap, info);
		alpha = M / prodScal(v, Ap, info);

		#pragma omp parallel for schedule(static)
		for(i = 0; i < info->nloc; i++)
		{
			x[i] += alpha * v[i];
			r[i] -= alpha * Ap[i];
		}

		Mold = M;

		#pragma omp parallel for schedule(static)
		for(i = 0; i < info->nloc; i++)
			Br[i] = 0.0;
			
		prodMatVect(B, r, Br, info);
		M = prodScal(r, Br, info);
		beta = M / Mold;

		#pragma omp parallel for schedule(static)
		for(i = 0; i < info->nloc; i++)
			v[i] = Br[i] + beta * v[i];
	}
	
	free(Ap);
	free(Br);
	free(r);
	free(v);
	return iterG;
}

/*
  PROCEDURE : vecteur_b
  DESCRIPTION : permet de construire le vecteur b issu du systeme lineaire Ax = b
  IN : vecteur à remplir et structure info_t
  OUT : vecteur b
*/
void vecteur_b(double *b, struct info_t *info)
{
	int i, j;

	#pragma omp parallel for schedule(static)
	for(i = 0; (i < info->nloc) && ((i + info->ideb) < n / 2); i++)
		b[i] = (double)(info->ideb + i);

	if(n / 2 < info->ifin)
	{
		if(n / 2 > info->ideb)
		{
			if(n % 2 == 0)
				j = n / 2 - 1;
			else j = n / 2;
			
			#pragma omp parallel for schedule(static)
			for(i = info->nloc / 2; i < info->nloc; i++, j--)
				b[i] = (double)j;
		}else{
			#pragma omp parallel for schedule(static)
			for(i = 0, j = n - info->ideb - 1; i < info->nloc; i++, j--)
				b[i] = (double)j;
		}
	}
}

/*
  PROCEDURE : matrixMarket
  DESCRIPTION : permet de construire une matrice à partir d'un fichier de MatrixMarket
  IN : matrice à remplir et nom du fichier
  OUT : matrice A
*/
void matrixMarket(double **A, char *nom)
{
	FILE *fichier;
	char mot[100], chemin[50] = "./Matrix/"; 
	int cpt = 0, end = 10, cpt2 = 0, i, j; 
	
	strcat(chemin, nom);
	fichier = fopen(chemin, "r");

	while(cpt != end) 
	{
		fscanf(fichier, "%s", mot);
		
		// permet d'obtenir le nombre d'elements non nuls se trouvant dans la matrice (ligne 2 mot 3)
		if(cpt == 7)
		{
			sscanf(mot, "%d", &end);
			end = end * 3 + 8;
		}
		
		// à partir de la ligne 3 incluse
		if(cpt > 7)
		{
			if(cpt2 == 0) // ligne "cpt + 1" mot 1
				sscanf(mot, "%d", &i);

			if(cpt2 == 1) // ligne "cpt + 1" mot 2
				sscanf(mot, "%d", &j);
			
			if(cpt2 == 2) // ligne "cpt + 1" mot 3
			{
				sscanf(mot, "%lf", &A[i-1][j-1]);
				sscanf(mot, "%lf", &A[j-1][i-1]);
			}	
			
			if(cpt2 < 2)
				cpt2++;
			else cpt2 = 0;
		}
		
		cpt++;
	}
	
	fclose(fichier);
}

/*
  FONCTION : CGR
  DESCRIPTION : permet d'effectuer le gradient conjugue residuel à partir d'une matrice preconditionne
  IN : matrice initiale, vecteur resultat, vecteur recherche et structure info_t
  OUT : vecteur resultat x
  RETOUR : scalaire contenant le nombre d'iteration
*/
int CGR(double **A, double *x, double *b, struct info_t *info)
{
	int iter = 0, i, j, k, maxiter = max(100.0, sqrt(n));
	double alpha, bnrm2 = norme(b, info);
	double *Ap, *Ar, *beta, *r, **v;
	
	Ap = (double*) malloc(info->nloc * sizeof(double));
	Ar = (double*) malloc(info->nloc * sizeof(double));
	r = (double*) calloc(info->nloc, sizeof(double));
	v = (double**) calloc(maxiter + 1, sizeof(double*));
	beta = (double*) malloc(maxiter * sizeof(double));
	
	for(i = 0; i <= maxiter; i++)
		v[i] = (double*) calloc(info->nloc, sizeof(double));

	prodMatVect(A, x, r, info);
	
	#pragma omp parallel for schedule(static)
	for(i = 0; i < info->nloc; i++)
	{
		r[i] = b[i] - r[i];
		v[0][i] = r[i];
	}
	
	if(bnrm2 == 0.0) 
		bnrm2 = 1.0;

	while(((norme(r, info) / bnrm2) > DBL_EPSILON) && (iter < maxiter))
	{
		iter++;

		#pragma omp parallel for schedule(static)
		for(i = 0; i < info->nloc; i++)
			Ap[i] = 0.0;
			
		prodMatVect(A, v[iter-1], Ap, info);
		alpha = prodScal(r, Ap, info) / prodScal(Ap, Ap, info);
		
		#pragma omp parallel for schedule(static)
		for(i = 0; i < info->nloc; i++)
		{
			x[i] += alpha * v[iter - 1][i];
			r[i] -= alpha * Ap[i];
		}
		
		#pragma omp parallel for schedule(static)
		for(k = 0; k < iter; k++)
		{
			for(i = 0; i < info->nloc; i++)
			{
				Ap[i] = 0.0;
				Ar[i] = 0.0;
			}
			
			prodMatVect(A, v[k], Ap, info);
			prodMatVect(A, r, Ar, info);
			beta[k] = - prodScal(Ar, Ap, info) / prodScal(Ap, Ap, info);
		}
		
		#pragma omp parallel for schedule(static)
		for(k = 0; k < iter; k++)
		{
			for(i = 0; i < n; i++)
				v[iter][i] += beta[k] * v[k][i];
		}
		
		#pragma omp parallel for schedule(static)
		for(i = 0; i < info->nloc; i++)
			v[iter][i] += r[i];
	}
	
	free(Ap);
	free(Ar);
	free(r);
	free(v);
	free(beta);
	
	return iter;
}

/*
  PROCEDURE : prodMatVect
  DESCRIPTION : permet de calculer le produit d'une matrice avec un vecteur
  IN : matrice, vecteur de calcul, vecteur resultat et structure info_t
  OUT : vecteur resultat r
*/
void prodMatVect(double **m, double *v, double *r, struct info_t *info)
{
	int i, j;
	double tmp[n];

	MPI_Allgather(v, info->nloc, MPI_DOUBLE, tmp, info->nloc, MPI_DOUBLE, MPI_COMM_WORLD);

	#pragma omp parallel for schedule(static) private(j, r)
	for(i = 0; i < info->nloc; i++)
		for(j = 0; j < n; j++)
			r[i] += m[i][j] * tmp[j];
}

/*
  FONCTION : prodScal
  DESCRIPTION : permet de calculer le produit scalaire d'un vecteur
  IN : vecteur 1, vecteur 2 et structure info_t
  OUT : /
  RETOUR : scalaire contenant le resultat du produit
*/
double prodScal(double *x, double *y, struct info_t *info)
{
	int i;
	double res_loc = 0.0, res_glob;

	#pragma omp parallel for schedule(static) reduction(+:res_loc)
	for(i = 0; i < info->nloc; i++)
		res_loc += x[i] * y[i];

	MPI_Reduce(&res_loc, &res_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	return res_glob;
}

/*
  FONCTION : CG
  DESCRIPTION : permet d'effectuer le gradient conjugue à partir d'une matrice 
  IN : matrice initiale, vecteur resultat, vecteur recherche, vecteur residu et structure info_t
  OUT : vecteur resultat x
  RETOUR : scalaire contenant le nombre d'iteration
*/
int CG(double **A, double *b, double *x, struct info_t *info)
{
	int i, j, iter = 0, maxiter = max(100.0, sqrt(n));
	double alpha, M, Mold, beta;
	double *v, *r;

	v = (double*) malloc(info->nloc * sizeof(double));
	r = (double*) malloc(info->nloc * sizeof(double));

	for(i = 0; i < info->nloc; i++)
	{
		r[i] = b[i];
		v[i] = r[i];
	}

	M = prodScal(r, r, info);

	while(((norme(r, info) / norme(b, info)) > DBL_EPSILON) && (iter < maxiter))
	{
		iter++;
		alpha = M / produitT(v, A, info);

		for(i = 0; i < info->nloc; i++)
			x[i] += alpha * v[i];

		prodMatVect(A, v, r, info);
		for(i = 0; i < info->nloc; i++)
			r[i] -= alpha * r[i];
		
		Mold = M;
		M = prodScal(r, r, info);
		beta = M / Mold;

		for(i = 0; i < info->nloc; i++)
			v[i] = r[i] + beta * v[i];
	}
	
	free(v);
	free(r);
	return iter;
}

/*
  PROCEDURE : MPI_initialize
  DESCRIPTION : permet d'initialiser l'environnement MPI
  IN : nombre d'arguments, vecteur d'arguments et structure info_t
  OUT : /
*/
void MPI_initialize(int argc, char **argv, struct info_t *info)
{
	int Q, R;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &(info->rang));
	MPI_Comm_size(MPI_COMM_WORLD, &(info->nproc));

	info->temps = 0.0;
	Q = n / info->nproc;
	R = n % info->nproc;

	if(info->rang < R)
	{
		info->nloc = Q + 1;
		info->ideb = info->rang * (Q + 1);
		info->ifin = info->ideb + info->nloc;
	}else{
		info->nloc = Q;
		info->ideb = R * (Q+1) + (info->rang - R) * Q;
		info->ifin = info->ideb + info->nloc;
	}
}

/*
  FONCTION : get_timer
  DESCRIPTION : permet d'initialiser le timer
  IN : /
  OUT : /
  RETOUR : retourne le temps actuel
*/
double get_timer()
{
	struct timeval t;

	gettimeofday(&t, NULL);

	return t.tv_sec * 1000000 + t.tv_usec;
}

/*
  FONCTION : diff_time
  DESCRIPTION : permet d'effectuer la difference entre le temps de fin et le temps de debut
  IN : temps de fin et temps de debut
  OUT : /
  RETOUR : retourne la difference de temps
*/
double diff_time(double end, double start)
{
	return end - start;
}

/*
  PROCEDURE : print_time
  DESCRIPTION : permet d'afficher le temps d'execution de chaque processus
  IN : structure info_t et temps d'execution
  OUT : /
*/
void print_time(struct info_t *info, double runtime)
{
	int i;
	double runtime_in_seconds = runtime / 1000000;
	double runtime_t[info->nproc];

	MPI_Gather(&runtime_in_seconds, 1, MPI_DOUBLE, &runtime_t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(info->rang == 0)
	{
		for(i = 0; i < info->nproc; i++)
		{
			if(info->temps < runtime_t[i])
				info->temps = runtime_t[i];
				
			fprintf(stdout, "[%d] Temps de l'execution : %e s\n", i, runtime_t[i]);
		}
	}
}
