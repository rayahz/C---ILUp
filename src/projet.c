/*
	AUTHORS : Walid HAMMACHE, Yanis KHORSI, Rayhana ZIARA
	Students in M1 IHPS 2014 / 2015
	Project regarding incomplet factorisation (ILU)
	LANGUAGE USED : C
*/
#include <projet.h>

/*
  PROCEDURE : affichageMat
  DESCRIPTION : permet l'affichage de la matrice sur la sortie standard
  IN : matrice
  OUT : /
*/
void affichageMat(double **M, struct info_t *info)
{
	int i, j;
	double A[n][n];

	for(i = 0; i < n; i++)
		MPI_Gather(M[i], info->nloc, MPI_DOUBLE, A[i] + info->rang * info->nloc, info->nloc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(info->rang == 0)
	{
		for(i = 0; i < n; i++)
		{
			for(j = 0; j < n; j++)
				printf("%.2lf   ", A[i][j]);
			printf("\n");
		}
		printf("\n");
	}
}

/*
  PROCEDURE : affichageMatSpy
  DESCRIPTION : permet l'affichage d'une matrice sous forme spy
  IN : matrice 
  OUT : /
*/
void affichageMatSpy(double **M, struct info_t *info)
{
	int i, j;
	double A[n][n];

	for(i = 0; i < n; i++)
		MPI_Gather(M[i], info->nloc, MPI_DOUBLE, A[i] + info->rang * info->nloc, info->nloc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(info->rang == 0)
	{
		for(i = 0; i < n; i++)
		{
			for(j = 0; j < n; j++)
				if(fabs(A[i][j]) > DBL_EPSILON)
					fputs("* ", stdout);
				else fputs("  ", stdout);
			printf("\n");
		}
		printf("\n");
	}
}

/*
  PROCEDURE : affichageVect
  DESCRIPTION : permet l'affichage d'un vecteur sur la sortie standard
  IN : vecteur
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
			printf("%.2lf   ",A[i]);
		printf("\n\n");
	}
}

/*
  PROCEDURE : poisson2D
  DESCRIPTION : permet de remplir une matrice à l'aide de l'equation de poisson 2D
  IN : matrice à remplir
  OUT : matrice 2D A
*/
void poisson2D(double **A, struct info_t *info)
{
	int i, v1 = 1, v2 = 0, v3 = 0;

	if(info->rang == 0)
		v1 = 0;
		
	if(info->rang == info->nproc - 1)
		v2 = 1;
		
	if((info->ifin - 1 - nx) >= 0)
		v3 = nx;

	#pragma omp parallel for schedule(static)
	for(i = info->ideb; i < info->ifin; i++){
		A[i][i - info->nloc * info->rang] = 4.0;
		A[i - v1][i-info->nloc * info->rang + 1 - v1] = -1.0;
	}

	#pragma omp parallel for schedule(static)
	for(i = info->ideb; i < info->ifin - v2; i++)
		A[i + 1][i - info->nloc * info->rang] = -1.0;

	#pragma omp parallel for schedule(static)
	for(i = info->ifin-1; (i >= info->ideb) && ((i - v3) >= 0); i--)
		A[i - v3][i - info->nloc * info->rang + nx - v3] = -1.0;

	#pragma omp parallel for schedule(static)
	for(i = info->ideb; (i < info->ifin) && (i < (n - nx)); i++)
		A[i + nx][i - info->nloc * info->rang] = -1.0;
}

/*
  PROCEDURE : ilup
  DESCRIPTION : permet d'effectuer la factorisation LU incomplete au rang p
  IN : matrice initiale, matrice de niveaux de remplissage et matrice LUi 
  OUT : matrice LUi sous forme LU incomplet et matrice lev contenant le niveau de remplissage pour chaque element
*/
void ilup(double **A, int **level, double **LUi, struct info_t *info)
{
	int i, j, k;
	
	// initialisation de la matrice de remplissage 
	// si A(i,j) != 0
		// alors lev(i,j) = 0
		// sinon lev(i,j) = inf
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < info->nloc; j++)
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
}

/*
  PROCEDURE : LUfact
  DESCRIPTION : permet d'effectuer la factorisation LU complete 
  IN : matrice à factoriser et matrice factorisée 
  OUT : matrice LU --> A = LU
*/
/*void LUfact(double **A, double **LU)
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
}*/

/*
  FONCTION : norme
  DESCRIPTION : permet de calculer la norme d'un vecteur
  IN : vecteur à calculer
  OUT : /
  RETOUR : scalaire contenant la norme du vecteur
*/
double norme(double *a, int N)
{
	return sqrt(fabs(prodScal(a, a, N)));
}

/*
  FONCTION : produitT
  DESCRIPTION : permet d'effectuer le produit d'un vecteur transpose avec une matrice puis avec le vecteur (v' * a * v)
  IN : vecteur et matrice 
  OUT : /
  RETOUR : scalaire contenant le resultat du produit
*/
/*double produitT(double *v, double **a)
{
	double *r1, r2 = 0.0;
	
	r1 = (double*) calloc(n, sizeof(double));
	prodMatVect(a, v, r1);
	r2 = prodScal(r1, v);
	free(r1);
	
	return r2;
}*/

/*
  PROCEDURE : PCG
  DESCRIPTION : permet d'effectuer le gradient conjugue preconditionne à partir d'une matrice conditionnée (LU incomplete)
  IN : matrice initiale, vecteur recherche, vecteur resultat, matrice A sous forme ILU et vecteur residu
  OUT : vecteur resultat x
*/
int PCG(double **A, double *x, double *b, double **B, struct info_t *info)
{
	int iter = 0, i, j, maxiter = max(100.0, sqrt(n));
	double alpha, beta, M, Mold, bnrm2 = norme(b, info->nloc);
	double *Ap, *Br, *r, *v;

	Ap = (double*) malloc(n * sizeof(double));
	Br = (double*) malloc(n * sizeof(double));
	r = (double*) malloc(n * sizeof(double));
	v = (double*) calloc(n, sizeof(double));

	for(i = 0; i < n; i++)
		r[i] = b[i];

	prodMatVect(B, r, v, info);
	M = prodScal(r, v, info->nloc);

	if(bnrm2 == 0.0) 
		bnrm2 = 1.0;

	while(((norme(r, info->nloc) / bnrm2) > DBL_EPSILON) && (iter < maxiter))
	{
		iter++;

		for(i = 0; i < n; i++)
			Ap[i] = 0.0;
			
		prodMatVect(A, v, Ap, info);
		alpha = M / prodScal(v, Ap, info->nloc);

		for(i = 0; i < n; i++)
		{
			x[i] += alpha * v[i];
			r[i] -= alpha * Ap[i];
		}

		Mold = M;

		for(i = 0; i < n; i++)
			Br[i] = 0.0;
			
		prodMatVect(B, r, Br, info);
		M = prodScal(r, Br, info->nloc);
		beta = M / Mold;

		for(i = 0; i < n; i++)
			v[i] = Br[i] + beta * v[i];
	}
	
	free(Ap);
	free(Br);
	free(r);
	free(v);

	return iter;
}

/*
  PROCEDURE : vecteur_b
  DESCRIPTION : permet de construire le vecteur b issu du systeme lineaire Ax = b
  IN : vecteur à remplir
  OUT : vecteur b
*/
void vecteur_b(double *b, struct info_t *info)
{
	int i, j;

	#pragma omp parallel for schedule(static)
	for(i = 0, j = info->ideb; (i < info->nloc) && (j < n / 2); i++, j++)
		b[i] = (double)j;

	if(n / 2 < info->ifin)
	{
		if(n / 2 >= info->ideb)
		{
			#pragma omp parallel for schedule(static)
			for(i = n / 2 - info->ideb, j = n / 2 - 1; i < info->nloc; i++, j--)
				b[i] = (double)j;
		}else{
			#pragma omp parallel for schedule(static)
			for(i = 0, j = n - info->nloc * info->rang - 1; i < info->nloc; i++, j--)
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
  IN : matrice initiale, vecteur resultat, vecteur recherche et vecteur residu
  OUT : vecteur resultat x
  RETOUR : scalaire contenant le nombre d'iteration
*/
int CGR(double **A, double *x, double *b, struct info_t *info)
{
	int iter = 0, i, j, k, maxiter = max(100.0, sqrt(n));
	double alpha, bnrm2 = norme(b, info->nloc);
	double *Ap, *Ar, *beta, *r, **v;
	
	Ap = (double*) malloc(n * sizeof(double));
	Ar = (double*) malloc(n * sizeof(double));
	r = (double*) calloc(n, sizeof(double));
	v = (double**) calloc(n, sizeof(double*));
	
	for(i = 0; i < n; i++)
		v[i] = (double*) calloc(maxiter, sizeof(double));

	prodMatVect(A, x, r, info);
	
	for(i = 0; i < n; i++)
	{
		r[i] = b[i] - r[i];
		v[i][0] = r[i];
	}

	if(bnrm2 == 0.0) 
		bnrm2 = 1.0;

	while(((norme(r, info->nloc) / bnrm2) > DBL_EPSILON) && (iter < maxiter))
	{
		iter++;

		for(i = 0; i < n; i++)
		{
			Ap[i] = 0.0;
			for(j = 0; j < n; j++)
				Ap[i] += A[i][j] * v[j][iter - 1];
		}
		
		alpha = prodScal(r, Ap, info->nloc) / prodScal(Ap, Ap, info->nloc);
		
		for(i = 0; i < n; i++)
		{
			x[i] += alpha * v[i][iter - 1];
			r[i] -= alpha * Ap[i];
		}

		beta = (double*) malloc(iter * sizeof(double));
		
		for(k = 0; k < iter; k++)
		{
			for(i = 0; i < n; i++)
			{
				Ap[i] = 0.0;
				Ar[i] = 0.0;
				for(j = 0; j < n; j++)
					Ap[i] += A[i][j] * v[j][k];
				prodMatVect(A, r, Ar, info);
			}
			beta[k] = (-1) * prodScal(Ar, Ap, info->nloc) / prodScal(Ap, Ap, info->nloc);
		}
		
		for(k = 0; k < iter; k++)
		{
			for(i = 0; i < n; i++)
				v[i][iter] += beta[k] * v[i][k];
		}
		
		for(i = 0; i < n; i++)
			v[i][iter] += r[i];

		free(beta);
	}
	
	free(Ap);
	free(Ar);
	free(r);
	free(v);
	
	return iter;
}

void prodMatVect(double **m, double *v, double *r, struct info_t *info)
{
	int i;
	double tmp[n], tmp2[n];

	#pragma omp parallel for schedule(static) private(tmp)
	for(i = 0; i < n; i++)
		tmp[i] = prodScal(m[i], v, info->nloc);

	MPI_Allreduce(tmp, tmp2, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Scatter(tmp2, info->nloc, MPI_DOUBLE, r, info->nloc, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
}

double prodScal(double *x, double *y, int N)
{
	int i;
	double res_loc = 0.0, res_glob;

	#pragma omp parallel for schedule(static) reduction(+:res_loc)
	for(i = 0; i < N; i++)
		res_loc += x[i] * y[i];

	MPI_Reduce(&res_loc, &res_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	return res_glob;
}

/*
  FONCTION : CG
  DESCRIPTION : permet d'effectuer le gradient conjugue à partir d'une matrice 
  IN : matrice initiale, vecteur resultat, vecteur recherche et vecteur residu
  OUT : vecteur resultat x
  RETOUR : scalaire contenant le nombre d'iteration
*/
/*int CG(double **A, double *b, double *x, struct info_t *info)
{
	int i, j, iter = 0, maxiter = max(100.0, sqrt(n));
	double alpha, M, Mold, beta;
	double *v, *r;

	v = (double*) malloc(n * sizeof(double));
	r = (double*) malloc(n * sizeof(double));

	for(i = 0; i < n; i++)
	{
		r[i] = b[i];
		v[i] = r[i];
	}

	M = prodScal(r, r);

	while(((norme(r) / norme(b)) > DBL_EPSILON) && (iter < maxiter))
	{
		iter++;
		alpha = M / produitT(v, A);

		for(i = 0; i < n; i++)
			x[i] += alpha * v[i];

		prodMatVect(A, v, r);
		for(i = 0; i < n; i++)
			r[i] -= alpha * r[i];
		
		Mold = M;
		M = prodScal(r, r);
		beta = M / Mold;

		for(i = 0; i < n; i++)
			v[i] = r[i] + beta * v[i];
	}
	
	free(v);
	free(r);
	return iter;
}*/

void MPI_initialize(int argc, char **argv, struct info_t *info)
{
	int Q, R;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &(info->rang));
	MPI_Comm_size(MPI_COMM_WORLD, &(info->nproc));

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

double get_timer()
{
	struct timeval t;

	gettimeofday(&t, NULL);

	return t.tv_sec * 1000000 + t.tv_usec;
}

double diff_time(double end, double start)
{
	return end - start;
}

void print_time(struct info_t *info, double runtime)
{
	int i;
	double runtime_in_seconds = runtime / 1000000;
	double runtime_t[info->nproc];
	double max = 0.0;

	MPI_Gather(&runtime_in_seconds, 1, MPI_DOUBLE, &runtime_t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if(info->rang == 0)
	{
		for(i = 0; i < info->nproc; i++)
		{
			fprintf(stderr, "[%d] Temps de l'execution : %e s\n", i, runtime_t[i]);
			//if(max < runtime_t[i])
				//max = runtime_t[i];
		}
		//printf("%d \t %e\n", info->nproc, max);
	}
}
