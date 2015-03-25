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
void affichageMat(double **M)
{
	int i, j;

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
			printf("%.2lf   ", M[i][j]);
		printf("\n");
	}
	printf("\n");
}

/*
  PROCEDURE : affichageMatSpy
  DESCRIPTION : permet l'affichage d'une matrice sous forme spy
  IN : matrice 
  OUT : /
*/
void affichageMatSpy(double **M)
{
	int i, j;
	
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			if(fabs(M[i][j]) > DBL_EPSILON)
				printf("*  ");
			else printf("   ");
		}
		printf("\n");
	}
	printf("\n");
}

/*
  PROCEDURE : affichageVect
  DESCRIPTION : permet l'affichage d'un vecteur sur la sortie standard
  IN : vecteur
  OUT : /
*/
void affichageVect(double *V)
{
	int i;

	for(i = 0; i < n; i++)
		printf("%.2lf   ",V[i]);
	printf("\n\n");
}

/*
  PROCEDURE : poisson2D
  DESCRIPTION : permet de remplir une matrice à l'aide de l'equation de poisson 2D
  IN : matrice à remplir
  OUT : matrice 2D A
*/
void poisson2D(double **A)
{
	int i;
	
	#pragma omp parallel for
	for(i = 0; i < n; i++)
		A[i][i] = 4.0;

	#pragma omp parallel for 
	for(i = 0; i < n - 1; i++)
	{
		A[i][i + 1] = -1.0;
		A[i + 1][i] = -1.0;
	}

	#pragma omp parallel for 
	for(i = 0; i < n - nx; i++)
	{
		A[i][i + nx] = -1.0;
		A[i + nx][i] = -1.0;
	}
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
	int *map = malloc(n * (sizeof * map));
	
	for(i = 0; i < n; i++)
		map[i] = (i) % (info->nproc);
	
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			if((fabs(A[i][j]) > DBL_EPSILON) || (i == j))
				level[i][j] = 0;
			else level[i][j] = inf;

			// copie de la matrice A dans la matrice LUi
			LUi[i][j] = A[i][j]; 
		}
	}
	
	for(k = 0; k < n; k++)
	{
		if(map[k] == (info->rang))
		{ 
			for(i = k+1; i < n; i++)
			{
				if(level[i][k] <= p)
					LUi[i][k] = LUi[i][k] / LUi[k][k];
			}
		}
		
		MPI_Bcast (&A[k][k], n-k, MPI_DOUBLE, map[k], MPI_COMM_WORLD);

		for(j = k + 1; j < n; j++)
		{
			if(map[j] == (info->rang))
			{
				for(i = k + 1; i < n; i++)
				{
					if(level[i][k] <= p)
					{
						LUi[i][j] = LUi[i][j] - LUi[i][k] * LUi[k][j]; 

						if(fabs(LUi[i][j]) > DBL_EPSILON)
							level[i][j] = min(level[i][j], level[i][k] + level[k][j] + 1);
					}
				}
			}
		}
	}
	
	for(j = 0 ; j < n; j++)	
	{
		for(i = 0 ; i < n; i++) 
		{
			if(level[i][j] > p)
				LUi[i][j] = 0.0;
		}
	}
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
  IN : vecteur à calculer
  OUT : /
  RETOUR : scalaire contenant la norme du vecteur
*/
double norme(double *a)
{
	return sqrt(fabs(prodScal(a, a)));
}

/*
  FONCTION : produitT
  DESCRIPTION : permet d'effectuer le produit d'un vecteur transpose avec une matrice puis avec le vecteur (v' * a * v)
  IN : vecteur et matrice 
  OUT : /
  RETOUR : scalaire contenant le resultat du produit
*/
double produitT(double *v, double **a)
{
	double *r1, r2 = 0.0;
	
	r1 = (double*) calloc(n, sizeof(double));
	prodVect(a, v, r1);
	r2 = prodScal(r1, v);
	free(r1);
	
	return r2;
}

/*
  PROCEDURE : PCG
  DESCRIPTION : permet d'effectuer le gradient conjugue preconditionne à partir d'une matrice conditionnée (LU incomplete)
  IN : matrice initiale, vecteur recherche, vecteur resultat, matrice A sous forme ILU et vecteur residu
  OUT : vecteur resultat x
*/
int PCG(double **A, double *x, double *b, double **B, struct info_t *info)
{
	int iter = 0, i, j, maxiter = max(100.0, sqrt(n));
	double alpha, beta, M, Mold, bnrm2 = norme(b);
	double *Ap, *Br, *r, *v;

	Ap = (double*) malloc(n * sizeof(double));
	Br = (double*) malloc(n * sizeof(double));
	r = (double*) malloc(n * sizeof(double));
	v = (double*) calloc(n, sizeof(double));

	for(i = 0; i < n; i++)
		r[i] = b[i];

	prodVect(B, r, v);
	M = prodScal(r, v);

	if(bnrm2 == 0.0) 
		bnrm2 = 1.0;

	while(((norme(r) / bnrm2) > DBL_EPSILON) && (iter < maxiter))
	{
		iter++;

		for(i = 0; i < n; i++)
			Ap[i] = 0.0;
			
		prodVect(A, v, Ap);
		alpha = M / prodScal(v, Ap);

		for(i = 0; i < n; i++)
		{
			x[i] += alpha * v[i];
			r[i] -= alpha * Ap[i];
		}

		Mold = M;

		for(i = 0; i < n; i++)
			Br[i] = 0.0;
			
		prodVect(B, r, Br);
		M = prodScal(r, Br);
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
void vecteur_b(double *b)
{
	int i;
	
	#pragma omp sections nowait
	{
		#pragma omp section 
		{
			for(i = 0; i < n / 2; i++)
				b[i] = (double)i;
		}
		
		#pragma omp section 
		{
			if(n % 2 == 0)
			{
				for(i = n / 2 - 1; i >= 0; i--)
					b[n - i - 1] = (double)i;
			}else{
				for(i = n / 2; i >= 0; i--)
					b[n - i - 1] = (double)i;
			}
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
	double alpha, bnrm2 = norme(b);
	double *Ap, *Ar, *beta, *r, **v;
	
	Ap = (double*) malloc(n * sizeof(double));
	Ar = (double*) malloc(n * sizeof(double));
	r = (double*) calloc(n, sizeof(double));
	v = (double**) calloc(n, sizeof(double*));
	
	for(i = 0; i < n; i++)
		v[i] = (double*) calloc(maxiter, sizeof(double));

	prodVect(A, x, r);
	
	for(i = 0; i < n; i++)
	{
		r[i] = b[i] - r[i];
		v[i][0] = r[i];
	}

	if(bnrm2 == 0.0) 
		bnrm2 = 1.0;

	while(((norme(r) / bnrm2) > DBL_EPSILON) && (iter < maxiter))
	{
		iter++;

		for(i = 0; i < n; i++)
		{
			Ap[i] = 0.0;
			for(j = 0; j < n; j++)
				Ap[i] += A[i][j] * v[j][iter - 1];
		}
		
		alpha = prodScal(r, Ap) / prodScal(Ap, Ap);
		
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
				prodVect(A, r, Ar);
			}
			beta[k] = (-1) * prodScal(Ar, Ap) / prodScal(Ap, Ap);
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

void prodVect(double **m, double *v, double *r)
{
	int i, j;

	#pragma omp parallel shared(m,r,v) private(i,j) 
	{
		#pragma omp for schedule(static)
		for(i = 0; i < n; i++)
		{
			for(j = 0; j < n; j++)
				r[i] += m[i][j] * v[j];
		}
	}
}

double prodScal(double *x, double *y)
{
	int i;
	double resultat = 0.0;

	#pragma omp parallel for reduction(+:resultat)
	for(i = 0; i < n; i++)
		resultat += x[i] * y[i];
	
	return resultat;
}

/*
  FONCTION : CG
  DESCRIPTION : permet d'effectuer le gradient conjugue à partir d'une matrice 
  IN : matrice initiale, vecteur resultat, vecteur recherche et vecteur residu
  OUT : vecteur resultat x
  RETOUR : scalaire contenant le nombre d'iteration
*/
int CG(double **A, double *b, double *x, struct info_t *info)
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

		prodVect(A, v, r);
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
}

void MPI_initialize(int argc, char **argv, struct info_t *info)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &(info->rang));
	MPI_Comm_size(MPI_COMM_WORLD, &(info->nproc));
}
