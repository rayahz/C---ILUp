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
	
	for(i = 0; i < n; i++)
		A[i][i] = 4.0;
	
	for(i = 0; i < n - 1; i++)
	{
		A[i][i + 1] = -1.0;
		A[i + 1][i] = -1.0;
	}
	
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
void ilup(double **A, int **level, double **LUi)
{
	int i, j, k;
	
	// initialisation de la matrice de remplissage 
	// si A(i,j) != 0
		// alors lev(i,j) = 0
		// sinon lev(i,j) = inf
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
	
	for(i = 1; i < n; i++)
	{
		for(k = 0; k <= i-1; k++)
		{
			if(level[i][k] <= p)
			{
				LUi[i][k] = LUi[i][k] / LUi[k][k]; // calcul
			
				for(j = k + 1; j < n; j++)	
				{
					LUi[i][j] = LUi[i][j] - LUi[i][k] * LUi[k][j]; // calcul
					
					// remplissage de la matrice de niveau de remplissage
						// factorisation symbolique : Le but de cette étape est de calculer la structure de remplissage des facteurs avant l’étape de factorisation numérique.
					if(fabs(LUi[i][j]) > DBL_EPSILON)
						level[i][j] = min(level[i][j], level[i][k] + level[k][j] + 1);
				}
			}
		}
		
		// tous les elements dont le niveau de remplissage est superieur à p sont remplaces par 0

		for(j = 0 ; j < n; j++)
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
	int i;
	
	LU[0][0] = A[0][0];
	
	for(i = 1; i < n; i++)
	{
		LU[i][i - 1] = A[i][i - 1] / LU[i - 1][i - 1];
		LU[i][i] = A[i][i] - LU[i][i - 1] * A[i - 1][i];
		LU[i - 1][i] = A[i - 1][i];
	}
	
	for(i = nx; i < n; i++)
	{
		LU[i][i - nx] = A[i][i - nx] / LU[i - nx][i - nx];
		LU[i][i] = A[i][i] - LU[i][i - nx] * A[i - nx][i];
		LU[i - nx][i] = A[i - nx][i];
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
	int i;
	double resultat = 0.0;
	
	for(i = 0; i < n; i++)
		resultat = resultat + a[i] * a[i];

	return sqrt(fabs(resultat));
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
	int i, j;
	
	r1 = (double*) malloc(n * sizeof(double));
	
	for(i = 0; i < n; i++)
	{
		r1[i] = 0.0;
		for(j = 0; j < n; j++)
			r1[i] = r1[i] + v[j] * a[i][j];
		
		r2 = r2 + r1[i] * v[i];
	}
	
	free(r1);
	
	return r2;
}

void InvdiagMat(int NbElement, double **Mat, double **Minv)
{
	double **temp;
	int i, j, k;
	
	temp = (double**) calloc(NbElement, sizeof(double*));
	
	for(i = 0; i < NbElement; i++)
		temp[i] = (double*) calloc(NbElement, sizeof(double));

	for(i = 0; i < NbElement; i++)
	{
		for(j = 0; j < NbElement; j++)
		{
			temp[i][i] = 1 / Mat[i][i];
		
			if(j != i)
				temp[i][j]=-Mat[i][j]/Mat[i][i];
			
			for(k=0;k<NbElement;k++)
			{
				if(k!=i)
					temp[k][i]=Mat[k][i]/Mat[i][i];
					
				if(j!=i &&k!=i)
					temp[k][j]=Mat[k][j]-Mat[i][j]*Mat[k][i]/Mat[i][i];
			}
		}
		
		for(i=0;i<NbElement;i++)
		{
			for(j=0;j<NbElement;j++)
				Minv[i][j]=temp[i][j];
		}
	}
	
	free(temp);
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
	int Nloc = info->nloc;
	Ap = (double*) malloc(n * sizeof(double));
	Br = (double*) malloc(n * sizeof(double));
	r = (double*) malloc(n * sizeof(double));
	v = (double*) calloc(n, sizeof(double));

	for(i=0; i<n; i++)
		r[i] = b[i];

	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			v[i] += B[i][j] * r[j];
	M = prodScal(r, v, Nloc);

	if(bnrm2 == 0.0) bnrm2 = 1.0;

	while(((norme(r) / bnrm2) > DBL_EPSILON) && (iter < maxiter)){
		iter++;

		for(i=0; i<n; i++){
			Ap[i] = 0.0;
			for(j=0; j<n; j++)
				Ap[i] += A[i][j] * v[j];
		}
		alpha = M / prodScal(v, Ap, Nloc);

		for(i=0; i<n; i++){
			x[i] += alpha * v[i];
			r[i] -= alpha * Ap[i];
		}

		Mold = M;

		for(i=0; i<n; i++){
			Br[i] = 0.0;
			for(j=0; j<n; j++)
				Br[i] += B[i][j] * r[j];
		}
		M = prodScal(r, Br, Nloc);

		beta = M / Mold;

		for(i=0; i<n; i++)
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
	
	for(i = 0; i < n / 2; i++)
		b[i] = (double)i;
		
	if(n % 2 == 0)
	{
		for(i = n / 2 - 1; i >= 0; i--)
			b[n - i - 1] = (double)i;
	}else{
		for(i = n / 2; i >= 0; i--)
			b[n - i - 1] = (double)i;
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

	int Nloc = info->nloc;
	
	Ap = (double*) malloc(Nloc * sizeof(double));
	Ar = (double*) malloc(Nloc * sizeof(double));
	r = (double*) calloc(Nloc, sizeof(double));
	v = (double**) calloc(Nloc, sizeof(double*));
	
	for(i = 0; i < n; i++)
		v[i] = (double*) calloc(maxiter, sizeof(double));

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
			r[i] += A[i][j] * x[j];

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
		
		alpha = prodScal(r, Ap, Nloc) / prodScal(Ap, Ap, Nloc);
		
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
				{
					Ap[i] += A[i][j] * v[j][k];
					Ar[i] += A[i][j] * r[j];
				}
			}
			beta[k] = (-1) * prodScal(Ar, Ap, Nloc) / prodScal(Ap, Ap, Nloc);
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

double prodScal(double *v, double *v2, int N)
{
	int i;
	double result_local = 0.0, result_global;

	for(i = 0; i < N; i++)
		result_local += v[i] * v2[i];
	
	MPI_Allreduce(&result_local, &result_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	return result_global;
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
	int i, j, iter = 0, maxiter;
	double alpha, M, Mold, beta;
	double *v, *r;
	int Nloc = info->nloc;
	maxiter = max(100.0, sqrt(n));
	v = (double*) malloc(Nloc * sizeof(double));
	r = (double*) malloc(Nloc * sizeof(double));

	for(i = 0; i < Nloc; i++)
	{
		r[i] = b[i];
		v[i] = r[i];
	}

	M = prodScal(r, r, Nloc);

	while(((norme(r) / norme(b)) > DBL_EPSILON) && (iter < maxiter))
	{
		iter++;
		alpha = M / produitT(v, A);

		for(i = 0; i < n; i++)
			x[i] += alpha * v[i];

		for(i = 0; i < n; i++)
		{
			for(j = 0; j < n; j++)
				r[i] -= alpha * A[i][j] * v[j];
		}
		
		Mold = M;
		M = prodScal(r, r, Nloc);
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
	int Q, R;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &(info->rang));
	MPI_Comm_size(MPI_COMM_WORLD, &(info->nproc));
	
	Q = n / info->nproc;
	R = n % info->nproc;

	info->ntot = n;

	if (info->rang < R) 
	{
		info->nloc = Q+1;
		info->ideb = info->rang * (Q+1);
		info->ifin = info->ideb + info->nloc;
	} else{
		info->nloc = Q;
		info->ideb = R * (Q+1) + (info->rang - R) * Q;
		info->ifin = info->ideb + info->nloc;
	}
}
