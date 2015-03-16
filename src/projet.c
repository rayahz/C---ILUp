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

/*
  PROCEDURE : PCG
  DESCRIPTION : permet d'effectuer le gradient conjugue preconditionne à partir d'une matrice conditionnée (LU incomplete)
  IN : matrice initiale, vecteur recherche, vecteur resultat, matrice A sous forme ILU et vecteur residu
  OUT : vecteur resultat x
*/
void PCG(double **A, double *x, double *b, double **B, double *r)
{
	int i, j, iter = 0, maxiter;
	double alpha, M, Mold, beta, *m;
	
	m = (double*) malloc(n * sizeof(double));
	
	// initialisation du maximum d'iterations
	maxiter = max(100.0, sqrt(n));
	
	// initialisation du vecteur residu
	for(i = 0; i < n; i++)
		r[i] = b[i];
	
	// calcul = r' * B * r
	M = produitT(r, B);
	
	for(i = 0; i < n; i++)
	{
		m[i] = 0.0;
		for(j = 0; j < n; j++)
			m[i] = m[i] + B[i][j] * r[j];
	}

	while(((norme(r) / norme(b)) > DBL_EPSILON) && (iter < maxiter))
	{
		iter++;
		// calcul = M / (m' * A * m)
		alpha = M / produitT(m, A); 
		
		for(i = 0; i < n; i++)
		{
			// calcul du vecteur solution
			x[i] = x[i] + alpha * m[i]; 
	
			// calcul du vecteur residu
			for(j = 0; j < n; j++)
				r[i] = r[i] - alpha * A[i][j] * m[j]; 
		}
		
		Mold = M;
		M = produitT(r, B); // calcul = r' * B * r
		beta = M / Mold;
		
		for(i = 0; i < n; i++)
		{
			for(j = 0; j < n; j++)
				m[i] = B[i][j] * r[j];
			m[i] = m[i] + beta * m[i];
		}
	}

	printf("nb iteration = %d\n", iter);
	free(m);
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
int CGR(double **A, double *x, double *b)
{
	int i, j, l, iter = 0, maxiter;
	double erreur, alpha, beta, tmp1 = 0, tmp2 = 0, bnrm2 = norme(b); 
	double *v, *Ap, *Ar, *residu, *betaP;
	
	//maxiter = max(100.0, sqrt(n));
	maxiter = 1000;
	residu = (double*) calloc(n, sizeof(double));
	v = (double*) calloc(n, sizeof(double));
	Ar = (double*) calloc(n, sizeof(double));
	Ap = (double*) calloc(n, sizeof(double));
	
	
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
			residu[i] += A[i][j] * x[j];
		
		residu[i] = b[i] - residu[i];
		v[i] = residu[i];
	}
	
	if(bnrm2 == 0.0)
		bnrm2 = 1.0;
		
	erreur = norme(residu) / bnrm2;
	
	while((erreur >= 0.000001) && (iter < maxiter))
	{
		for(i = 0; i < n; i++)
		{
			for(j = 0; j < n; j++)
				Ap[i] += A[i][j] * v[j];
		}
		
		for(j = 0; j < n; j++)
			tmp1 += residu[j] * Ap[j];

		alpha = tmp1 / prodScal(Ap);

		for(j = 0; j < n; j++)
		{	
			x[j] += alpha * v[j];
			residu[j] -= alpha * Ap[j];
		}
		betaP = (double*) calloc(iter+1, sizeof(double));
		for(l = 0; l < iter+1; l++)
		{
			for(i = 0; i < n; i++)
			{
				Ap[i] = 0.0;
				Ar[i] = 0.0;
				
				for(j = 0; j < n; j++)
				{
					Ar[i] += A[i][j] * residu[j];
					Ap[i] += A[i][j] * v[j];
				}
			}
		
			for(j = 0; j < n; j++)
				tmp2 += Ar[j] * Ap[j];

			betaP[l] = tmp2 / prodScal(Ap);
		}
		
		for(l = 0; l < iter+1; l++)
		{
			for(i = 0; i < n; i++)
				v[i] += betaP[l] * v[i];
		}
		
		for(i = 0; i < n; i++)
			v[i] += residu[i];

		iter++;
		erreur = norme(residu) / bnrm2;
		free(betaP);
	}
	
	free(v);
	free(Ap);
	free(Ar);
	free(residu);
	//free(betaP);
	
	return iter;
}

/*
  FONCTION : prodScal
  DESCRIPTION : permet de calculer le produit scalaire d'un vecteur
  IN : vecteur à calculproduiter
  OUT : /
  RETOUR : scalaire contenant le resultat
*/
double prodScal(double *v)
{
	int i;
	double resultat = 0.0;

	for(i = 0; i < n; i++)
		resultat += v[i] * v[i];

	return resultat;
}

/*
  FONCTION : CG
  DESCRIPTION : permet d'effectuer le gradient conjugue à partir d'une matrice 
  IN : matrice initiale, vecteur resultat, vecteur recherche et vecteur residu
  OUT : vecteur resultat x
  RETOUR : scalaire contenant le nombre d'iteration
*/
int CG(double **A, double *b, double *x, double *r)
{
	int i, j, iter = 0, maxiter;
	double alpha, M, Mold, beta;
	double *v;

	maxiter = max(100.0, sqrt(n));
	v = (double*) malloc(n * sizeof(double));

	for(i = 0; i < n; i++)
	{
		r[i] = b[i];
		v[i] = r[i];
	}

	M = prodScal(r);

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
		M = prodScal(r);
		beta = M / Mold;

		for(i = 0; i < n; i++)
			v[i] = r[i] + beta * v[i];
	}
	
	free(v);
	return iter;
}
