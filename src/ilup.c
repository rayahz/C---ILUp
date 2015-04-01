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

