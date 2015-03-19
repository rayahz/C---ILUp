/*
	AUTHORS : Walid HAMMACHE, Yanis KHORSI, Rayhana ZIARA
	Students in M1 IHPS 2014 / 2015
	Project regarding incomplet factorisation (ILU)
	LANGUAGE USED : C
*/
#ifndef __ProjetFunctions__
#define __ProjetFunctions__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <mpi.h>

#define min(X, Y) (((X) < (Y)) ? (X) : (Y))
#define max(X, Y) (((X) < (Y)) ? (Y) : (X))
#define nx 8
#define ny 4
#define n (nx*ny)
#define p 2
#define inf 1024

struct info_t 
{
	int nproc;
	int rang;

	int ntot;
	int nloc;

	int ideb;
	int ifin;
};


////////////////////////////////////
/* *** PROCEDURES D'AFFICHAGE *** */
////////////////////////////////////
// affichage d'une matrice
void affichageMat(double**); 

// affichage d'une matrice sous forme spy
void affichageMatSpy(double**); 

// affichage d'un vecteur
void affichageVect(double*); 

//////////////////////////////////
/* *** PROCEDURES DE CALCUL *** */
//////////////////////////////////
// permet de calculer la norme d'un vecteur
double norme(double*); 

// permet d'effectuer le produit d'un vecteur transpose avec une matrice puis avec le vecteur (v' * a * v)
double produitT(double*, double**); 

// permet de calculer le produit scalaire d'un vecteur
double prodScal(double *, double *, int);

//////////////////////////////////
/* *** PROCEDURES DU PROJET *** */
//////////////////////////////////
// permet d'effecteur le gradient conjugue preconditionne
int PCG(double**, double*, double*, double**, struct info_t *); 

// permet de calculer de la factorisation lu incompléte au niveau p
void ilup(double**, int**, double**); 

// permet d'effectuer la factorisation LU complete
void LUfact(double**, double**);

// permet d'effectuer le gradient conjugue 
int CG(double **, double *, double *, struct info_t *);

// permet d'effectuer le gradient conjugue residuel
int CGR(double **, double *, double *, struct info_t *);

///////////////////////////////////////
/* *** PROCEDURES DE REMPLISSAGE *** */
///////////////////////////////////////
// création de la matrice de poisson2D
void poisson2D(double**); 

// création du vecteur b
void vecteur_b(double*); 

// permet d'initialiser une matrice provenant de matrix market
void matrixMarket(double **, char*);

///////////////////////////////////////
/* *** PROCEDURES DE REMPLISSAGE *** */
///////////////////////////////////////
void MPI_initialize(int, char **, struct info_t *);

#endif
