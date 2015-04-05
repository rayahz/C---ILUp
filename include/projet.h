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
#include <time.h>
#include <mpi.h>

#define min(X, Y) (((X) < (Y)) ? (X) : (Y))
#define max(X, Y) (((X) < (Y)) ? (Y) : (X))
#define nx 2
#define ny 3
#define n (nx*ny)
#define p 4
#define inf 1024

struct info_t 
{
	int nproc;
	int rang;
	int nloc;
	int ideb;
	int ifin;
	double temps;
};

////////////////////////////////////
/* *** PROCEDURES D'AFFICHAGE *** */
////////////////////////////////////
// affichage d'une matrice
void affichageMat(double**, struct info_t *); 

// affichage d'une matrice sous forme spy
void affichageMatSpy(double**, struct info_t *); 

// affichage d'un vecteur
void affichageVect(double*, struct info_t *); 

// permet d'ecrire dans un fichier de sortie les données du programme
void ecrireFichier(struct info_t *, int, int, double*);

//////////////////////////////////
/* *** PROCEDURES DE CALCUL *** */
//////////////////////////////////
// permet de calculer la norme d'un vecteur
double norme(double*, struct info_t *); 

// permet d'effectuer le produit d'un vecteur transpose avec une matrice puis avec le vecteur (v' * a * v)
double produitT(double*, double**, struct info_t *); 

// permet de calculer le produit scalaire d'un vecteur
double prodScal(double *, double *, struct info_t *);

// permet de calculer le produit d'une matrice avec un vecteur
void prodMatVect(double **, double *, double *, struct info_t *);

//////////////////////////////////
/* *** PROCEDURES DU PROJET *** */
//////////////////////////////////
// permet d'effectuer le gradient conjugue preconditionne
int PCG(double**, double*, double*, double**, struct info_t *); 

// permet de calculer de la factorisation lu incompléte au niveau p
void ilup(double**, double**, struct info_t *); 

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
void poisson2D(double**, struct info_t *); 

// création du vecteur b
void vecteur_b(double*, struct info_t *); 

// permet d'initialiser une matrice provenant de matrix market
void matrixMarket(double **, char*);

////////////////////////////
/* *** PROCEDURES MPI *** */
////////////////////////////
// permet d'initialiser l'environnement MPI
void MPI_initialize(int, char **, struct info_t *);

/////////////////////////////////
/* *** PROCEDURES DE TEMPS *** */
/////////////////////////////////
// permet d'initialiser le timer
double get_timer();

// permet d'effectuer la difference entre le temps de fin et le temps de debut
double diff_time(double, double);

// permet d'afficher le temps d'execution de chaque processus
void print_time(struct info_t *, double );

#endif

