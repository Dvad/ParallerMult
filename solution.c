
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

#define N 8

int main (int argc, char* argv[])
{
	int rang, nbProcs, nbLignes, iLig, iCol, iProc, i;
	int A[N*N], B[N*N], C[N*N], verif[N*N];
	int *lignesA, *colonnesB, *monResultat;
	
	srand(time(NULL));
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &nbProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rang);
	
	if (rang == 0)
	{
		/* Initialisation des matrices */
		for (iLig = 0; iLig < N; iLig++)
		{
			for (iCol = 0; iCol < N; iCol++)
			{
				A[iLig * N + iCol] = rand();
				B[iLig * N + iCol] = rand();
			}
		}
		
		/* Calcul du produit, pour vérification ultérieure */
		for (iLig = 0; iLig < N; iLig++)
		{
			for (iCol = 0; iCol < N; iCol++)
			{
				verif[iLig * N + iCol] = 0;
				for (i = 0; i < N; i++)
				{
					verif[iLig * N + iCol] += A[iLig * N + i] * B[i * N + iCol];
				}
			}
		}
	}
	
	/* Combien de lignes et colonnes par processus ? */
	nbLignes = N / nbProcs;
	lignesA = malloc(nbLignes * N * sizeof(int));
	colonnesB = malloc(nbLignes * N * sizeof(int));
	monResultat = malloc(nbLignes * N * sizeof(int));
	
	/* Transmission de colonnes de B. Pour cela, transposée de B dans C. */
	if (rang == 0)
	{
		for (iLig = 0; iLig < N; iLig++)
		{
			for (iCol = 0; iCol < N; iCol++)
			{
				C[iLig * N + iCol] = B[iCol * N + iLig];
			}
		}
	}
	MPI_Scatter(C, N * nbLignes, MPI_INT, colonnesB, N * nbLignes, MPI_INT, 0, MPI_COMM_WORLD);
	
	/* Calcul des produits */
	for (iProc = 0; iProc < nbProcs; iProc++)
	{
		/* Le proc 0 transmet les nbLignes lignes suivantes de A */
		if (rang == 0)
		{
			memcpy(lignesA, &(A[iProc * nbLignes * N]), N * nbLignes * sizeof(int));
		}
		MPI_Bcast(lignesA, N*nbLignes, MPI_INT, 0, MPI_COMM_WORLD);
		
		/* Chaque proc calcul le bloc (iProc,rang) */
		for (iLig = 0; iLig < nbLignes; iLig++)
		{
			for (iCol = 0; iCol < nbLignes; iCol++)
			{
				monResultat[iProc * nbLignes + iCol * N + iLig] = 0;
				for (i = 0; i < N; i++)
				{
					monResultat[iProc * nbLignes + iCol * N + iLig] += lignesA[iLig * N + i] * colonnesB[iCol * N + i];
				}
			}
		}
	}
	
	/* Récupération du résultat : chaque proc dispose de la colonne rang dans la variable monResultat */
	MPI_Gather(monResultat, nbLignes*N, MPI_INT, C, nbLignes*N, MPI_INT, 0, MPI_COMM_WORLD);
	
	/* Le résultat est transposé, il faut le réinverser */
	for (iLig = 0; iLig < N; iLig++)
	{
		for (iCol = 0; iCol < iLig; iCol++)
		{
			i = C[iLig*N + iCol];
			C[iLig*N + iCol] = C[iCol*N + iLig];
			C[iCol*N + iLig] = i;
		}
	}
	
	/* Vérification du résultat */
	if (rang == 0)
	{
		printf("Résultat calculé :");
		for (i = 0; i < N*N; i++)
		{
				if (i % 8 == 0) printf("\n");
				printf("%d ", C[i]);
		}
		
		printf("\n\nSolution recherchée :");
		for (i = 0; i < N*N; i++)
		{
			if (i % N == 0) printf("\n");
			printf("%d ", verif[i]);
		}
		
		printf("\n\nErreurs :\n");
		for (iLig = 0; iLig < N; iLig++)
		{
			for (iCol = 0; iCol < N; iCol++)
			{
				if (C[iLig * N + iCol] != verif[iLig * N + iCol])
				{
					printf("x ");
				}
				else
				{
					printf("o ");
				}
			}
			printf("\n");
		}
	}
	
	free(lignesA);
	free(colonnesB);
	free(monResultat);
	
	MPI_Finalize();
	return EXIT_SUCCESS;
}

