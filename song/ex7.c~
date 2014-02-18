#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define N  1024
#define NPRINT 4 

int main (int argc, char* argv[])
{
        const int PAS=N*N/NPRINT;
	int rang, nbProcs, nbLignes, iLig, iCol, iProc, i, aux;
	float * A = malloc(N*N*sizeof(float));
	float * B = malloc(N*N*sizeof(float));
	float * C = malloc(N*N*sizeof(float));
	float * verif = malloc(N*N*sizeof(float));
	float *lignesA, *colonnesB, *monResultat, *m, *lignesA1[3];//*lignesA1, *lignesA2, *lignesA3;
	double inittime,totaltime,recvtime;
	
	MPI_Status status[3];
	MPI_Request send_request[3], recv_request[3];
	
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
				A[iLig * N + iCol] = rand()%100;
				B[iLig * N + iCol] = rand()%100;
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
	lignesA = malloc(nbLignes * N * sizeof(float));
	//lignesA1 = malloc(nbLignes * N * sizeof(float));
	lignesA1[0] = malloc(nbLignes * N * sizeof(float));
	lignesA1[1] = malloc(nbLignes * N * sizeof(float));
	lignesA1[2] = malloc(nbLignes * N * sizeof(float));
	colonnesB = malloc(nbLignes * N * sizeof(float));
	monResultat = malloc(nbLignes * N * sizeof(float));
	
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
	
	
	MPI_Scatter(C, N * nbLignes, MPI_FLOAT, colonnesB, N * nbLignes, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Scatter(A, N * nbLignes, MPI_FLOAT, lignesA, N * nbLignes, MPI_FLOAT, 0, MPI_COMM_WORLD);
	
	MPI_Isend(lignesA,N * nbLignes, MPI_FLOAT, (rang+1)%4,0,MPI_COMM_WORLD,&send_request[0]);
	MPI_Isend(lignesA,N * nbLignes, MPI_FLOAT, (rang+2)%4,0,MPI_COMM_WORLD,&send_request[1]); 
	MPI_Isend(lignesA,N * nbLignes, MPI_FLOAT, (rang+3)%4,0,MPI_COMM_WORLD,&send_request[2]);
	
 	MPI_Irecv(lignesA1[0],N * nbLignes, MPI_FLOAT, (rang+3)%4,0,MPI_COMM_WORLD,&recv_request[0]);
	MPI_Irecv(lignesA1[1],N * nbLignes, MPI_FLOAT, (rang+2)%4,0,MPI_COMM_WORLD,&recv_request[1]);
	MPI_Irecv(lignesA1[2],N * nbLignes, MPI_FLOAT, (rang+1)%4,0,MPI_COMM_WORLD,&recv_request[2]);
	
	inittime = MPI_Wtime();
	
	aux = rang;
	for (iLig = 0; iLig < nbLignes; iLig++)
		{
			for (iCol = 0; iCol < nbLignes; iCol++)
			{
				monResultat[aux * nbLignes + iCol * N + iLig] = 0;
				for (i = 0; i < N; i++)
				{
					monResultat[aux * nbLignes + iCol * N + iLig] += lignesA[iLig * N + i] * colonnesB[iCol * N + i];
				}
			}
		}
	
	
	/* Calcul des produits */
	for (iProc = 0; iProc < nbProcs-1; iProc++)
	{	
		  MPI_Wait(&send_request[iProc],&status[iProc]);
		  MPI_Wait(&recv_request[iProc],&status[iProc]);
		  m = lignesA;
 		  lignesA = lignesA1[iProc];
 		  lignesA1[iProc] = m;
		  
		aux=(rang + 3 - iProc )%4;
		/* Chaque proc calcul le bloc (iProc,rang) */
		for (iLig = 0; iLig < nbLignes; iLig++)
		{
			for (iCol = 0; iCol < nbLignes; iCol++)
			{
				monResultat[aux * nbLignes + iCol * N + iLig] = 0;
				for (i = 0; i < N; i++)
				{
					monResultat[aux * nbLignes + iCol * N + iLig] += lignesA[iLig * N + i] * colonnesB[iCol * N + i];
				}
			}
		}
		
// 		if(iProc < 3)
// 		{
// 		  m = lignesA;
// 		  lignesA = lignesA1[iProc];
// 		  lignesA1[iProc] = m;
// 		}
	}
	/* Récupération du résultat : chaque proc dispose de la colonne rang dans la variable monResultat */
	MPI_Barrier(MPI_COMM_WORLD);	
	recvtime = MPI_Wtime();
	totaltime = recvtime - inittime;
	if(rang==0)
		printf ("Your calculations took %.6lf in seconds to run, rANK: %d.\n", totaltime, rang );
		
	MPI_Gather(monResultat, nbLignes*N, MPI_FLOAT, C, nbLignes*N, MPI_FLOAT, 0, MPI_COMM_WORLD);
	
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
		for (i = 0; i < N*N; i+=PAS)
		{
				if (i % 8 == 0) printf("\n");
				printf("%d ", C[i]);
		}
		
		printf("\n\nSolution recherchée :");
		for (i = 0; i < N*N; i+=PAS)
		{
			if (i % N == 0) printf("\n");
			printf("%d ", verif[i]);
		}
		
		printf("\n\nErreurs :\n");
		for (iLig = 0; iLig < N; iLig+=PAS/N)
		{
			for (iCol = 0; iCol < N; iCol+=PAS/N)
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

