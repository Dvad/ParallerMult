
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define N  1024 
#define NPRINT 4 

int main (int argc, char* argv[])
{
        const int PAS=N*N/NPRINT;
	int rang, nbProcs, nbLignes, iLig, iCol, iProc, i;
	float * A = malloc(N*N*sizeof(float));
	float * B = malloc(N*N*sizeof(float));
	float * C = malloc(N*N*sizeof(float));
	float * verif = malloc(N*N*sizeof(float));
	float* lignesA[2], *colonnesB, *monResultat;
	double temps_debut, temps_fin;
	MPI_Status statut;
	
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
				A[iLig * N + iCol] = rand()%16;
				B[iLig * N + iCol] = rand()%16;
				//A[iLig * N + iCol] = 1;
				//B[iLig * N + iCol] = 1;
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
	lignesA[0] = (float*) malloc(nbLignes * N * sizeof(float));
	lignesA[1] = (float*) malloc(nbLignes * N * sizeof(float));
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
	//MPI_Scatter(A, N * nbLignes, MPI_FLOAT, lignesA, N * nbLignes, MPI_FLOAT, 0, MPI_COMM_WORLD);
	if(rang==0){	
		memcpy(lignesA[1], &(A[0]), N * nbLignes * sizeof(float));
	}	
	MPI_Bcast(lignesA[1], N*nbLignes, MPI_FLOAT, 0, MPI_COMM_WORLD);
	//MPI_Barrier(MPI_COMM_WORLD);

	/* Calcul des produits */
	
	int const precedent = (rang + nbProcs - 1) % nbProcs;
   	int const suivant = (rang + 1) % nbProcs;
	int const dernier = nbProcs - 1;
	int nbr_element = nbLignes * N ;
	MPI_Request Reqenv, Reqrecept;
   	int actuel = 0;
	temps_debut= MPI_Wtime();
	for (iProc = 0; iProc < nbProcs; iProc++)
	{
		/* Le proc 0 transmet les nbLignes lignes suivantes de A */
		if (rang == 0)
		{
			memcpy(lignesA[actuel], &(A[nbLignes * N * (iProc+1)]), N * nbLignes * sizeof(float));
			//MPI_Issend(lignesA, N * nbLignes, MPI_FLOAT, 0, 1000, MPI_COMM_WORLD, &Reqenv);
			MPI_Isend(lignesA[actuel], N * nbLignes, MPI_FLOAT, 0, 1000, MPI_COMM_WORLD, &Reqenv);
			MPI_Isend(lignesA[actuel], N * nbLignes, MPI_FLOAT, 1, 1000, MPI_COMM_WORLD, &Reqenv);
			MPI_Isend(lignesA[actuel], N * nbLignes, MPI_FLOAT, 2, 1000, MPI_COMM_WORLD, &Reqenv);
			MPI_Isend(lignesA[actuel], N * nbLignes, MPI_FLOAT, 3, 1000, MPI_COMM_WORLD, &Reqenv);
		}
		//else
		//{
			//MPI_Bcast(lignesA, N*nbLignes, MPI_INT, 0, MPI_COMM_WORLD);
			//MPI_Issend((lignesA+(actuel*nbr_element)), N * nbLignes, MPI_FLOAT, precedent, 1000, MPI_COMM_WORLD, &Reqenv);
			MPI_Irecv(lignesA[actuel] , N * nbLignes, MPI_FLOAT, 0, 1000, MPI_COMM_WORLD, &Reqrecept);
		//}
	
		
		/* Chaque proc calcul le bloc (iProc,rang) */

			for (iLig = 0; iLig < nbLignes; iLig++)
			{
				for (iCol = 0; iCol < nbLignes; iCol++)
				{
					monResultat[iProc * nbLignes + iCol * N + iLig] = 0;
					for (i = 0; i < N; i++)
					{	
						//if ((rang == 0)&&(iProc!=0))
						//{
						//	{monResultat[(iProc) * nbLignes + iCol * N + iLig] += lignesA[actuel][iLig * N + i] * colonnesB[iCol * N + i];}							
							
						//}
						//else
						//{
							monResultat[(iProc) * nbLignes + iCol * N + iLig] += lignesA[!actuel][iLig * N + i] * colonnesB[iCol * N + i];
						//}
					}
				}
			}
		if (rang == 0)
		{
			MPI_Wait(&Reqenv, &statut);
		}
		if (rang != 0)
		{
			MPI_Wait(&Reqrecept, &statut);
			actuel=!actuel;
		}
	}
	/*
	for (iProc = 0; iProc < nbProcs; iProc++)
	{
		 MPI_Issend((lignesA+(actuel*nbr_element)), N * nbLignes, MPI_FLOAT, precedent, 1000, MPI_COMM_WORLD, &Reqenv);
        	 MPI_Irecv(lignesA +((1-actuel)*nbr_element), N * nbLignes, MPI_FLOAT, suivant, 1000, MPI_COMM_WORLD, &Reqrecept);
		//printf("\n\nle précedent du processus %d est %d  :",rang,precedent);
		//printf("\n\nle suivant du processus %d est %d  :",rang,suivant);
		
		
		for (iLig = 0; iLig < nbLignes; iLig++)
		{
			for (iCol = 0; iCol < nbLignes; iCol++)
			{
				monResultat[((iProc+rang)%4) * nbLignes + iCol * N + iLig] = 0;
				for (i = 0; i < N; i++)
				{
					monResultat[((iProc+rang)%4)  * nbLignes + iCol * N + iLig] += lignesA[iLig * N + i+actuel*nbr_element] * colonnesB[iCol * N + i];
				}
			}
		}
	MPI_Wait(&Reqenv, &statut);
	MPI_Wait(&Reqrecept, &statut);
        actuel = 1 - actuel;
	}*/
	
	temps_fin= MPI_Wtime()- temps_debut;
	MPI_Barrier(MPI_COMM_WORLD);
	if(rang==0)
		printf("Temps mis pour calculer est  :%6.3f secondes \n", rang, temps_fin);
	
	/* Récupération du résultat : chaque proc dispose de la colonne rang dans la variable monResultat */
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
				printf("%f ", C[i]);
		}
		
		printf("\n\nSolution recherchée :");
		for (i = 0; i < N*N; i+=PAS)
		{
			if (i % N == 0) printf("\n");
			printf("%f ", verif[i]);
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
	
	free(lignesA[0]);
	free(lignesA[1]);
	free(colonnesB);
	free(monResultat);
	
	MPI_Finalize();
	return EXIT_SUCCESS;
}

