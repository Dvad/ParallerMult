
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "time.h"
#include <stdbool.h>

 
#define NPRINT 8
#define ETIQUETTE 100

#define N  1024
#define ASYNCHRONE_AVEC_TEST


int main (int argc, char* argv[])
{
        const int PAS=N*N/NPRINT;
	int rang, nbProcs, nbLignes, iLig, iCol, iProc, i;

	float * A = malloc(N*N*sizeof(float));
	float * B = malloc(N*N*sizeof(float));
	float * C = malloc(N*N*sizeof(float));
	float * verif = malloc(N*N*sizeof(float));

	float *lignesA, *lignesAbis, *colonnesB, *monResultat;
	double temps_debut,temps_ecoule,temps_max;
	MPI_Request* send_requests;
	MPI_Request* receive_requests;
	MPI_Status statut;
	bool* processed;
	int nbProcessed=0;
	int flag;
	
	srand(time(NULL));
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &nbProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rang);

	send_requests = (MPI_Request*) malloc(nbProcs*sizeof(MPI_Request));
	receive_requests = (MPI_Request*) malloc(nbProcs*sizeof(MPI_Request));
	processed = (bool*) malloc(nbProcs*sizeof(bool));

	for(i=0; i<nbProcs; i++)
	  processed[i]=false;

	if (rang == 0)
	{
	  printf("Matrice %dx%d\n",N,N);

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
	lignesA = malloc(nbLignes * N * sizeof(float));
	lignesAbis = malloc(nbLignes * N * sizeof(float));
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
	
	/* Transmission de colonnes de A : */
	MPI_Scatter(A, N * nbLignes, MPI_FLOAT, lignesA, N * nbLignes, MPI_FLOAT, 0, MPI_COMM_WORLD);


	temps_debut = MPI_Wtime();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef SYNCHRONE
	if (rang == 0)
	  printf("Synchrone\n");

	/* Calcul des produits */
	for (iProc = 0; iProc < nbProcs; iProc++)
	{

		/* Le proc concerné transmet les nbLignes lignes suivantes de A */
		if (rang == iProc)
		{
			memcpy(lignesAbis, lignesA, N* nbLignes * sizeof(float));
		}
		MPI_Bcast(lignesAbis, N*nbLignes, MPI_FLOAT, iProc, MPI_COMM_WORLD);
		
		/* Chaque proc calcul le bloc (iProc,rang) */
		for (iLig = 0; iLig < nbLignes; iLig++)
		{
			for (iCol = 0; iCol < nbLignes; iCol++)
			{
				monResultat[iProc * nbLignes + iCol * N + iLig] = 0;
				for (i = 0; i < N; i++)
				{
					monResultat[iProc * nbLignes + iCol * N + iLig] += lignesAbis[iLig * N + i] * colonnesB[iCol * N + i];
				}
			}
		}
	}

#endif

	////////////////////////////////////////////////////////////////////////////////////



#ifdef ASYNCHRONE

	if (rang == 0)
	  printf("Asynchrone\n");

	/* Chaque processus envoie sa ligne de A à tous les autres de manière asynchrone : */
	
	for(iProc=0; iProc<nbProcs; iProc++)
	  {
	    if(rang != iProc)
	      {
		MPI_Issend(lignesA,N*nbLignes,MPI_FLOAT,iProc,ETIQUETTE,MPI_COMM_WORLD,&send_requests[iProc]);	
		MPI_Irecv(&(A[N*nbLignes*iProc]),N*nbLignes,MPI_FLOAT,iProc,ETIQUETTE,MPI_COMM_WORLD,&receive_requests[iProc]);
	      }
	  }


	/* Calcul du bloc diagonal */
	for (iLig = 0; iLig < nbLignes; iLig++)
	  {
	    for (iCol = 0; iCol < nbLignes; iCol++)
	      {
		monResultat[rang * nbLignes + iCol * N + iLig] = 0;
		for (i = 0; i < N; i++)
		  {
		    monResultat[rang * nbLignes + iCol * N + iLig] += lignesA[iLig * N + i] * colonnesB[iCol * N + i];
		  }
	      }
	  }

	/* Calcul du reste de la matrice */

	for(iProc=0; iProc<nbProcs; iProc++)
	  {
	    if(rang != iProc)
	      {
						
		MPI_Wait(&receive_requests[iProc], &statut);
							
		for (iLig = 0; iLig < nbLignes; iLig++)
		  {
		    for (iCol = 0; iCol < nbLignes; iCol++)
		      {
			monResultat[iProc * nbLignes + iCol * N + iLig] = 0;
			for (i = 0; i < N; i++)
			  {
			    monResultat[iProc * nbLignes + iCol * N + iLig] += A[(nbLignes*iProc + iLig) * N + i] * colonnesB[iCol * N + i];
			  }
		      }
		  }		
	      }
	  }
	
#endif
	
	////////////////////////////////////////////////////////////////////////////////////



#ifdef ASYNCHRONE_AVEC_TEST
	if (rang == 0)
	  printf("Asynchrone avec test \n");

	/* Chaque processus envoie sa ligne de A à tous les autres de manière asynchrone : */
	
	for(iProc=0; iProc<nbProcs; iProc++)
	  {
	    if(rang != iProc)
	      {
		MPI_Issend(lignesA,N*nbLignes,MPI_FLOAT,iProc,ETIQUETTE,MPI_COMM_WORLD,&send_requests[iProc]);	
		MPI_Irecv(&(A[N*nbLignes*iProc]),N*nbLignes,MPI_FLOAT,iProc,ETIQUETTE,MPI_COMM_WORLD,&receive_requests[iProc]);
	      }
	  }


	/* Calcul du bloc diagonal */
	for (iLig = 0; iLig < nbLignes; iLig++)
	  {
	    for (iCol = 0; iCol < nbLignes; iCol++)
	      {
		monResultat[rang * nbLignes + iCol * N + iLig] = 0;
		for (i = 0; i < N; i++)
		  {
		    monResultat[rang * nbLignes + iCol * N + iLig] += lignesA[iLig * N + i] * colonnesB[iCol * N + i];
		  }
	      }
	  }

	/* Calcul du reste de la matrice */


	while(nbProcessed < (nbProcs-1))
	  {
	    for(iProc=0; iProc<nbProcs; iProc++)
	      {
		if(rang != iProc && ! processed[iProc])
		  {
		    MPI_Test(&receive_requests[iProc], &flag, &statut);

		    if(flag == 1)
		      {
			MPI_Wait(&receive_requests[iProc], &statut);
			
			processed[iProc]=true;
			nbProcessed++;

			for (iLig = 0; iLig < nbLignes; iLig++)
			  {
			    for (iCol = 0; iCol < nbLignes; iCol++)
			      {
				monResultat[iProc * nbLignes + iCol * N + iLig] = 0;
				for (i = 0; i < N; i++)
				  {
				    monResultat[iProc * nbLignes + iCol * N + iLig] += A[(nbLignes*iProc + iLig) * N + i] * colonnesB[iCol * N + i];
				  }
			      }
			  }		
		      }
		  }
		
	      }
	  }
#endif
	

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	temps_ecoule = MPI_Wtime() - temps_debut;

	/* Récupération du résultat : chaque proc dispose de la colonne rang dans la variable monResultat */
	MPI_Gather(monResultat, nbLignes*N, MPI_FLOAT, C, nbLignes*N, MPI_FLOAT, 0, MPI_COMM_WORLD);

	/* Récupération du temps max d'exécution pour chaque processus : */
	MPI_Reduce(&temps_ecoule,&temps_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	if(rang == 0)
	  printf("Temps écoulé = %f sec\n", temps_ecoule);


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
	/*
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
	*/
	free(lignesA);
	free(colonnesB);
	free(monResultat);
	free(send_requests);
	free(receive_requests);
	free(processed);
	
	MPI_Finalize();
	return EXIT_SUCCESS;
}

