/*This program optimizes five benchmark functions using swarm algorithm
  Asynchronous version
  Yuhui Shi, May 15, 1998
  Modified by Hiram Firpi, Nov. 14, 2005
  Modified by A.F. Aug. 11, 2009*/

#include "headfile.h"
#include "global.h"
#include "mem_loc.h"
#include <signal.h>
#include "NeuralNetLib.h"

int main() 
{
	int NUMBER_OF_AGENTS, MAXITER, DIMENSION, run_no, PpW, NoOfFeatures;  /*number of runs*/
	int a=0, b=0, i=0, j=0, iter, gbest, firsttime, finish, M, N, NoEn, NoBg;
	float E_CUTOFF,  MAXV, MAXX;
	float weight, weight_up, R = RAND_MAX;
	float IRang_L, IRang_R;  /* initialization rang: left and right range */
	double minval = 0.0;
	double avep, *fopt, *bestpart, *Eevolution;
    double *Data;
    float c[] = {1.4916, 1.4916};
    int *NNC, *TimeCopies = NULL, *Delays = NULL;
    int iTimeCopies[2]={1,1};
    int iDelays[2]={9,0};
    int iNNC[3]={10,2,1};
    double iData[10][2588];

	time_t tt;   
   /* Creating outputs */
    //plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    //fopt = mxGetPr(plhs[0]);
       
   NUMBER_OF_AGENTS = 10;
   MAXITER = 400;
   DIMENSION = 22;
   MAXX = 10;
   MAXV = 0.5;
   E_CUTOFF = 1.0000e-04;
   weight = 0.9;
   NoBg = 2510;
   NoEn = 251;
   NoOfFeatures = 1;
   PpW = 10;
   M = 10;
   TimeCopies = (int *)iTimeCopies;
   Delays = (int *)iDelays;
   NNC = (int *)iNNC;
   for(i=0;i<10;i++)
      for(j=0;j<2588;j++)
         iData[i][j]=0;
   Data = (double *)iData;
   N = NoBg + NoEn;
   bestpart = malloc(DIMENSION * sizeof(float)); 
   Eevolution = malloc(MAXITER * sizeof(float));        




	IRang_L = -5;
	IRang_R = 5;
	run_no = 1;

    LoadingData(NoBg, NoEn, NNC);


	/*Allocating Memory*/
	FVectorAllocate(&pbest,             NUMBER_OF_AGENTS);
	FVectorAllocate(&maxx,   DIMENSION);
	FMatrixAllocate(&vx,     DIMENSION, NUMBER_OF_AGENTS);
	FMatrixAllocate(&xx,     DIMENSION, NUMBER_OF_AGENTS);
	FMatrixAllocate(&tx,     DIMENSION, NUMBER_OF_AGENTS);
	FMatrixAllocate(&pbestx, DIMENSION, NUMBER_OF_AGENTS);
	
	for (a = 0; a < DIMENSION; a++)
		maxx[a] = MAXX; /* range of xx[] */
	
	time(&tt);
	/*printf("begin time: %s\n", ctime(&tt));*/
	
	/*loop for runs*/
	for (i=0; i<run_no; i++)
	{
		firsttime = 1; /*first iteration of this run*/
		iter = 0;
		gbest = 0;               /*initialy assume the first particle as the gbest*/
	
		/* This loop initializes the individual agents  for each run */
		for (a = 0; a < NUMBER_OF_AGENTS; a++)
		{
			for (b = 0; b < DIMENSION; b++)
			{
				xx[b][a] = ((IRang_R - IRang_L)*(rand()/R) + IRang_L);
				pbestx[b][a] = xx[b][a];
				vx[b][a] = MAXV*(rand()/R);

				if ((rand()/R) > 0.5)
					vx[b][a] = -vx[b][a];
			}
		}

		finish = 0;

		/* Main Work Loop for each run here */
		do
		{
			iter++;
			if (iter > R) iter = 0;   /* so it doesn't crash the data type */
		
			/*update inertia weight*/
			weight_up = (weight-0.4)*(MAXITER - iter)/MAXITER + 0.4;    /*time variant weight, linear from weight to 0.4*/

			/*weight_up=weight;		//constant inertia weight*/
			
			for(a = 0; a < NUMBER_OF_AGENTS; a++)
			{
				/*Evaluating Neural Network*/
				minval = OneLayerTDNN(a, Data, N, M, NoBg, NoEn, NNC, PpW, NoOfFeatures, TimeCopies, Delays);
				
				if (firsttime == 1)
                {
					pbest[a] = minval;
                    *fopt = minval;
                }

				if (minval < pbest[a])
				{
					pbest[a] = minval;
					for (b = 0; b < DIMENSION; b++)
						pbestx[b][a] = xx[b][a];

					if (pbest[a] < pbest[gbest])
                    {
						gbest = a;
                        *fopt = minval;
                    }
				}
				
				/* asynchronous version */
				for (b = 0; b < DIMENSION; b++)
				{
					vx[b][a] = weight_up*vx[b][a] + c[0]*(rand()/R)*(pbestx[b][a] - xx[b][a]) +	c[1]*(rand()/R)*(pbestx[b][gbest] - xx[b][a]);
					
					if (vx[b][a] > MAXV)
						vx[b][a] = MAXV;
					else if (vx[b][a] < -MAXV)
						vx[b][a] = -MAXV;
				}

				/* Tx allows simultaneous updates */
				for (b = 0; b < DIMENSION; b++)
					tx[b][a]=xx[b][a] + vx[b][a];
			}/******** END OF a LOOP ************************* */

			/********************************************************
				Update positions
			*********************************************************/
			for (a = 0;a < NUMBER_OF_AGENTS; a++)
			{
				/* Define new coordinates */
				for (b = 0; b < DIMENSION; b++)
					xx[b][a] = tx[b][a];

			}/* end a loop */

			/* In case iterations become greater than 32767 */

			if (firsttime != 1)
				if (iter == R-1)
					iter = 0;

			/*Calculate average of the fitness*/
			avep = 0.0;
			for(j = 0; j < NUMBER_OF_AGENTS; j++)
				avep += pbest[j]*((float)1/NUMBER_OF_AGENTS);

			
			/*printf("Generation: %d   Best Fitness: %.5f\n", iter, pbest[gbest]);*/
            Eevolution[iter-1] = pbest[gbest];
             
			/*fprintf(fp, "%d %.5f %.5f\n", iter, pbest[gbest], avep);*/
			
			/* Terminate on criterion */
			if ((pbest[gbest] <= E_CUTOFF) || (iter >= MAXITER))
			{
				/*printf("%d run finished!\n", i);*/
				
				for (j = 0; j < DIMENSION; j++)
					bestpart[j] = pbestx[j][gbest];
				finish = 1;
            }
			firsttime = 0;

		}while(! finish); /* End of do-loop */
	} /*Run loop*/
	time(&tt);
	/*printf("end time: %s\n",ctime(&tt));*/

	FreeMem();
    free(pbest);
    free(maxx);
    FMatrixFree(xx, DIMENSION);
    FMatrixFree(vx, DIMENSION);
    FMatrixFree(pbestx, DIMENSION);
    FMatrixFree(tx, DIMENSION);

    return;
    
}

