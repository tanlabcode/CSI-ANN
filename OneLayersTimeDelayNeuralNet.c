/*This program optimizes five benchmark functions using swarm algorithm
  Asynchronous version
  Yuhui Shi, May 15, 1998
  Modified by Hiram Firpi, Nov. 2005
  Modified by Li Teng Aug. 2010*/

#include "headfile.h"
#include "global.h"
#include "mem_loc.h"
#include <signal.h>
#include "NeuralNetLib.h"

#define WindowSize 10

double OneLayersTimeDelayNeuralNet(double *fopt, double *bestpart, double *Eevolution, int iNagents, int iitermax, int indim, double imaxx, double imaxv, double ierror, double iw, int iNoBg, int iNoEn, int iNoOfFeatures, int iPpW, int *iTimeCopies, int *iDelays, int * iNNArch, double *Data, int iM )
{

   int NUMBER_OF_AGENTS, MAXITER, DIMENSION, run_no, PpW, NoOfFeatures;  /*number of runs*/
   int a=0, b=0, i=0, j=0, iter, gbest, firsttime, finish, M, N, NoEn, NoBg;
   double E_CUTOFF,  MAXV, MAXX;
   double weight, weight_up, R = RAND_MAX;
   double IRang_L, IRang_R;  /* initialization rang: left and right range */
   double minval = 0.0;
   double avep;
   double c[] = {1.4916, 1.4916};
   int *NNC, *TimeCopies = NULL, *Delays = NULL;

   time_t tt;

   NUMBER_OF_AGENTS = iNagents;
   MAXITER = iitermax;
   DIMENSION = indim;
   MAXX = imaxx;
   MAXV = imaxv;
   E_CUTOFF = ierror;
   weight = iw;
   NoBg = iNoBg;
   NoEn = iNoEn;
   NoOfFeatures = iNoOfFeatures;
   PpW = iPpW;
   TimeCopies = iTimeCopies;
   Delays = iDelays;
   NNC = iNNArch; 

   IRang_L = -5;
   IRang_R = 5;
   run_no = 1;

    N = NoBg + NoEn;
    M = iM; //M=10
   

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
            //printf("Run OneLayerTDNN\n"); // on watch
            minval = OneLayerTDNN(a, Data, N, M, NoBg, NoEn, NNC, PpW, NoOfFeatures, TimeCopies, Delays);
            //printf("Finish running OneLayerTDNN\n"); // on watch
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
               vx[b][a] = weight_up*vx[b][a] + c[0]*(rand()/R)*(pbestx[b][a] - xx[b][a]) + c[1]*(rand()/R)*(pbestx[b][gbest] - xx[b][a]);

               if (vx[b][a] > MAXV)
                  vx[b][a] = MAXV;
               else if (vx[b][a] < -MAXV)
                  vx[b][a] = -MAXV;
            }
            /* Tx allows simultaneous updates */
            for (b = 0; b < DIMENSION; b++)
               tx[b][a]=xx[b][a] + vx[b][a];
         } /******** END OF a LOOP ************************* */

/********************************************************
        Update positions
*********************************************************/
         for (a = 0;a < NUMBER_OF_AGENTS; a++)
         {
             /* Define new coordinates */
            for (b = 0; b < DIMENSION; b++)
               xx[b][a] = tx[b][a];
         } /* end a loop */

         /* In case iterations become greater than 32767 */
         if (firsttime != 1)
            if (iter == R-1)
               iter = 0;
         /*Calculate average of the fitness*/

         avep = 0.0;
         for(j = 0; j < NUMBER_OF_AGENTS; j++)
            /*avep += pbest[j]*((float)1/NUMBER_OF_AGENTS);*/
            avep += pbest[j]*((double)1/NUMBER_OF_AGENTS);
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


