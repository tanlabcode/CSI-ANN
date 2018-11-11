#include <stdio.h>
#include <string.h> 
#include <stdlib.h>
#include <math.h>


#define MaxHistNo 25 // permitted maximum number of histones would be MaxHistNo = 25
#define MaxSize MaxHistNo*2 
#define WindowSize 10 //size of the window is 10
#define MaxTrainingSample 10000  //maximum number of training sample = 10000
#define Fold 16    // fold of leave-one-out crossvaliation 

#include "inverse.h"  // only works on matrix with size less than MaxSize * MaxSize
#include "svd2.h"
#include "randperm.h"


#include "global.h"
#include "mem_loc.h" 
#include "extern.h"
#include "NeuralNetLib.h"
#include "OneLayersTimeDelayNeuralNet.h"
#include "TDNNEva.h"
 
#define PRINTN printf("\n");

main(){
   double x;
   int no_row=0, no_column=0, i, j, k,n,flag,wei;
   double xs;
   int class_size[2]={0};
   double  Pi[2];
   int *Label;
   double **histvalue;
   double **cov; // for calculation of covariance matrix;
   double y, *Y; // result of FDA transfer

   double Sb[MaxSize][MaxSize], temp[MaxSize][MaxSize], Jwb[MaxSize][MaxSize], Sw[MaxSize][MaxSize];  // Sw is used for calculation of Sw and inverse Sw;
   double D[2*MaxSize][MaxSize], *S2;  //for calculation of SVD decomposition
   double d; 
   double *means, *deviations, *mean_class1, *mean_class2, *Mo;
   char buffer[300];
   char filename[50];
   char *ptr;
   FILE *fp;
   
   int iNagents = 10;
   int iitermax = 400;
   int indim = 22;
   double imaxx = 10;
   double imaxv = 0.5;
   double ierror = 1.0000e-04;
   double iw = 0.9;
   int iNoBg, iNoEn;
   int iNoOfFeatures = 1;
   int iPpW = 10;
   int iTimeCopies[2]={1,1};
   int iDelays[2]={9,0};
   int iNNArch[3]={10,2,1};
   int iM = 10;
   double *fopt, *bestpart, *Eevolution;
   int no_feature =0;

   double Data[MaxTrainingSample][WindowSize];
   double TrainData[MaxTrainingSample][WindowSize], TestData[MaxTrainingSample][WindowSize];
   int subBg, subEn, startBg, startEn, indexTrainBg, indexTestBg, indexTrainEn, indexTestEn;
   double output[MaxTrainingSample];
   int *randvBg, *randvEn;
   int endsBg[Fold], endsEn[Fold];
  
   double awindow[WindowSize];
   double errors[Fold];
   double particalweights[Fold][indim];

   //double *SS;
   for(i=0;i<MaxSize;i++)
      for(j=0;j<MaxSize;j++)
         Sw[i][j] = 0;
   for(i=0;i<MaxSize;i++)
      for(j=0;j<MaxSize;j++)
         Sb[i][j] = 0;
   for(i=0;i<MaxSize;i++)
      for(j=0;j<MaxSize;j++)
         temp[i][j] = 0;
   for(i=0;i<MaxSize;i++)
      for(j=0;j<MaxSize;j++)
         Jwb[i][j] = 0;
   for(i=0;i<2*MaxSize;i++)
      for(j=0;j<MaxSize;j++)
         D[i][j] = 0;  // initialize Sw, Sb, temp, Jwb and D
   for(i=0;i<Fold;i++)
      errors[i] = 0;

/* Step 1:======================================= Read training data from file ======================================*/
   printf("PLEASE INPUT THE FILE NAME FOR TRAINING DATA: ");
   scanf("%s",filename);
   //strcpy(filename,"training.txt");
   fp = fopen(filename,"r");
   if(fp == NULL) {perror(filename);}
   while(!feof(fp) && (fgets(buffer,300,fp)!= NULL)){
      if( (ptr = strstr(buffer, "Label") ) != NULL){ 
         // handle the header of the training file to find the dimension the data
         while(*ptr != '\0'){
            if((*ptr == '\t') && (*(ptr+1) != '\0') &&(*(ptr+1) != '\n') ){ no_column++;}
            ptr++;
         } 
      }
      else{ no_row++;
      }
   }
   //printf("There are %d columns and %d rows in the training data. \n", no_column, no_row);
   fclose(fp);
   if(no_column > MaxSize){printf("\n The code is currently programmed to process data with no more than %d histones.\n Please change the value of MaxSize in file CSIANNTrain.c then try again!\n", MaxSize/2); return;}
   if( (no_row/WindowSize) > MaxTrainingSample){printf("\n The code is currently programmed to process no more than %d training samples.\n Please change the value of MaxTrainingSample then try again!\n", MaxTrainingSample); return;}
 
   Label = malloc(no_row * sizeof(int));
   if(Label == NULL){printf("Allocate %d bytes for Label failed! \n", no_row * sizeof(int)); return;}
   histvalue = malloc(no_row * sizeof(double *));
   if(histvalue == NULL){printf("Allocate %d bytes for histvalue failed! \n", no_row * sizeof(double *)); return;}
   for(i=0;i<no_row;i++){
      histvalue[i] = malloc(no_column * sizeof(double));
      if(histvalue[i]==NULL){ printf("Allocate %d bytes for histvalue[%d] failed! \n", no_column * sizeof(double), i); return;}
   }
   for(i=0;i<no_row;i++)
      for(j=0;j<no_column;j++)
         histvalue[i][j] = 0;

   fp = fopen(filename,"r");
   if(fp == NULL) {perror(filename);}
   i =0;
   while(!feof(fp) && (fgets(buffer,300,fp)!= NULL)){
      if( (ptr = strstr(buffer, "Label") ) != NULL){ 
        // i=0;   // handle the header of the training file to find the dimension the data
      }
      else{
         Label[i] = (int) buffer[0] -48;
         if(Label[i] == 1){class_size[0]++;}
         else if(Label[i] == 2) {class_size[1]++;}
         else{printf("There are more than 2 classes! \nThe ANN model only works for 2 classes.\n"); return;}
         ptr = buffer;
         ptr++; ptr++; //ptr pointing to the first number
         k =0; flag =1; 
         while(flag!= -1){
            if(*ptr == '\0' || *ptr == '\n'){ flag = -1;} //'\0' and '\n'
            else if(*ptr=='.'){
               flag=0;
               xs = 0;
               wei = 0;
               ptr++;
            } 
            else if(*ptr == '\t'){
               k++;
               flag = 1;
               ptr++;
            }
            else { 
               if(flag==1) {          
                  histvalue[i][k]= histvalue[i][k]*10 + *ptr-48;
                  ptr++;
               }
               if(flag==0){
                  wei++;
                  xs=*ptr-48;
                  for(j=0;j<wei;j++)
                     xs = xs/10;
                  histvalue[i][k]= histvalue[i][k] + xs;
                  ptr++;
               }
            }
         }
         i++;
      }
   }
   //printf("Finish reading training data of %d rows.\nThe two class have size %d and %d respectively!\n", i, class_size[0], class_size[1]);
   fclose(fp); //*/
   iNoBg = class_size[0] / WindowSize;
   iNoEn = class_size[1] / WindowSize;

/* Step 2:======================================= Normalized the columns ======================================*/
// normalized the columns so that the columns have zero mean and 1 standard deviation 
// d[i][j] = (d[i][j] - mean[j])/std[j]
   means = malloc(no_column * sizeof(double));
   if(means == NULL) {printf("Allocate %d bytes for means failed!", no_column*sizeof(double)); return;}
   for(i=0;i<no_column;i++) { means[i] = 0;}

   deviations = malloc(no_column * sizeof(double));
   if(deviations == NULL) {printf("Allocate %d bytes for deviations failed!", no_column*sizeof(double)); return;}
   for(i=0;i<no_column;i++) { deviations[i] = 0;}
 
   for(i=0;i<no_column;i++){
      for(j=0;j<no_row;j++){ 
         means[i] = means[i]+histvalue[j][i];
      }
   }
   
   for(i=0;i<no_column;i++){
      means[i] = means[i] / no_row;
      //printf("%f\t", means[i]);
   }
   
   for(i=0;i<no_column;i++){
      for(j=0;j<no_row;j++){ 
         deviations[i] = deviations[i] + (histvalue[j][i] - means[i]) * (histvalue[j][i] - means[i]);
      }
   } 
   for(i=0;i<no_column;i++){
      deviations[i] = sqrt(deviations[i] / (no_row-1));
      //printf("%f\t",deviations[i]);
   }
   /*fp = fopen("NormalizedTraining.txt","w");
   if(fp == NULL) {perror("NormalizedTraining.txt");return;}*/
   for(j=0;j<no_row;j++){
      for(i=0;i<no_column;i++){ 
         histvalue[j][i] = (histvalue[j][i] - means[i] )/ deviations[i];
         //fprintf(fp, "%lf\t",histvalue[j][i]);
      }
     // fprintf(fp, "\n");
   }
   //fclose(fp);*/
   //free(means);  
   //free(deviations);
   /*printf("Writing Normalized Data to File NormalizedTraining.txt finished! \n");
   printf("Normalization finished!\n");*/

/* Step 3:=======================================     FDA transfer      ======================================*/
   Pi[0] = (double) class_size[0] / (double) no_row;
   Pi[1] = (double) class_size[1] / (double) no_row;
   //printf("%d\t%d\t%d\t%lf\t%lf\n",class_size[0],class_size[1], no_row, Pi[0], Pi[1]);
   //------calculate intraclass covariance matrix SW-------------------------//
 
   cov = malloc(no_row * sizeof(double *));   //start initialize cov
   if(cov == NULL){printf("Allocate %d bytes for cov failed! \n", no_row * sizeof(double *)); return;}
   for(i=0;i<no_row;i++){
      cov[i] = malloc(no_column * sizeof(double));
      if(cov==NULL){ printf("Allocate %d bytes for cov[%d] failed! \n", no_column * sizeof(double), i); return;}
   }   
   mean_class1 = malloc(no_column * sizeof(double));
   if(mean_class1 == NULL) {printf("Allocate %d bytes for mean_class1 failed!", no_column*sizeof(double)); return;}
   for(j=0;j<no_column;j++) { mean_class1[j] = 0;}
   for(i=0;i<class_size[0];i++){
      for(j=0;j<no_column;j++){ 
         mean_class1[j] = mean_class1[j]+histvalue[i][j];
      }
   }
   for(j=0;j<no_column;j++){
      mean_class1[j] = mean_class1[j] / class_size[0];   // mean for the first class
      //printf("%f\t", means[i]);
   }
   for(i=0;i<class_size[0];i++){
      for(j=0;j<no_column;j++){
         cov[i][j] = histvalue[i][j] - mean_class1[j];
      }
   }

   mean_class2 = malloc(no_column * sizeof(double));
   if(mean_class2 == NULL) {printf("Allocate % bytes for mean_class2 failed!", no_column*sizeof(double)); return;}
   for(j=0;j<no_column;j++) { mean_class2[j] = 0;}
   for(i=class_size[0];i<no_row;i++){
      for(j=0;j<no_column;j++){ 
         mean_class2[j] = mean_class2[j]+histvalue[i][j];
      }
   }
   for(j=0;j<no_column;j++){
      mean_class2[j] = mean_class2[j] / class_size[1];  // mean for the second class
      //printf("%f\t", means[i]);
   }
   for(i=class_size[0];i<no_row;i++){
      for(j=0;j<no_column;j++){
         cov[i][j] = histvalue[i][j] - mean_class2[j];
      }
   }                // end of cov

   for(i=0;i<no_column;i++){
      for(j=i;j<no_column;j++){
         for(k=0;k<class_size[0];k++){
            Sw[i][j] = Sw[i][j] + cov[k][i]*cov[k][j];
         }
         Sw[i][j] = Sw[i][j] / (class_size[0] -1);
         Sw[j][i] = Sw[i][j];
      }
   }
   for(i=0;i<no_column;i++){
      for(j=i;j<no_column;j++){
         Sw[i][j] = Sw[i][j] / Pi[0];
         Sw[j][i] = Sw[i][j];
      }
   }   // Sw = 1/Pi(1)*cov(X(1:dim(1),:));

   for(i=0;i<no_column;i++){
      for(j=i;j<no_column;j++){
         for(k=class_size[0];k<no_row;k++){
            temp[i][j] = temp[i][j] + cov[k][i]*cov[k][j];
         }
         temp[i][j] = temp[i][j] / (class_size[1] -1);
         temp[j][i] = temp[i][j];
      }
   }
   for(i=0;i<no_column;i++){
      for(j=i;j<no_column;j++){
         temp[i][j] = temp[i][j] / Pi[1];
         temp[j][i] = temp[i][j];
      }
   }  // temp = 1/Pi[2] * cov(X(sum(dim(1:i-1) +1: sum(dim(1:i)),:))

   //printf("\n Sw is :\n");
   for(i=0;i<no_column;i++){
      for(j=0;j<no_column;j++){
         Sw[i][j] = Sw[i][j] + temp[i][j];
         //printf("%lf\t", Sw[i][j]);
      }
      //printf("\n");
   } //  checked and right! */
   
   //-----------------------Classes Mean Vector-----------------------------------------------//
   Mo = malloc(no_column * sizeof(double));
   if(Mo == NULL) {printf("Allocate %d bytes for Mo failed!", no_column*sizeof(double)); return;}
   for(j=0;j<no_column;j++) { Mo[j] = 0;}
   for(j=0;j<no_column;j++){
      Mo[j] = mean_class1[j] * Pi[0] + mean_class2[j] * Pi[1]; 
      //printf("%f\t",Mo[j]);
   }  // checked and right!

   //-----------------------calculate interclass covariance matrix Sb-------------------------//
   for(i=0;i<no_column;i++) {mean_class1[i] = mean_class1[i] - Mo[i];}
   for(i=0;i<no_column;i++) {mean_class2[i] = mean_class2[i] - Mo[i];}
   for(i=0;i<no_column;i++){
      for(j=i;j<no_column;j++){
         Sb[i][j] = Pi[0] * mean_class1[i] * mean_class1[j]  + Pi[1] * mean_class2[i] * mean_class2[j]; 
         Sb[j][i]=Sb[i][j];
      }
   }  //checked and right!
   
   /*printf("\nSb is: \n");
   for(i=0;i<no_column;i++){
      for(j=0;j<no_column;j++){
         printf("%lf\t", Sb[i][j]);
      }
      printf("\n");
   } */ 

   //free(temp);
   //free(mean_class1);
   //free(mean_class2);
   //free(Mo);//*/ //when the pointers are freeed the value of Sw changed, wried!! 

   //-------------------------inverse(matrix) Sw------------------------
   inv((double *) Sw, no_column, MaxSize);
   /*printf("\nThe inverse of Sw is: \n");
   for(i=0;i<no_column;i++){
      for(j=0;j<no_column;j++){
         printf("%lf\t", Sw[i][j]);
      }
      printf("\n");
   }  //inverse(Sw) checked and right!*/

   for(i=0;i<no_column;i++)
      for(j=0;j<no_column;j++)
         for(k=0;k<no_column;k++)
            Jwb[i][j] = Jwb[i][j]+ Sw[i][k] * Sb[j][k];  // Jwb=Sw\Sb  checked and right!

   /*printf("\nJwb is: \n");
   for(i=0;i<no_column;i++){
      for(j=0;j<no_column;j++){
         printf("%lf\t", Jwb[i][j]);
      }
      printf("\n");
   } */ 
   //svd
   S2 = malloc(no_column * sizeof(double));
   for (i=0;i<no_column;i++)
      S2[i] = 0;

   for(i=0;i<no_column;i++)
      for(j=0;j<no_column;j++)
         D[i][j] = Jwb[i][j];
   for(i=no_column;i<2*no_column;i++)
      for(j=0;j<no_column;j++)
         D[i][j] = 0;
   svd2((double *)D, S2, no_column, MaxSize); //checked and right!
   /*printf("The result of SVD(D) is \n");
   for(i=0;i<2*no_column;i++){
      for(j=0;j<no_column;j++)
         printf("%lf\t", D[i][j]);
      printf("\n");
   }*/

   Y = malloc( no_row  * sizeof(double));
   if(Y == NULL) {
      printf("Allocate space for Y failed!\n");
      return;
   } //*/
   //fp = fopen("FDAresultTraining.txt","w");
   //if(fp == NULL) {perror("FDAresultTraining.txt");return;}
   for(i=0;i<no_row;i++){
      y = 0;
      for (j=0;j<no_column;j++)
         y = y + histvalue[i][j] * D[no_column+j][0];
      Y[i] = y;
     // fprintf(fp,"%lf\n",Y[i]);
   } 
   //fclose(fp);
   //printf("FDA transfer finished!\n");
   fp = fopen("Features.txt","w");
   if(fp == NULL) perror("Features.txt");
   for(i=0;i<no_column;i++)
      fprintf(fp,"%lf\n",D[no_column+i][0]);
    fclose(fp);
   printf("\nFDA result is saved in file Features.txt\n");
/* Step 4:======================================= train and evaluate TDNN models ======================================*/
   fopt = malloc(sizeof(double));
   fopt[0] = 0;
   /*bestpart = malloc(indim * sizeof(double)); 
   for(i=1;i<indim;i++)
      bestpart[i] = 0;*/  
   Eevolution = malloc(iitermax * sizeof(double)); 
   for(i=1;i<iitermax;i++)
      Eevolution[i] = 0;

   i=0; j=0;
   for (k=0;k<no_row;k++){
      if(i< WindowSize){
         Data[j][i] = Y[k];
         i++;
      }
      else{
         j++;
         i = 0;
         Data[j][i] = Y[k];
         i++;
      }
   }
   // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  Cross-one Valiation  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   randvBg = malloc(iNoBg * sizeof(int));
   randperm(randvBg, iNoBg);  //randomly permuted  index of training dataset  [0~ iNoBg-1]
   /*for(j=0;j<iNoBg;j++) printf("%d\t", randvBg[j]);
   PRINTN*/

   randvEn = malloc(iNoEn * sizeof(int));
   randperm(randvEn, iNoEn); //randomly permuted  index of testing dataset  [0 ~ iNoEn-1]
   /*for(j=0;j<iNoEn;j++) printf("%d\t", randvEn[j]);
   PRINTN*/

   for(i=0;i<Fold;i++){
      endsBg[i] = (int) ((double) iNoBg / (double) Fold * (double)(i+1));   //index end of each fold of training data, endsBg[Fold-1] = iNoBg 
      endsEn[i] = (int) ((double) iNoEn / (double) Fold * (double)(i+1));   //index end of each fold of testing data, endsEn[Fold -1] = iNoEn
   }
      

   if(Fold == 1){
      i=0;
      OneLayersTimeDelayNeuralNet(fopt, (double *)particalweights[i], Eevolution, iNagents, iitermax, indim, imaxx, imaxv, ierror, iw, iNoBg, iNoEn, iNoOfFeatures, iPpW, iTimeCopies, iDelays, iNNArch, (double *)Data, iM);

      // Evaluation the TDNN using the TestData
      //// using the output to test all the cutoff and construct a ROC curve, calculate the AUC-area and find the best cutoff of current run
      for(j=0;j<(iNoBg + iNoEn);j++){
         for(k=0;k<WindowSize;k++)
            awindow[k] = Data[j][k]; 
         output[j] = TDNNEva((double *)awindow, WindowSize, (double *) particalweights[i]);
      }
      for(j=0;j<(iNoBg + iNoEn);j++){
         if(j<iNoBg){
            errors[i] = errors[i] + output[j] * output[j];
         }
         else{
            errors[i] = errors[i] + (1-output[j]) * (1-output[j]);
         }
      }
      errors[i] = errors[i]/(iNoBg+iNoEn);

      //printf("\nThe Training for Fold %d finished!\n Mean Squared Error on testing data is %lf \n",i+1, errors[i]);
   }
   else{
      printf("\nStart training the ANN model...\n");
      for(i=0;i<Fold;i++){
         if(i==0){
            subBg = endsBg[0];
            subEn = endsEn[0];
            startBg = 0;
            startEn = 0;
           // printf("\n In Fold %d, subBg =%d, subEn =%d, startBg = %d, startEn =%d\n",i+1,subBg,subEn,startBg,startEn);
         }
         else{
            subBg = endsBg[i] - endsBg[i-1];  // number of Bg samples in the training subset 
            subEn = endsEn[i] - endsEn[i-1];  // number of En samples in the training subset
            startBg = endsBg[i-1];
            startEn = endsEn[i-1];
            //printf("\n In Fold %d, subBg =%d, subEn =%d, startBg = %d, startEn =%d\n",i+1,subBg,subEn,startBg,startEn);
         }
   // put the data into TrainData and TestData respectively
         indexTrainBg = 0; indexTestBg =0;  // put the backgroup into the Train and Test subset
         for(j=0; j<iNoBg;j++)
         {  
            if((j>=startBg) && (j< (startBg + subBg)))
            {  
               for(k=0;k<WindowSize;k++)
                  TestData[indexTestBg][k]= Data[randvBg[indexTestBg]][k];
               indexTestBg++;
            }  
            else
            {
               for(k=0;k<WindowSize;k++)
                  TrainData[indexTrainBg][k]= Data[randvBg[indexTrainBg]][k];
               indexTrainBg++;
            }
         }
        // printf("\nCheck: subBg = %d vs indexTestBg=%d; iNoBg - subBg = %d vs indexTrainBg = %d\n", subBg, indexTestBg, iNoBg-subBg, indexTrainBg);

         indexTrainEn = 0; indexTestEn =0; // put the Enhancer into the Train and Test subset
         for(j=0; j<iNoEn;j++)
         {  
            if((j>=startEn) && (j< (startEn + subEn)))
            {  
               for(k=0;k<WindowSize;k++)
                  TestData[indexTestBg + indexTestEn][k]= Data[iNoBg + randvEn[indexTestEn]][k];
               indexTestEn++;
            }  
            else
            {
               for(k=0;k<WindowSize;k++)
                  TrainData[indexTrainBg + indexTrainEn][k]= Data[iNoBg + randvEn[indexTrainEn]][k];
               indexTrainEn++;
            }
         }
         //printf("\nCheck: subEn = %d vs indexTestEn=%d; iNoEn - subEn = %d vs indexTrainEn = %d\n", subEn, indexTestEn, iNoEn-subEn, indexTrainEn);
      
         // training the model using the TrainData 
         OneLayersTimeDelayNeuralNet(fopt, (double *)particalweights[i], Eevolution, iNagents, iitermax, indim, imaxx, imaxv, ierror, iw, indexTrainBg, indexTrainEn, iNoOfFeatures, iPpW, iTimeCopies, iDelays, iNNArch, (double *)TrainData, iM);

         // Evaluation the TDNN using the TestData
         //// using the output to test all the cutoff and construct a ROC curve, calculate the AUC-area and find the best cutoff of current run
         for(j=0;j<(subBg + subEn);j++){
            for(k=0;k<WindowSize;k++)
               awindow[k] = TestData[j][k]; 
            output[j] = TDNNEva((double *)awindow, WindowSize, (double *) particalweights[i]);
         }
         for(j=0;j<(subBg + subEn);j++){
            if(j<subBg){
               errors[i] = errors[i] + output[j] * output[j];
            }
            else{
               errors[i] = errors[i] + (1-output[j]) * (1-output[j]);
            }
         }
         errors[i] = errors[i]/(subBg+subEn);

         //printf("\n%dth out of %d training finished!",i+1, Fold);
         //printf("\nMean Squared Error on testing data is %lf \n",errors[i]);
      }
   }

  // ----------------------------------------------------------------------------------------------------------------------------------//
  // found the smallest error:
   xs = errors[0]; //
   k = 0;
   for(i=1;i<Fold;i++){
      if(errors[i] < xs){
         xs = errors[i];
         k = i;
      }
   }
   printf("\nSetting of the best trained ANN model is saved in file partical_weights.txt\n");
  // same the bestparticals
   fp = fopen("partical_weights.txt","w");
   if(fp == NULL) {perror("partical_weights.txt"); return;}
   for(j=0;j<indim;j++){
      fprintf(fp,"%lf\n", particalweights[k][j]);
   }
   fclose(fp);
  // when all the Fold are run, choose the smalles AUC-area TDNN model and choose the corresponding cutoff, save the bestpartical and cutoff
}



