// gcc-o c CSIANNPrediction.c -lm
/*input:
testing.txt
Features.txt  // file for the FDA parameters
partical_weights.txt  // file for the TDNN parameters
  output:
Prediction.txt   the prediction result
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>

#define MaxHistNo 25 // permitted maximum number of histones would be MaxHistNo = 25
#define MaxSize MaxHistNo*2 
#define WindowSize 10 //size of the window is 10, the window for enhancer prediction is 2kb, so 2000/Resolution = 10
#define PeakWindowSize 2500 
#define Resolution 200

#define TimeofCutting 100  // TimeofCutting should be smaller than 39, otherwise you'll have to add some code
#define Center2Center 2500
#define OUTLIE 3 // the window would be predicted positive only when there're more than 3 items have value>threshold. If OUTLIE =-1 there's no filtering
#define FDRCUTOFF 0.05 // 
//#define FootStep 0.02
#define FootStep 0.01

#define MaxEnhancerNo  100000
typedef struct {
   char chrID[10];
   int positionID;
   float histvalues[MaxSize/2];
} struct_position;


#include "TDNNEva.h"
#include "randperm.h"
#include "randperm2.h"
#include "RandomPerPrediction3.h"


int main(int argc, char *argv[]){

   int i,j,k,flag,q;
   int no_column=0, no_row=0;
   
   int PeakWindow = (int) PeakWindowSize / Resolution;
   int offset;   

   double cutoff[TimeofCutting]={0};
   int RandomFP[TimeofCutting]={0};
   int predictedP[TimeofCutting];
   double FDR[TimeofCutting];

   char filename[200], fileout[200];
   char buffer[300];
   double weights[22];
   
   struct_position *Position;

   FILE *fp;
   char *ptr;

   float xs, wei;
   float *means, *deviations;   
   
   double *features, *fdaresult, *output, awindow[WindowSize];
   int *outputflag; // a flag to show whether a window is OK to be a prediction

   double space;
   int predictions[MaxEnhancerNo];
   double maxx,randmaxx;
   double minn,randminn;
   double fdamean;

   float *positiontemp; // for scheme 1 of generating randomlized data

   double *randomfda, *randomoutput; // for randomlize the data
   int *randomoutputflag; // a flag to show whether a window is OK to be a prediction

   int x, randrow, randcolumn, *randindex;

   int  position; 
   float curdistance, mindistance;
   int FDRindex;

// get input : testing file
    if(argc != 3) {
           printf("Input wrong! ./XXX testing.txt   output\n");
	   exit(1);
    }	
   strcpy(filename, argv[1] );  
   strcpy(fileout, argv[2] );  
   printf("reading testing file: %s, footstep = %e, no.cutting = %d \n", filename, FootStep, TimeofCutting);

/* Step 1:======================================= Read testing data from file ======================================*/
//   printf("PLEASE INPUT THE FILE NAME FOR TESTING DATA: ");
//   scanf("%s",filename);
   fp = fopen(filename,"r");
   if(fp == NULL) {perror(filename);}
   while(!feof(fp) && (fgets(buffer,300,fp)!= NULL)){
      if( (ptr = strstr(buffer, "ID") ) != NULL){ 
         // handle the header of the training file to find the dimension the data
         while(*ptr != '\0'){
            if((*ptr == '\t') && (*(ptr+1) != '\0') &&(*(ptr+1) != '\n') ){ no_column++;}
            ptr++;
         } 
      }
      else{ no_row++;
      }
      bzero(buffer,300);
   }
   fclose(fp);
   no_column--; // the number of histions value columns

   //printf("There are %d columns and %d rows in the Testing data. \n", no_column, no_row);
   if(no_column -1> MaxSize){printf("\n The code is currently programmed to process data with no more than %d histones. Please change the value of MaxSize then try again!\n", MaxSize/2); return;}

   Position = malloc(no_row * sizeof(struct_position));
   space = no_row * sizeof(struct_position);
   if(Position == NULL){printf("Allocate %lf MB for Position failed! \n", space/(1024*1024)); return;}
   //else{printf("Allocate %lf MB for Position succeeded! \n", space /(1024 * 1024));}
   else{printf("\nAllocate memory for input data succeeded!\nReading file %s...\n",filename);}

   for(i=0;i<no_row;i++){
      bzero(Position[i].chrID,10);
      Position[i].positionID = 0;
      for(j=0;j<no_column;j++)
         Position[i].histvalues[j]=0;
   }

   fp = fopen(filename,"r");
   if(fp == NULL) {perror(filename);}
   i =0; bzero(buffer,300);
   while(!feof(fp) && (fgets(buffer,300,fp)!= NULL) ){
      if( (ptr = strstr(buffer, "ID") ) != NULL){  // i=0;   // handle the header of the training file to find the dimension the data
      }
      else{
         ptr = (char *)buffer;
         j=0;
         while(*ptr != '\t'){
            Position[i].chrID[j] = *ptr;
            ptr++;
            j++;
         }
         Position[i].chrID[j] = '\0';
         ptr++; //ptr pointing to positions now
         while(*ptr != '\t'){
            Position[i].positionID = Position[i].positionID *10 + *ptr -48;
            ptr++;
         }
         ptr++; //ptr pointing to the histone values now

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
                     Position[i].histvalues[k]= Position[i].histvalues[k]*10 + *ptr-48;
                     ptr++;
                  }
                  if(flag==0){
                     wei++;
                     xs=*ptr-48;
                     for(j=0;j<wei;j++)
                        xs = xs/10;
                     Position[i].histvalues[k]= Position[i].histvalues[k] + xs;
                     ptr++;
                  }
               }
         }
         bzero(buffer,300);
         i++;
      }
   }
   fclose(fp); //*/   

/* Step 2:======================================= Normalized the columns ======================================*/
// normalized the columns so that the columns have zero mean and 1 standard deviation 
// d[i][j] = (d[i][j] - mean[j])/std[j]
   means = malloc(no_column * sizeof(float));
   if(means == NULL) {printf("Allocate %d bytes for means failed!", no_column*sizeof(float)); return;}
   for(i=0;i<no_column;i++) { means[i] = 0;}
   deviations = malloc(no_column * sizeof(float));
   if(deviations == NULL) {printf("Allocate %d bytes for deviations failed!", no_column*sizeof(float)); return;}
   for(i=0;i<no_column;i++) { deviations[i] = 0;}
 
   for(i=0;i<no_column;i++){
      for(j=0;j<no_row;j++){ 
         means[i] = means[i]+Position[j].histvalues[i];
      }
   }
   for(i=0;i<no_column;i++){
      means[i] = means[i] / no_row;
   }

   for(i=0;i<no_column;i++){
      for(j=0;j<no_row;j++){ 
         deviations[i] = deviations[i] + (Position[j].histvalues[i] - means[i]) * (Position[j].histvalues[i] - means[i]);
      }
   } 
   for(i=0;i<no_column;i++){
      deviations[i] = sqrt(deviations[i] / (no_row-1));
   }
   for(j=0;j<no_row;j++){
      for(i=0;i<no_column;i++){ 
         Position[j].histvalues[i] = (Position[j].histvalues[i] - means[i] )/ deviations[i];
      }
   }
   free(deviations);

/* Step 3:=======================================     FDA transfer      ======================================*/
   fp = fopen("Features.txt","r");
   if(fp == NULL) {perror("Features.txt");return;}
   features = malloc(no_column*sizeof(double));
   if(features == NULL){printf("Allocate %d bytes for features failed!\n", no_column*sizeof(float)); return;}
   for(i=0;i<no_column;i++)
      fscanf(fp,"%lf",&features[i]);
   fclose(fp);

   fdaresult = malloc(no_row * sizeof(double));
   if(fdaresult == NULL){printf("Allocate %d bytes for features failed!\n", no_row*sizeof(float)); return;}
   for(i=0;i<no_row;i++)
      fdaresult[i] = 0;
   for(i=0;i<no_row;i++){
      for (j=0;j<no_column;j++)
         fdaresult[i] = fdaresult[i] + (double) Position[i].histvalues[j] * features[j];
   } 

   printf("\nTesting Data FDA transfering finished!\n");
  
/* Step 4:============================     Calculate output of ANN model on testing data======================================*/
   
   /*4-1 ===============================    Reading weight setting for the ANN model ===============*/
   fp = fopen("partical_weights.txt","r");
   if(fp == NULL) {perror("partical_weights.txt");return;}
   for(i=0;i<22;i++)
      fscanf(fp,"%lf",&weights[i]);
   fclose(fp);
   
   /*4-2===========  calculate the output of ANN model of each 2kb windows with delay 1 and save the output in file  ========*/
   output = malloc( (no_row - WindowSize + 1) * sizeof (double));
   if(output == NULL) {printf("Allocate space for output failed!\n"); return;}
   for(i=0;i< no_row - WindowSize + 1;i++)
      output[i] = 0;
   outputflag = malloc((no_row - WindowSize +1) * sizeof(int));
   if(outputflag == NULL) {printf("Allocate space for outputflag failed!\n");return;}
   for(i=0;i<no_row - WindowSize +1;i++)
      outputflag[i] = 1;
   fdamean=0;
   for(i=0;i<no_row;i++)
      fdamean = fdamean + fdaresult[i];
   fdamean = fdamean/no_row;

   maxx = 0;
   minn = 1;
   for(i=0;i< no_row - WindowSize + 1;i++){
      k=0;
      for(j=0;j<WindowSize;j++){
         awindow[j] = fdaresult[i+j];
         if(awindow[j] > fdamean) k++;
      }
      output[i] = TDNNEva((double *)awindow, WindowSize, (double *) weights);
      if(k<= OUTLIE)
         outputflag[i] =0; 
      if(output[i] > maxx) maxx = output[i];
      if(output[i]< minn) minn = output[i];
   }

   printf("\nWorking with the ANN model...\n");

   fp=fopen("ANNoutput","w");
   if(fp == NULL){perror("ANNoutput");}
   offset = WindowSize/2;
   for (i=0 ; i < no_row - WindowSize + 1; i++){
       fprintf(fp,"%s\t%d\t%f\n", Position[i + offset].chrID, Position[i + offset].positionID, output[ i ]);
   }
   fclose(fp);

/* Step 5: ============================  Randomlize the data and calculate ANN output of randomlized data =====================*/
   randomfda = malloc(no_row * sizeof(double));
   if(randomfda == NULL) {printf("\nAllocate space for randomfda failed!"); return;}
   for(i=0;i<no_row;i++) 
      randomfda[i]=0;    

  // generating randomlized data Scheme 2
   positiontemp = malloc (no_row * sizeof (float));
   for(i = 0;i<no_column;i++){
      for(j=0;j<no_row;j++)
         positiontemp[j] = Position[j].histvalues[i];

      randperm2(positiontemp, no_row );

      for(j=0;j<no_row;j++)
         Position[j].histvalues[i] = positiontemp[j];     
      
   }
   for(i=0;i<no_row;i++){  // random permuate all of the items in Position
      for(j=0;j<no_column;j++){ 
         randomfda[i] = randomfda[i] + (double) Position[i].histvalues[j] * features[j]; 
      }
   }  
//

/*  // generating randomlized data Scheme 2
   randindex = malloc( no_row * no_column * sizeof(int));
   if(randindex == NULL) {printf("\nAllocate space for randindex failed!\n"); return;}
   randperm(randindex, no_row * no_column);  // generate a random index for all of the items in Position
      
   //printf("\nStarting permutation!");
   for(i=0;i<no_row;i++){  // random permuate all of the items in Position
      for(j=0;j<no_column;j++){ 
         x = randindex[i*no_column + j];
         randrow = (int) x / no_column;  
         randcolumn = x % no_column;
         randomfda[i] = randomfda[i] + (double) Position[randrow].histvalues[randcolumn] * features[j]; 
      }
   }
   //printf("\nRandomly permutation finished !"); 

*/

   randomoutput = malloc( (no_row - WindowSize + 1) * sizeof (double));
   if(randomoutput == NULL) {printf("Allocate space for randomoutput failed!\n");return;}
   for(i=0;i< no_row - WindowSize + 1;i++)
      randomoutput[i] = 0;
   randomoutputflag = malloc((no_row - WindowSize +1) * sizeof(int));
   if(randomoutputflag == NULL) {printf("Allocate space for randomoutputflag failed!\n");return;}
   for(i=0;i<no_row - WindowSize +1;i++)
      randomoutputflag[i] = 1;

   randmaxx = 0;
   randminn = 1;
   for(i=0;i< no_row - WindowSize + 1;i++){
      k=0;
      for(j=0;j<WindowSize;j++){
         awindow[j] = randomfda[i+j];
         if(awindow[j] > fdamean) k++;  
      }
      randomoutput[i] = TDNNEva((double *)awindow, WindowSize, (double *) weights);
      if(k<=OUTLIE)
         randomoutputflag[i] = 0;
      if(randomoutput[i] > randmaxx) randmaxx = randomoutput[i];
      if(randomoutput[i]< randminn) randminn = randomoutput[i];
   }
   free(randomfda);

/* Step 6 ===========================  calculate prediction based on each possible cutoff =========================*/
    printf("maxx = %.5f, minn = %.5f, footstep = %.5f \n", maxx, minn, FootStep);
   // 6-1 ==== set possible cutoffs ===========================================================//
   offset = WindowSize/2;
   for(i=0;i<TimeofCutting;i++) {
//      cutoff[i] = maxx - 0.00125 * (i+1);
//      cutoff[i] = 0.5 - 0.025 * (i+1);
      cutoff[i] = maxx - FootStep * (i+1);
//      cutoff[i] = 0.5 - 0.01 * (i+1);
   }
//   cutoff[TimeofCutting -1 ] = 0.5;

   // 6-2 ==== run the prediction using each cutoff  ===========================================
   for(q=0;q<TimeofCutting;q++){
      //6-2-1 ========  find the peaks in each 25Kb windows ====================================
      j=0;
      if(no_row - WindowSize + 1 <= PeakWindow){  // when the searching space is no larger than one PeakWindow
         maxx = 0;
         for(i=0;i<no_row-WindowSize+1;i++){
            if(output[i]> cutoff[q] && outputflag[i]){
               if(output[i] > maxx){
                  maxx = output[i];
                  j=1;
               }
            }
         }
      }
      else{
         for (i=0;i<no_row - WindowSize +1; i++){
            if(output[i]> cutoff[q] && outputflag[i]) {
               flag =1;
               k = 1;
               if( i < PeakWindow / 2 ){  // the first a few positions
                  while(flag && (k<= PeakWindow/2)){
                     if(output[i] <= output[i+k]) flag =0;
                     k++;  
                  }
                  k=1;
                  while(flag && k<=i){
                     if(output[i] < output[i-k]) flag =0;
                     k++;
                  } 
               }
               else if(i > no_row-WindowSize+1 - PeakWindow/2){  // the last a few positions
                  while(flag && k<= PeakWindow/2){
                     if(output[i]<output[i-k]) flag =0;
                        k++;
                  }
                  k = no_row-WindowSize+1 -1;
                  while(flag && k>i){
                     if(output[i]<= output[k]) flag = 0;
                     k--;
                  }
               }
               else{      //We check [-PeakWindow/2 ~ +PeakWindow/2] positions (that's 12 positions in total when PeakWinwod = 12) to find a peak in the PeakWindow
                  while(flag && k<=PeakWindow/2){
                     if(output[i]<output[i-k] || output[i]<=output[i+k]) flag =0;
                     k++;
                  }
               }
               if(flag){  // if output[i] is the peak of its PeakWindow  
                  if(j < MaxEnhancerNo)
                     j++;   
                  i = i + PeakWindow/2; //move forward to the next position which is out of the current window, making the minimum distance bwteen two predictions to be 7*200=1400
               }
            }
         } 
      }
      predictedP[q] = j;

      //6-2.2 ==== calculate FDR ================================================================
      RandomFP[q] = RandomPerPrediction3(randomoutput, randomoutputflag, no_row, weights, cutoff[q]);
      FDR[q] = (double) RandomFP[q]/predictedP[q];
      printf("FDR[%d] = %.5f, at cutoff = %.5f , with %d predictions.\n", q + 1, FDR[q], cutoff[q], predictedP[q] );

   }
/* Step 7 ========  Choose a cutoff when FDR is closest to 0.05 =================================== */
   for(i=0;i<TimeofCutting;i++){
      if(i==0){
         mindistance = FDR[i] - FDRCUTOFF;
         if(mindistance <0) mindistance = - mindistance;
         FDRindex = i;
      }
      else {
         curdistance = FDR[i] - FDRCUTOFF;
         if(curdistance <0) curdistance = - curdistance;
         
         if(curdistance < mindistance){ 
            mindistance = curdistance;
            FDRindex = i;
         }
      }
   }
   printf("\nCutoff was chosen when FDR = %f, with %d predictions. \n",FDR[FDRindex], predictedP[FDRindex]);

//   FDRindex = 40;

   // ------------------------------Set the cutoff and Save the Predictions when FDR is closest to 0.05 -----------------//
      j=0;
      if(no_row - WindowSize + 1 <= PeakWindow){  // when the searching space is no larger than one PeakWindow
         maxx = 0;
         for(i=0;i<no_row-WindowSize+1;i++){
            if(output[i]> cutoff[FDRindex] && outputflag[i]){
               if(output[i] > maxx){
                  maxx = output[i];
                  predictions[0] = i;  // there's only one peak
                  j=1;
               }
            }
         }
      }
      else{
         for (i=0;i<no_row - WindowSize +1; i++){
            if(output[i]> cutoff[FDRindex] && outputflag[i]) {
               flag =1;
               k = 1;
               if( i < PeakWindow / 2 ){  // the first a few positions
                  while(flag && (k<= PeakWindow/2)){
                     if(output[i] <= output[i+k]) flag =0;
                     k++;  
                  }
                  k=1;
                  while(flag && k<=i){
                     if(output[i] < output[i-k]) flag =0;
                     k++;
                  } 
               }
               else if(i > no_row-WindowSize+1 - PeakWindow/2){  // the last a few positions
                  while(flag && k<= PeakWindow/2){
                     if(output[i]<output[i-k]) flag =0;
                        k++;
                  }
                  k = no_row-WindowSize+1 -1;
                  while(flag && k>i){
                     if(output[i]<= output[k]) flag = 0;
                     k--;
                  }
               }
               else{      //We check [-PeakWindow/2 ~ +PeakWindow/2] positions (that's 12 positions in total when PeakWinwod = 12) to find a peak in the PeakWindow
                  while(flag && k<=PeakWindow/2){
                     if(output[i]<output[i-k] || output[i]<=output[i+k]) flag =0;
                     k++;
                  }
               }
               if(flag){  // if output[i] is the peak of its PeakWindow  
                  predictions[j]=i;
                  if(j < MaxEnhancerNo)
                     j++;   
                  i = i + PeakWindow/2; //move forward to the next position which is out of the current window, making the minimum distance bwteen two predictions to be 7*200=1400
               }
            }
         } 
      }
      if (j== MaxEnhancerNo)
         printf("There are more than %d predictions. Some predictions are disrecarded. Please increase cutoff.\n", MaxEnhancerNo);

   j = predictedP[FDRindex];
   fp=fopen(fileout,"w");
   if(fp == NULL){perror(fileout);}
   for (i=0;i<j;i++){
       fprintf(fp,"%s\t%d\t%f\n", Position[predictions[i]+offset].chrID, Position[predictions[i]+offset].positionID, output[predictions[i]]);
//       fprintf(fp,"%s\t%d\t%f\n", Position[predictions[i]+offset].chrID, Position[predictions[i]+offset].positionID, output[predictions[i]]/maxx ); // normalize ANN output by maxx value
   }
   fclose(fp);

   printf("Prediction result is saving in %s!\n", fileout);

}
