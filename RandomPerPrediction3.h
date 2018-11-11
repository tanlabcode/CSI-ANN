/*
  Randomly Permuatate the FDAresult and do prediction based on given CSIANN model and given cutoff 
*/

// gcc-o c CSIANNPrediction.c -lm
/*input:
  output
the predicted positions for the randomly permutaed data 
*/


int RandomPerPrediction3(randomoutput, randomoutputflag, no_row, weights, cutoff)
double *randomoutput;
int *randomoutputflag;
int no_row;
double *weights;
double cutoff;
{ 

   int i,j,k,flag;
   
   int PeakWindow = (int) PeakWindowSize / Resolution;
   int offset;   

   double maxx;
   double minn;

   offset = WindowSize/2;

   j=0;
   if(no_row - WindowSize + 1 <= PeakWindow){  // when the searching space is no larger than one PeakWindow
      maxx = 0;
      for(i=0;i<no_row-WindowSize+1;i++){
         if(randomoutput[i]> cutoff && randomoutputflag[i]){
            if(randomoutput[i] > maxx){
               maxx = randomoutput[i];
               //predictions[0] = i;  // there's only one peak
               j=1;
            }
         }
      }
   }
   else{
      for (i=0;i<no_row - WindowSize +1; i++){
         if(randomoutput[i]> cutoff && randomoutputflag[i]) { 
            flag =1;
            k = 1;
            if( i < PeakWindow / 2 ){  // the first a few positions
               while(flag && (k<= PeakWindow/2)){
                  if(randomoutput[i] <= randomoutput[i+k]) flag =0;
                  k++;  
               }
               k=1;
               while(flag && k<=i){
                  if(randomoutput[i] < randomoutput[i-k]) flag =0;
                  k++;
               } 
            }
            else if(i > no_row-WindowSize+1 - PeakWindow/2){  // the last a few positions
               while(flag && k<= PeakWindow/2){
                  if(randomoutput[i]<randomoutput[i-k]) flag =0;
                     k++;
               }
               k = no_row-WindowSize+1 -1;
               while(flag && k>i){
                  if(randomoutput[i]<= randomoutput[k]) flag = 0;
                  k--;
               }
            }
            else{      //We check [-PeakWindow/2 ~ +PeakWindow/2] positions (that's 12 positions in total when PeakWinwod = 12) to find a peak in the PeakWindow
               while(flag && k<=PeakWindow/2){
                  if(randomoutput[i]<randomoutput[i-k] || randomoutput[i]<=randomoutput[i+k]) flag =0;
                  k++;
               }
            }
            if(flag){  // if output[i] is the peak of its PeakWindow  
               //predictions[j]=i;
               j++;   
               i = i + PeakWindow/2; //move forward to the next position which is out of the current window, making the minimum distance bwteen two predictions to be 7*200=1400
            }
         }
      } 
   }

   return(j);


}
