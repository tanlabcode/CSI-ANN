#include "randperm.h"

#define Window 2000
#define HalfWindowSize 1000
#define PromoterSize 2500
#define Resolution 200
#define WindowSize 10   //WindowSize = Window/Resolution
#define SignalWindow 10000 

typedef struct {
   char chrID[10];
   long center;
   long enhancer_start;
   long enhancer_end;
} struct_enhancer;

typedef struct {
   char chrID[10];
   long positionID;
   float histvalues[20];
   float histenergy[20];
   int  enhancer_flag; // identify whether this position belongs to an enhancer or not, 0 means background, 1 means in an enhancer Window, 2 means an enhancer center 
   int  promoter_flag; // identify whether this position belongs to a promoter or not, 0 means background, 1 means in a promoter region
} struct_position;

typedef struct {
   char chrID[10];
   char strand;
   long txStart;
   long txEnd;
} struct_gene;

void markenhancer( chroms, no_chr, positions, no_position, enhancers, no_enhancer)
long chroms[];
int no_chr;
struct_position *positions;
long no_position;
struct_enhancer *enhancers;
int no_enhancer;
{
   int i,flag,j, mindistance, e;
   long k,q, minindex;
   long enhancer_start,enhancer_end,enhancer_center;
   char chrID[10];
   long chrstart;

   for(i=0,q=0,e=0;i<no_enhancer;i++){
      enhancer_start = enhancers[i].enhancer_start;
      enhancer_end = enhancers[i].enhancer_end;
      enhancer_center = enhancers[i].center;
      strcpy(chrID,enhancers[i].chrID);

      flag = 1; // match has not been reached
      j=0;
      chrstart = 0;
      while(flag && (j<no_chr)){
         if(!strcmp(chrID,positions[chrstart].chrID)){

            flag = 0;// match found
            for(k=0;k<chroms[j];k++){ 
               if(k==0){
                  mindistance = abs(positions[chrstart+k].positionID - enhancer_center);
                  minindex = k;
               }
               else if( abs(positions[chrstart+k].positionID - enhancer_center) < mindistance){
                  mindistance = abs(positions[chrstart+k].positionID - enhancer_center);
                  minindex = k;
               }
               if( (positions[chrstart+k].positionID<=enhancer_end) && (positions[chrstart+k].positionID > enhancer_start) ){
                  if(positions[chrstart+k].enhancer_flag == 0){  // if the position hasn't been claimed to be an enhancer-center nor within an enhancer region
                     positions[chrstart+k].enhancer_flag = 1; //claim the position as in an enhancer region
                     e++;
                  }
               } 
            }
            q++;
            positions[chrstart+minindex].enhancer_flag = 2;   //mark the enhancer centers

         }
         else{

            chrstart = chrstart + chroms[j];
            j++;
         }
      }

   }
   //fclose(fp);
   printf("Enhancer regions have been marked.\n");

}

void markpromoter( chroms, no_chr, positions, no_position, genes, no_gene)
long chroms[];
int no_chr;
struct_position *positions;
long no_position;
struct_gene *genes;
int no_gene;
{
   int i,j,flag;
   long k, q=0, promoter_start, promoter_end, chrstart;
   char chrID[10];


   for(i=0,q=0;i<no_gene;i++){
      strcpy(chrID, genes[i].chrID);

      if(genes[i].strand == '+'){
         promoter_end = genes[i].txStart;
         promoter_start = promoter_end - PromoterSize;
      }
      else{
         promoter_start = genes[i].txEnd;
         promoter_end = promoter_start + PromoterSize;
      }

      flag = 1; // match has not been reached
      j=0;
      chrstart = 0;
      while(flag && (j < no_chr)){
         if(!strcmp(chrID, positions[chrstart].chrID)){
            flag = 0;// match found
            for(k=0;k<chroms[j];k++){ 
               if( (positions[chrstart+k].positionID<promoter_end) && (positions[chrstart+k].positionID > promoter_start) ){
               positions[chrstart+k].promoter_flag = 1; 
               q++;
               }
            }
         }
         else{

            chrstart = chrstart + chroms[j];
            j++;
         }
      }

   }

  printf("Promoter regions have been marked.\n");
}

void writetraining( positions, no_position, no_hist, no_enhancer)
struct_position *positions;
long no_position;
int no_hist;
int no_enhancer;
{
   int i,j,r,f,t=0,e,offset, sw;
   long k;
   char chrID[10];
   int flag = 1;
   int *randv, randrow, randcolumn;
   FILE *ftraining, *fp;

   ftraining = fopen("training.txt","w");
   if(ftraining == NULL){perror("training.txt");}
   //printf("\nStart writing file training.txt...\n");

   fprintf(ftraining, "Label");
   for(i=0;i<no_hist;i++){fprintf(ftraining,"\tHist%d\tHist%d_Energy",i+1, i+1);} 
   fprintf(ftraining,"\n"); //done writing the header
   
   randv = malloc(no_position * sizeof(int));
   if(randv == NULL) {printf("Allocate space for randv failed! \n"); abort();}
   randperm(randv, no_position); // randperm value [ 0~ no_position-1];
   i=0;
   for(e=0;e< 10 * no_enhancer;e++){  //picking out (10 * no_enhancer) non-enhancer, non-promoter positions as background for training
      r = randv[i];
      i++;
      if((positions[r].enhancer_flag == 0) && (positions[r].promoter_flag ==0) ){  //avoid picking ehancers or promoters
         strcpy(chrID,positions[r].chrID);
         if(r< WindowSize/2){ r = r + WindowSize/2;}  //avoid picking positions at the end or the begining of a chromosome
         if(r> no_position - WindowSize/2){r = r - WindowSize/2;}
         f = 1;
         flag = 1;
         while(flag && (f<= WindowSize/2)){ 
            if (strcmp(chrID,positions[r+f].chrID) || strcmp(chrID,positions[r-f].chrID) ) // if the closest WindowSize positions doesn't belong to the same chromosome
               {flag=0;}
            else  // the closest [-WindowSize/2 ~ +WindowSize/2] positions belong to the same chromosome,
               {f++;}
         }
         if(flag){ // current position is a non-enhancer, non-promoter position
            for(f=0;f< WindowSize;f++){
               fprintf(ftraining,"1");
               for(j=0;j<no_hist;j++){
                  fprintf(ftraining,"\t%f\t%f",positions[r+f-WindowSize/2].histvalues[j],positions[r+f-WindowSize/2].histenergy[j]);
               }
               fprintf(ftraining,"\n");
            }
            t++;
         }
         else {e--;} //pick up another position
      }
      else{e--;} //pick up another position
   } // write the non-enhancer positions
   //printf("Randomly choose %d non-enhancer positions for training.\n",t);
   free(randv);

   e=0;
   for(k=0;k<no_position;k++){
      //t = no_enhancer * 10;
      if(positions[k].enhancer_flag == 2){  // k is the enhancer center
         for(f=0;f<WindowSize;f++){
            fprintf(ftraining,"2");
            for(j=0;j<no_hist;j++)
               fprintf(ftraining,"\t%f\t%f",positions[k+f-WindowSize/2].histvalues[j],positions[k+f-WindowSize/2].histenergy[j]);
            
            fprintf(ftraining,"\n");
         }
         e++;
      }
   } // write the enhancer positions
   fclose(ftraining);
   printf("Done with writing file training.txt!\n");

}

void writetesting( positions, no_position, no_hist) // remember to zero the histone value for promoter positions
struct_position *positions;
int no_position;
int no_hist;
{
   int i;
   long k;
   FILE *ftesting;
   ftesting = fopen("testing.txt","w"); 
   if(ftesting==NULL){perror("testing.txt");}
   //printf("\nStart to write file testing.txt...\n");
   fprintf(ftesting,"ID\tposition");
   for(i=0;i<no_hist;i++){fprintf(ftesting,"\tHist%d\tHist%d_Energy",i+1, i+1);} 
   fprintf(ftesting,"\n"); //done writing the header
   for(k=0;k<no_position;k++){ //write each row when there's a promoter position write a whole zero row
      if(positions[k].promoter_flag == 1){ 
         fprintf(ftesting,"%s\t%ld",positions[k].chrID,positions[k].positionID);
         for(i=0;i<no_hist;i++){fprintf(ftesting,"\t0\t0");}
         fprintf(ftesting,"\n");
      }
      else{
         fprintf(ftesting,"%s\t%ld",positions[k].chrID,positions[k].positionID);
         for(i=0;i<no_hist;i++){fprintf(ftesting,"\t%f\t%f",positions[k].histvalues[i],positions[k].histenergy[i]);}
         fprintf(ftesting,"\n");
      }
   }
   printf("Done with writing file testing.txt!\n");
} 

