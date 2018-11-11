#include <stdio.h>
#include <string.h> 
#include <stdlib.h>

#include "step1.h"


main() {
   int i,j, flag, no_enhancer = 0, no_hist = 0, no_gene =0, no_chr;
   long no_position = 0; 
   char enhancer_fn[50], gene_fn[50], histone_fn[50], histone[50], chrID[10], buffer[300];// when using fgets+buffer to read a whole line, buffer has to be larger than a row
   char *ptr;
   long txStart,txEnd; //temporay para
   char bit; //temporary para
   long chroms[30];
   long no_within; // temp para
   char convert[10];

   FILE *fenhancer, *fgene, *fhistone, *fhistvalue;
   struct_enhancer enhancers[2000];
   struct_position *positions;
   struct_gene *genes;
   
//----------------------------------------- Read Enhancer File ----------------------------------------//
   printf("PLEASE INPUT THE FILE NAME FOR THE ENHANCER MARKERS: ");
   scanf("%s", enhancer_fn);
   fenhancer=fopen(enhancer_fn, "r");
   if (fenhancer == NULL){
      perror(enhancer_fn);
   }
   else{
      i = 0;
      while (!feof(fenhancer)) {
         fscanf(fenhancer,"%s %ld", &enhancers[i].chrID, &enhancers[i].center); 
 
         if(enhancers[i].chrID[0] != 'c'){  // handling diffenent format of enhancer marker file
            strcpy(convert, "chr\0");
            strcat(convert, enhancers[i].chrID);
            strcpy(enhancers[i].chrID, convert);
            //printf("convert = %s, enhancers[%d].chrID = %s\n", convert, i, enhancers[i].chrID);
         }

         if (fgets(buffer,500,fenhancer)!=NULL) {
            enhancers[i].enhancer_start = enhancers[i].center - HalfWindowSize;
            enhancers[i].enhancer_end = enhancers[i].center + HalfWindowSize;
            i++;
         }
         else{ break;}
      }
      fclose(fenhancer);
      no_enhancer = i;
      printf ("Done with reading information for %d enhancers.\n\n", no_enhancer);
   }

//------------------- Read Histone Modification Files ------------------------------//
/* To construct the structure array we have to first get the number of the positions first
   So the first file of reading the first file, it to get the number of the positions.
   Then we read the first file again and also the other files to construct the value of the 
   structures.*/
   no_chr = 0; // count the number of chromos
   flag = 1; // a symbol of first reading
   printf("PLEASE INPUT THE FILE NAME FOR THE HISTONE MODIFICATIONS: ");
   scanf("%s", histone_fn);
   //------get the number of positions ---------------------------------------------//
   fhistone = fopen(histone_fn,"r");
   if (fhistone == NULL){
      perror(histone_fn);
   }
   else{ // first time read the histvalue files, need to get the number of the positions
      fscanf(fhistone, "%s", &histone);   // get the first file name stored in histone
      if (fgets(buffer,300,fhistone) == NULL){ printf("The file is empty.\n", histone_fn);}
      else {
         fhistvalue = fopen(histone,"r"); // open the first file 
         if(fhistvalue==NULL){ perror(histone);}
         else{
            while(!feof(fhistvalue) && (fgets(buffer,300,fhistvalue)!= NULL)){
               if ((ptr = strstr(buffer, "track") ) != NULL ){/*description row no need to handle*/      }
               else if( (ptr = strstr(buffer, "chrom=") ) != NULL){
                  /*Start a new chromosome*/
                  if(flag) {
                     no_within = 0; 
                     flag = 0;
                  } /*description rows no need to handle */
                  else {
                     chroms[no_chr]= no_within; 
                     no_chr++; 
                     no_within=0; 
                  }
               }
               else { no_position++; no_within++;} // count the rows of position information
            }
            chroms[no_chr] = no_within;
            no_chr++; //record the last chromos
            fclose(fhistvalue);  // close the first file
         }
      }
      fclose(fhistone);
   }

   //------------construct the position structure and read information for the positions.-----------//
   positions = (struct_position *)malloc(no_position * sizeof(struct_position));
   if (positions == NULL) {printf("Allocate space for positions failed! \n");abort();}  // allocation space failed!
   else{
      fhistone = fopen(histone_fn,"r");
      if (fhistone == NULL){ perror(histone_fn); }  // open the file again to get the position information
      else {
         while(!feof(fhistone)){
            fscanf(fhistone, "%s", &histone);   // get the lines of the file name from histone_fn
            if (fgets(buffer,100,fhistone) != NULL){
               fhistvalue = fopen(histone,"r"); // open the first file 
               if(fhistvalue==NULL){
                  fclose(fhistone);
                  perror(histone);
               }
               else{
                  printf("Now reading file %s ... \n", histone);
                  i = 0; //index of the position
                  while(!feof(fhistvalue) && (fgets(buffer,200,fhistvalue)!= NULL)){
                     if( (ptr = strstr(buffer, "track") ) != NULL){ /*description rows no need to handle */ /*printf("Check! We are here.\n");*/}
                     else if( (ptr = strstr(buffer, "chrom=") ) != NULL){  //Start a new chromosome
                        ptr += 6; // increment the ip address following the colon
                        j=0;
                        while( *ptr != 32){  // 32 stands for blank
                           //printf("%d\t",*ptr);
                           chrID[j] = *ptr;
                           j++;
                           ptr ++;
                        }
                        chrID[j]='\0';
                     }
                     else { 
                        strcpy(positions[i].chrID,chrID);
                        fscanf(fhistvalue, "%ld %f", &positions[i].positionID, &positions[i].histvalues[no_hist]);//&positions[no_positions].histvalues[0] ); 
                        positions[i].histenergy[no_hist] = positions[i].histvalues[no_hist] * positions[i].histvalues[no_hist]; 
                        positions[i].enhancer_flag = 0;
                        positions[i].promoter_flag = 0;
                        i++;   //index of the positions
                     } 
                  }  // read one row of the histvalue file
                  fclose(fhistvalue);// finish reading one file
                  no_hist++; 
               }
            }
         }
         printf("Done with reading information from %d histones files for %ld positions.\n\n",no_hist,no_position);
         fclose(fhistone);
      }

   }

//---------------------------------  Read Gene Annotation File --------------------------------------------//
  //---------- first get the number of genes than allocate a space for store the genes---//
   printf("PLEASE INPUT THE FILE NAME FOR GENE ANNOTATION: ");
   scanf("%s", gene_fn);
   fgene = fopen(gene_fn,"r");
   if(fgene == NULL){
      perror(gene_fn);
   }
   else{
      while(!feof(fgene) && (fgets(buffer,300,fgene)!= NULL)){
         if( (ptr = strstr(buffer, "chr") ) != NULL){ no_gene++;}
      }
      fclose(fgene);
   }
   //------ build the gene structure and read the data into gene structures----------//
   genes = (struct_gene *)malloc(no_gene * sizeof(struct_gene));
   if (genes ==0) {abort();}  // allocation space failed!
   else{
      fgene = fopen(gene_fn,"r");
      if(fgene == NULL){
         perror(gene_fn);
      }
      else{
         i = 0;
         while(!feof(fgene) && (fgets(buffer,300,fgene)!=NULL)){
            if((ptr = strstr(buffer,"chr"))!=NULL){  //move ptr to the position of chrXX
               flag = 4; // we need to read four items from the buffer
               j=0;
               txStart = 0;
               txEnd = 0;
               while(flag){
                  if( (flag==4) && (*ptr != 9) ){  //for blank ASCII=32; for TAB ASCII =9;  if *ptr is not a number or a-Z symbol, we exclude it
                     if (*ptr < 90 || *ptr > 97){
                        chrID[j] = *ptr;
                        j++;
                        ptr++;  //move the pointer
                     }
                     else{
                        chrID[j]='\0';
                        j++;
                        ptr++;
                     } // these are some special processing
                  }
                  else if((flag==4)&& (*ptr == 9)){
                     chrID[j] = '\0';
                     strcpy(genes[i].chrID,chrID);
                     flag = 3;
                      ptr++;
                 }
                 else if((flag == 3) && (*ptr !=9)){
                    genes[i].strand = *ptr;
                    ptr++;
                 }
                 else if((flag == 3) && (*ptr ==9)){
                    flag = 2;
                    ptr++;
                 }
                 else if((flag == 2) && (*ptr !=9)){
                    txStart = txStart * 10 + (*ptr-48);
                    ptr++; 
                 }     
                 else if((flag == 2) && (*ptr ==9)){
                    genes[i].txStart = txStart;
                    flag = 1; 
                    ptr++;
                 }
                 else if((flag == 1) && (*ptr !=9)){
                    txEnd = txEnd * 10 + (*ptr -48) ;
                    ptr++; 
                 }
                 else if((flag == 1) && (*ptr ==9)){
                    genes[i].txEnd = txEnd;
                    flag = 0;
                 }  
              }
              i++;
            }
         }
      fclose(fgene);
      printf("Done with reading annotation information for %d genes!\n\n",no_gene);
      }
   }   // finish reading the gene information
   printf("Start genertating training and testing data ...\n");
   markpromoter( chroms, no_chr, positions, no_position, genes, no_gene);
   markenhancer( chroms, no_chr, positions, no_position, enhancers, no_enhancer);
   writetraining( positions, no_position, no_hist, no_enhancer);
   writetesting( positions, no_position, no_hist);

}






