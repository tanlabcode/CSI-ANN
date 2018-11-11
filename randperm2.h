//This program uses time function to randomly permute the items in vector

void randperm2(vector, size_vector)
float *vector;
int size_vector;
{
   int i,r;
   float x;
   
   srand((unsigned)time(NULL));
 
   for(i=0;i<size_vector;i++){
      r = rand() % (size_vector-i);  // get a random value smaller than i
      x = vector[i];
      vector[i] = vector[r+i];
      vector[r+i] = x;  // exchange the value of vector[i] and vector[r];
   } 
}

