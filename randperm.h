//This program uses time function to seed random number generator
//and then generates random number to generate randperm vector of value size_vector

void randperm(vector, size_vector)
int *vector;
int size_vector;
{
   int i,r;
   int x;
   for(i=0;i<size_vector;i++)
      vector[i] = i;
   
   srand((unsigned)time(NULL));
 
   for(i=0;i<size_vector;i++){
      r = rand() % (size_vector-i);  // get a random value smaller than i
      x = vector[i];
      vector[i] = vector[r+i];
      vector[r+i] = x;  // exchange the value of vector[i] and vector[r];
   } 
}

