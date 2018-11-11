/* Evaluate the TDNN Model */
/* input/output
   Data        the input data set sequence,the dimension is 10/--
   size_window size of the enhancer window/ 
   parweights  weights for the particles, currently there are 22 weights/--

*/

// size_window = 10; which is actually fixed in the code

double TDNNEva(double *Data, int size_window, double *parweights) 
{
   int i,j,k;
   double node1=0, node2=0, output = 0;

   for(i=0;i<size_window;i++){
      node1 = node1 + parweights[i] * Data[i];
      node2 = node2 + parweights[size_window + i] * Data[i];
   }
   
   node1 = 1/ (1 + exp(-1 * node1));
   node2 = 1/ (1 + exp(-1 * node2));
   
   output = node1 * parweights[2*size_window] + node2* parweights[2*size_window+1];

   return ( 1/(1 + exp(-1 * output)) );
 
}


