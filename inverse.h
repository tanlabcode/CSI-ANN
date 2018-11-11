/* inverse matrix Sw with maximum size of MaxSize * MaxSize */
/* Modified by Li Teng 2010/10 */

void arg(double *a,double *b,int *n,int x,int y)
{
   int k,l,i,j;
   for(i=0,k=0;i<*n;i++,k++)
   {
      for(j=0,l=0;j<*n;j++,l++)
      {
         if(i==x)
            i++;
         if(j==y)
            j++;
         *(b+ MaxSize*k+l)=*(a+MaxSize*i+j);

      }
   }
   *n=*n-1;
}

float det(double *p,int *n)
{
   double d[MaxSize][MaxSize],sum=0;
   int i,j,m;
   m=*n;
   if(*n==2)
   return(*p**(p+MaxSize+1)-*(p+1)**(p+MaxSize));
   for(i=0,j=0;j<m;j++)
   {
      *n=m;
      arg(p,&d[0][0],n,i,j);
      sum=sum+*(p+MaxSize*i+j)*pow(-1,(i+j))*det(&d[0][0],n);
   }

   return(sum);
}

/* input/output:
 *
 * float *Sw    the data matrix to be inversed is in Sw[n][n] / the result of inverse is in Sw[n][n]   
 * int n        the size of the input data matrix
 * int maxsize  the dimension of Sw is [maxsize][maxsize]
 */
void inv(double *Sw, int n, int maxsize)
{
   double d;
   int i,j,m;
   double a[MaxSize][MaxSize],b[MaxSize][MaxSize],c[MaxSize][MaxSize]; // for calculation of the inverse matrix of dimension lower than 10
   
   for(i=0;i<n;i++)
     for(j=0;j<n;j++)
        a[i][j] = Sw[i*maxsize + j];

  
   if(n==2)
   {
      c[0][0]=a[1][1];
      c[1][1]=a[0][0];
      c[0][1]=-a[0][1];
      c[1][0]=-a[1][0];
      d=a[0][0]*a[1][1]-a[0][1]*a[1][0];

      if(d==0)
      {
         getchar();
         exit(d-'0');
      }
   }
   else
   {
      m=n;
      for(i=0;i<m;i++)
      {
         for(j=0;j<m;j++)
         {
            n=m;
            arg(&a[0][0],&b[0][0],&n,i,j);
            c[j][i]=pow(-1,(i+j))*det(&b[0][0],&n);
         }
      }
      n=m;
      d=det(&a[0][0],&n);

      if(d==0)
      {
         printf("INVERSE DOES NOT EXIST");
         exit(d-'0');
      }
      for(i=0;i<m;i++)
      {
         for(j=0;j<m;j++)
            Sw[i*maxsize + j] = c[i][j]/d;
      }
   } 
}

