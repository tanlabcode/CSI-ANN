#include "headfile.h"
#include "extern.h"
#include "NeuralNetLib.h"

double *res;
double *yh1;
double *yh2;
double *x;
double *yout;


double OneLayersFFNN(int a, double *S, int N, int M, int NB, int NE, int *NNC)
{
	int j, b, c, i;
	double invdcoeff = 0;
	
	/*Initializing Arrays*/
	yh1[0] = 1;
    /*yout = 0;*/

    for(c = 0; c <= NNC[0]; c++)
		x[c] = 0;
	for(c = 1; c <= NNC[1]; c++)
		yh1[c] = 0;
	for(i = 0; i < (NB + NE); i++)
    {
		res[i] = 0;
        yout[i] = 0;
    }
	
	/*--> Calculation for Normal Patterns <--*/
	for(i = 0; i < N; i++)
	{
        /*Inputs*/
        x[0] = 1;
        for(j = 1; j <= NNC[0]; j++)
            x[j] = S[(i*NNC[0])+j-1];

        /*Evaluate NN*/
        /* --> First Hidden Layer <--*/
        for(b = 1; b <= NNC[1]; b++) /*row*/
            for(c = 0; c <= NNC[0]; c++)  /*column*/
                yh1[b] += xx[(b-1)*(NNC[0]+1)+c][a]*x[c]; /* before x[(b-1)*(NNC[0]+1)+c]*/
        /*Activation function*/
        for(c = 1; c <= NNC[1]; c++)
            yh1[c] = logsig(yh1[c]);

        /*--> Output Layer <--*/
        b = (NNC[0]+1)*NNC[1];
        for(c = 0; c <= NNC[1]; c++)
            yout[i] += xx[b+c][a]*yh1[c];
        res[i] = logsig(yout[i]);

        /*Clearing arrays*/
        for(c = 1; c <= NNC[1]; c++)
            yh1[c] = 0;
	}

	/*Classification*/
    /* No. of Background, No. of Enhancers */
	invdcoeff = error(NB, NE);
	
	return invdcoeff;
}

double TwoLayersFFNN(int a, double *S, int N, int M, int NB, int NE, int *NNC)
{
	int j, b, c, i;
	double invdcoeff = 0;
	
	/*Initializing Arrays*/
    x[0] = 1;
    yh1[0] = 1;
	yh2[0] = 1;

    for(c = 1; c <= NNC[0]; c++)
		x[c] = 0;
	for(c = 1; c <= NNC[1]; c++)
		yh1[c] = 0;
	for(c = 1; c <= NNC[2]; c++)
		yh2[c] = 0;
	for(i = 0; i < (NB + NE); i++)
    {
		res[i] = 0;
        yout[i] = 0;
    }
	
	/*--> Calculation for Normal Patterns <--*/
	for(i = 0; i < N; i++)
	{
        /*Inputs*/
        for(j = 1; j <= NNC[0]; j++)
            x[j] = S[(i*NNC[0])+j-1];

        /*Evaluate NN*/
        /* --> First Hidden Layer <--*/
        for(b = 1; b <= NNC[1]; b++) /*row*/
            for(c = 0; c <= NNC[0]; c++)  /*column*/
                yh1[b] += xx[(b-1)*(NNC[0]+1)+c][a]*x[c];
        /*Activation function*/
        for(c = 1; c <= NNC[1]; c++)
            yh1[c] = logsig(yh1[c]);

        /*--> Second Hidden Layer <--*/
        for(b = 1; b <= NNC[2]; b++)
            for(c = 0; c <= NNC[1]; c++) 
                yh2[b] += xx[(b-1)*(NNC[1]+1)+c+(NNC[0]+1)*NNC[1]][a]*yh1[c];

        /*Activation function*/
        for(c = 1; c <= NNC[2]; c++)
            yh2[c] = logsig(yh2[c]);

        /*--> Output Layer <--*/
        b = (NNC[0]+1)*NNC[1] + (NNC[1]+1)*NNC[2];
        for(c = 0; c <= NNC[2]; c++)
            yout[i] += xx[b+c][a]*yh2[c];
        res[i] = logsig(yout[i]);

        /*Clearing arrays*/
        for(c = 1; c <= NNC[1]; c++)
            yh1[c] = 0;
        for(c = 1; c <= NNC[2]; c++)
            yh2[c] = 0;
	}

	/*Classification*/
    /* No. of Background, No. of Enhancers */
	invdcoeff = error(NB, NE);
	
	return invdcoeff;
}


double OneLayerTDNN(int a, double *S, int N, int M, int NB, int NE, int *NNC, int PpW, int NoOfFeatures, int *TimeCopies, int *Delays)
{
   int j, b, c, i, IpN, Node, RepectFields = 0, count, TCL, PartIndex1 = 0, d;
   double invdcoeff = 0, *yh1, yout;
   double  aa =1.716, bb = (double) 2/3;

   /*Initializing Arrays*/
   yh1 = (double *)malloc((TimeCopies[0]*NNC[1])*sizeof(*yh1));

   yout = 0;

   for(c = 0; c < NNC[0]; c++)
      x[c] = 0;
   for(c = 0; c < TimeCopies[0]*NNC[1]; c++)
      yh1[c] = 0;
   for(i = 0; i < (NB + NE); i++)
      res[i] = 0;

   /*--> Calculation for Normal Patterns <--*/
   for(i = 0; i < N; i++)  //?? why segmentation fault here??
   //for(i=0;i<2588;i++)
   {
      /*Inputs*/
      for(j = 0; j < NNC[0]; j++)
         x[j] = S[(i*NNC[0])+j];

      /*================> Evaluate NN <================*/
      /* --> First Hidden Layer <--*/
      //printf(" On Watch, i = %d\n",i); // on watch

      RepectFields = Delays[0] + 1;
      IpN = RepectFields*NoOfFeatures;

      for(TCL = 0; TCL < TimeCopies[0]; TCL++)  /*Time-Copy Layer*/
         for(Node = 0; Node < NNC[1]; Node++)    /*Nodes*/
            for(d = 0; d < IpN; d++)
               yh1[TCL*NNC[1] + Node] += xx[Node*IpN + d][a]*x[NoOfFeatures*TCL+d];
   
        /*Activation function*/
      for(c = 0; c < TimeCopies[0]*NNC[1]; c++)
         yh1[c] = logsig(yh1[c]); ////$$$$$$$$$$$        
         //yh1[c] = aa * (exp(bb * yh1[c]) - exp(-bb *yh1[c])) /(exp(bb * yh1[c]) + exp(-bb *yh1[c]));

        /*--> Output Layer <--*/
      PartIndex1 = NoOfFeatures*RepectFields*NNC[1];
      RepectFields = Delays[1] + 1;
      IpN = RepectFields*NNC[1];

      for(c = 0; c < IpN; c++)
         yout += xx[PartIndex1+c][a]*yh1[c];

      res[i] = logsig(yout);//$$$$$$$$$$$$$$$$$
      //res[i] = aa * (exp(bb * yout) - exp(-bb * yout)) /(exp(bb * yout) + exp(-bb * yout));

      yout = 0;

      for(c = 0; c < TimeCopies[0]*NNC[1]; c++)
         yh1[c] = 0;
   }

   /*Classification*/
   /* No. of Background, No. of Enhancers */
   invdcoeff = error(NB, NE);
   free(yh1);
   return invdcoeff;
}

double TwoLayerTDNN(int a, double *S, int N, int M, int NB, int NE, int *NNC, int PpW, int NoOfFeatures, int *TimeCopies, int *Delays)
{
	int j, b, c, i, IpN, Node, RepectFields, Count, TCL, PartIndex1 = 0, d;
	double invdcoeff = 0, *yh1, yout;
	
	/*Initializing Arrays*/
 	yh1 = (double *)malloc((TimeCopies[0]*NNC[1])*sizeof(*yh1));
	/*yout = (double *)malloc(2*sizeof(*yout));*/

    
	/*yh1[0] = 1;*/
	/*yh2[0] = 1;*/
    yout = 0;
    /*yout[1] = 0;*/
    
    for(c = 0; c < NNC[0]; c++)
		x[c] = 0;
	for(c = 0; c < TimeCopies[0]*NNC[1]; c++)
		yh1[c] = 0;
/*	for(c = 0; c < TimeCopies[1]*NNC[2]; c++)
		yh2[c] = 0;*/
	for(i = 0; i < (NB + NE); i++)
		res[i] = 0;
	
	/*--> Calculation for Normal Patterns <--*/
	for(i = 0; i < N; i++)
	{
        /*Inputs*/
        for(j = 0; j < NNC[0]; j++)
            x[j] = S[(i*NNC[0])+j];

        
        /*================> Evaluate NN <================*/
            
        /* --> First Hidden Layer <--*/
        
        RepectFields = Delays[0] + 1;
        IpN = RepectFields*NoOfFeatures;
        /*Count = 0;*/

        for(TCL = 0; TCL < TimeCopies[0]; TCL++)  /*Time-Copy Layer*/
            for(Node = 0; Node < NNC[1]; Node++)    /*Nodes*/
                for(d = 0; d < IpN; d++)
                {
                    yh1[TCL*NNC[1] + Node] += xx[Node*IpN + d][a]*x[NoOfFeatures*TCL+d];
                }

        /*Activation function*/
        for(c = 0; c < TimeCopies[0]*NNC[1]; c++)
        {
            yh1[c] = logsig(yh1[c]);
/*                    if(i < 10)
                    printf("yh1[%d] = %f  ", c, yh1[c]);*/
            
        }
        
        /*--> Second Hidden Layer <--*/
        /*PartIndex1 = NoOfFeatures*RepectFields*NNC[1];
        RepectFields = Delays[1] + 1;
        IpN = RepectFields*NNC[1];

        Count = 0;
        for(TCL = 0; TCL < TimeCopies[1]; TCL++)  /*Time-Copy Layer*/
            /*for(Node = 0; Node < NNC[2]; Node++)    /*Nodes*/
            /*    for(d = 0; d < IpN; d++)
                {
                    yh2[TCL*NNC[2] + Node] += xx[PartIndex1 + Node*IpN + d][a]*yh1[NNC[1]*TCL+d];
                   if(Node == 0 && TCL == 0 && i == 0)
                        printf("yh2[%d] = %f ", d, yh1[NNC[1]*TCL+d]);
                }

        /*Activation function
        for(c = 0; c < TimeCopies[1]*NNC[2]; c++)
        {
            if(i == 0)
                printf("yh2[%d] = %f  ", d, yh2[c]);
            yh2[c] = logsig(yh2[c]);
        }*/

        /*--> Output Layer <--*/
        PartIndex1 = NoOfFeatures*RepectFields*NNC[1];
        RepectFields = Delays[1] + 1;
        IpN = RepectFields*NNC[1];

/*        if(i == 0)
            printf("\n IpN = %d", IpN);*/

        
        for(c = 0; c < IpN; c++)
        {
            yout += xx[PartIndex1+c][a]*yh1[c];
/*            yout[c] = xx[PartIndex1+c][a]*yh1[c];*/
/*            if(i == 0)*/
/*            if(i < 10)
                printf("yh2[%d] = %f  ", c, yh2[c]);*/
        }
/*        if(i < 10)
        printf("\n");*/
        res[i] = logsig(yout);
/*        printf("yout[%d] = %f  ", i, yout);*/
        /*Clearing arrays*/
        yout = 0;
        /*yout[1] = 0;        */
        for(c = 0; c < TimeCopies[0]*NNC[1]; c++)
            yh1[c] = 0;
/*        for(c = 0; c < TimeCopies[1]*NNC[2]; c++)
            yh2[c] = 0;*/
	}

	/*Classification*/
    /* No. of Background, No. of Enhancers */
	invdcoeff = error(NB, NE);

    free(yh1);
    /*free(yout);*/
/*    free(yh2);*/
	
	return invdcoeff;

}

void LoadingData(int NoBg, int NoEn, int *NNC) /*Function that read data for from files for the problem*/
{
	int n;

	n = NoBg + NoEn;
	res = (double *)malloc(n*sizeof(*res));
	x = (double *)malloc((NNC[0]*NNC[1]+1)*sizeof(*x));
	yh1 = (double *)malloc((NNC[1]+1)*sizeof(*yh1));
    yh2 = (double *)malloc((NNC[2]+1)*sizeof(*yh2));
    yout = (double *)malloc(n*sizeof(*yout));

    
}

void FreeMem(void)
{
	free(res);
    free(x);
    free(yh1);
    free(yout);
    free(yh2);
}

double ifdr(int nnp, int nfp)
{
	int i;
	double mean1 = 0, mean2 = 0, var1 = 0, var2 = 0, FDRI = 0;
	
	/*--> Means <--*/
	for (i = 0; i < nnp; i++)
		mean1 += res[i]/nnp;
	
	for (i = nnp; i < (nnp+nfp); i++)
		mean2 += res[i]/nfp;
	
	/*--> Unbiased Variances <--*/
	for (i = 0; i < nnp; i++)
		var1 = var1 + pow((res[i] - mean1), 2)/(nnp - 1);
	
	for (i = nnp; i < (nnp+nfp); i++)
		var2 = var2 + pow((res[i] - mean2), 2)/(nfp - 1);
	
	/*Calculating the Inverse of Fisher Discriminant Ratio: (FDR)^-1*/
	FDRI = sqrt((var1 + var2 + .00000000001)/2)/(fabs(mean1 - mean2));
    /*printf("FDRI: %.4f\n", FDRI);*/
	return FDRI;
}

double error(int NoBg, int NoEn)
{
    int i, N = 0;
    double Error = 0;
    //double aa = 1.716;

    N = NoBg + NoEn;
    for(i = 0; i < NoBg+NoEn; i++)
    {
        if(i < NoBg)
            Error += pow(res[i],2)/N;   // backgrand target is 0 
        else
            Error += pow(1 - res[i],2)/N;  // enhancer target is 1  
    }

    /*for(i = 0; i < NoBg+NoEn; i++)
    {
        if(i < NoBg)
            Error += pow(-aa -res[i],2)/N;   // backgrand target is -aa 
        else
            Error += pow(aa - res[i],2)/N;  // enhancer target is aa  
    }*/

    /*ErrorRisk = ClasfError[0]/NoBg + ClasfError[1]/NoEn;*/
    return Error;
}

double logsig(double v)
{
	double a, d;
	double f;

/* 	a = 1.716;
 	d = 0.666666666666666667;
    f = a*tanh(d*v);*/
    f = 1/(1+exp(-v));

	/*df = a*b*sech(b*x).^2;*/

	return f;
}

double act(double v)
{
	double a, d;
	double f;

 	a = 1.716;
 	d = 0.666666666666666667;
    f = a*tanh(d*v);
/*    f = 1/(1+exp(-v));*/

	/*df = a*b*sech(b*x).^2;*/

	return f;
}

double hardlim(double v)
{
	double f;

    if (v >= 0)
        f = 1;
    else
        f = 0;
    
	return f;
}

double hardlim2(double v)
{
	double f, x;
    x = abs(v)/10;

    if (x >= 0.5)
        f = 1;
    else
        f = 0;
    
	return f;
}





