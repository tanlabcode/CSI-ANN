#ifndef __NEURALNETLIB_H__
#define __NEURALNETLIB_H__

extern void LoadingData(int, int, int *);
extern void FreeMem();
extern double OneLayersFFNN(int, double *, int, int, int, int, int *);
extern double TwoLayersFFNN(int, double *, int, int, int, int, int *);
extern double OneLayerTDNN(int, double *, int, int, int, int, int *, int, int, int *, int *);
extern double TwoLayerTDNN(int, double *, int, int, int, int, int *, int, int, int *, int *);
double ifdr(int , int);
double logsig(double);
double hardlim(double);
double hardlim2(double);
double act(double v);
double error(int, int);


#endif

