#ifndef __DEFINE_H__
#define __DEFINE_H__

#define DESCENT 1.0

/* locate viewport                      */
#define RADIUS 1
#define VPLEFT 50
#define VPTOP 50
#define VPRIGHT 600
#define VPBOTTOM 400

/* to place objects within viewport     */
#define WIDTH VPRIGHT-VPLEFT
#define HEIGHT VPBOTTOM - VPTOP

typedef int *P_INT;
typedef P_INT IVECTOR;
typedef P_INT *IMATRIX;

typedef double *P_FLOAT;
typedef P_FLOAT FVECTOR;
typedef P_FLOAT *FMATRIX;

#endif
