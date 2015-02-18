#include "model.h"
#include <qpOases_public.h>
#ifndef MPCCONTROL_H
#define MPCCONTROL_H

typedef struct
{
gsl_matrix *A;
gsl_matrix *B;
gsl_matrix *C;
gsl_matrix *D;
gsl_matrix *Q; //state weights
gsl_matrix *R; //input weights
gsl_matrix *H;
gsl_matrix *F;
gsl_matrix *G;
gsl_matrix *constraints;
double *hval; /** size= cHorizon*Nu*cHorizon*Nu*/
double *gval;
double *fval;
double *nV;
double *nC;
double *lb;
double *ub;
double *lbA;
double *ubA;
double *lbAMX;
double *ubAMX;
gsl_matrix *CM;
gsl_matrix *M;
double *CMval;
int contHor;
double Ts;
gsl_matrix *xdata;
gsl_matrix *statedata;
double *xss;
double *uss;
}MPC_struc;

typedef struct
{
    int size1;
    double *data;
}double_size;

void InitMPC(MPC_struc *mpcptr,Model *m,int cHorizon,gsl_matrix *Q,gsl_matrix *P,gsl_matrix *R,double *lb,double *ub,double *lbA,double *ubA);
void MPC_Step(MPC_struc *mpcptr,gsl_matrix *xdata, gsl_matrix *u);
/**returns the steady state values*/
double* MPCcalcSS(MPC_struc *mpcptr, double *refr,double *input_dist, double *output_dist,gsl_matrix *Bd,gsl_matrix *Cref);
#endif
