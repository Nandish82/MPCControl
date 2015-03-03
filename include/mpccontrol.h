#include "model.h"
#include <qpOases_public.h>
#ifndef MPCCONTROL_H
#define MPCCONTROL_H

typedef enum
{
    DELTA=0,
    NORMAL
}MPCType; /** Determines whether to construct Delta Matrices or Normal Matrices*/

typedef enum
{
    STATE=0,
    OUTPUT,
    NOTDEFINED
}MPCPredictionType;/** Determines whether Prediction types are on the states or the output*/

typedef struct
{
gsl_matrix *A;
gsl_matrix *B;
gsl_matrix *C;
gsl_matrix *D;
gsl_matrix *Q; //state weights
gsl_matrix *P;//terminal weights
gsl_matrix *R;
gsl_matrix *H; //
gsl_matrix *F;
gsl_matrix *G;
gsl_matrix *SteadyState;
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
int predHor; ///prediction horizon
int contHor; ///control horizon
double Ts;
gsl_matrix *xdata;
gsl_matrix *statedata;
double *xss;
double *uss;
MPCType type;
MPCPredictionType predtype;
}structMPC;

typedef struct
{
    int size1;
    double *data;
}double_size;




void InitMPCType(structMPC *mpcptr,Model *modeldiscrete,MPCType type);
void InitSteadyState(structMPC *mpcptr,gsl_matrix *Cref);
void InitMPC(structMPC *mpcptr,int cHorizon,gsl_matrix *Q,gsl_matrix *P,gsl_matrix *R,double *lb,double *ub,double *lbA,double *ubA);
void MPC_Step(structMPC *mpcptr,gsl_matrix *xdata, gsl_matrix *u);
/**returns the steady state values*/
double* MPCcalcSS(structMPC *mpcptr, double *refr,double *input_dist, double *output_dist,gsl_matrix *Bd,gsl_matrix *Cref);
void MPCpredmat(structMPC *mpc,MPCPredictionType predtype);
#endif
