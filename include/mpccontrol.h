#include "model.h"
#include "qpOases_public.h"
#ifndef MPCCONTROL_H
#define MPCCONTROL_H

typedef enum
{
    DELTA=2,
    NORMAL=1
}MPCType; /** Determines whether to construct Delta Matrices or Normal Matrices*/

typedef enum
{
    STATE=2,
    OUTPUT=1
}MPCPredictionType;/** Determines whether Prediction types are on the states or the output*/

typedef struct
{
gsl_matrix *A;
gsl_matrix *B;
gsl_matrix *C;
gsl_matrix *D;
gsl_matrix *Q; //state weights
gsl_matrix *P;//terminal weights
gsl_matrix *R; //input wieights
gsl_matrix *Rrate; ///rate inputs
gsl_matrix *H; //
gsl_matrix *F;
gsl_matrix *G;
gsl_matrix *SteadyState;
gsl_matrix *constraints;
gsl_matrix *CM;
gsl_matrix *M;
gsl_matrix *Su;
gsl_matrix *Sx;
gsl_matrix *xdata;
gsl_matrix *statedata;
double *hval; /** size= cHorizon*Nu*cHorizon*Nu*/
double *gval;
double *fval;
int nVar; ///number of variables to optimise [used by qpoases]....usually the [control horizon*control varialbles] Nu*Nc
int nCons; ///number of constraints of the form lbX<u<lbU Np*Ns
double *lb,*lbuss; ///lbuss adjusts the values after steady state has been calculated.
double *ub,*ubuss;
double *lbA,*lbAxss;
double *ubA,*ubAxss;
double *lbAMX;
double *ubAMX;
double *CMval;
double *suval;
int predHor; ///prediction horizon
int contHor; ///control horizon
double Ts;
double *xss;
double *uss;
double *yss;
MPCType type; ///fomulation type
MPCPredictionType predtype; ///predition type
}structMPC;

typedef struct
{
    int size1;
    double *data;
}double_size;




void InitMPCType(structMPC *mpcptr,Model *modelptr,MPCType type,MPCPredictionType predtype);
int InitSteadyState(structMPC *mpcptr,double *Cref,int Ntr);
void InitMPC(structMPC *mpcptr,int cHorizon,gsl_matrix *Q,gsl_matrix *P,gsl_matrix *R,double *lb,double *ub,double *lbA,double *ubA);
void MPC_Step(structMPC *mpcptr,gsl_matrix *xdata, gsl_matrix *u);
/**returns the steady state values*/
double* MPCcalcSS(structMPC *mpcptr, double *refr,double *input_dist, double *output_dist,gsl_matrix *Bd,gsl_matrix *Cref);
int MPCpredmat(structMPC *mpcptr,int Np,int Nc);
int InitMPCconstraints(structMPC *mpcptr,double *lbu,double *ubu, double *lbxy, double *ubxy,double *lbdelta,double *ubdelta);

int AssignMPCWeights(structMPC *mpcptr,double *qxy,double *rinput,double *rrate);
int StepMPCconstraints(structMPC *mpcptr,double *xdata);
int StepSteadyState(structMPC *mpcptr,double *refr,double *inputdist,double *outputdist, double *Bd, int Nid);
int StepMPC(structMPC *mpcptr,double *x,double *u);
void print2FileMPC(structMPC *mpcptr,char *filename);

#endif
