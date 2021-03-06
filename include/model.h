

#include "helper.h"

#ifndef MODEL_H
#define MODEL_H

typedef struct{
gsl_matrix *A;
gsl_matrix *B;
gsl_matrix *C;
gsl_matrix *D;
gsl_matrix *X0; //initial conditions.
gsl_matrix *statedata;
gsl_matrix *outputdata;
gsl_matrix *inputdata;
gsl_matrix *Bd;
int initfalg; //if model already initialised return 3
double Ts; //sample time
int step_value; //which time step we are
gsl_matrix *currxdata; //the value of the currentdata
gsl_matrix *currydata;

long Ndatapoints;
}Model;


void InitModel(Model *m,gsl_matrix *A,gsl_matrix *B,gsl_matrix *C,gsl_matrix *D,gsl_matrix *X0);
void InitSimModel(Model *m,gsl_matrix *X0,double Ts,double tsim,gsl_matrix *Bd);
void ModelStep(Model *m,int step_time,gsl_matrix *u,double dist);
void printModeldata(Model *m,int s,char *filename);
void ReadJac(Model *m,char *s);
void discretize_model(Model *mc,Model *md,double Ts);
void InitSimModelAme(Model *m,gsl_matrix *X0,double Ts,double tsim);
void LoadDoubles(Model *m,double *a,double *b, double *c, double *d, double *x0,int Ns,int Nu,int Ny);
void DiscrModel(Model *src,Model *dest, double Ts);
void DiscrDisturbance(Model *cont,double* Bdc,double* Bdk,int numDisturbances,double Ts);
#endif
