#include "model.h"
#include "kalman.h"
typedef struct{
gsl_matrix *state; //matrix to store simulation data.....dimensions are [Ndatapoints*statevectors]
gsl_matrix *output;
double sampling_time;   //sampling time
double simulation_time; //simulation time
int Ndatapoints;       //number of datapoints to be simulated
int initflag;
double *time_vec;       //flag to indicate is simulation has been initialised.
}Simdata;

void InitSim(Model *m,Simdata *sim,double Ts, double tsim);
void Simulate(Model *m, Simdata *msim, gsl_matrix *u,double Ts,double tsim);
