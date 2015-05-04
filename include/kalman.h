#include "model.h"

#ifndef KALMAN_H
#define KALMAN_H

typedef struct
{
gsl_matrix *A;
gsl_matrix *B;
gsl_matrix *C;
gsl_matrix *D;
gsl_matrix *Q; //process covariance
gsl_matrix *R; //measurement covariance
gsl_matrix *P; //error covariance
gsl_matrix *K; //Kalman Gain
int initflag;
double Ts;
gsl_matrix *xdata;
gsl_matrix *statedata;
}Kalman_struc;

void Kalman_Init(Kalman_struc *Kalman,Model *m,gsl_matrix *Qcov, gsl_matrix *Rcov);
void Kalman_Step(Kalman_struc *Kalman,gsl_matrix *ymeas,gsl_matrix *u);
void Kalman_Init_d(Kalman_struc *Kalman,Model *m);
void Kalman_Init_d_o(Kalman_struc *Kalman,Model *m);
void Kalman_Step2(Kalman_struc *Kalman,gsl_matrix *ymeas,gsl_matrix *u);
void printKalmanData(Kalman_struc *m,int s,char *filename);
void Kalman_Init_d_servo(Kalman_struc *Kalman,Model *m);
void Kalman_Init_Ame_servo(Kalman_struc *Kalman,Model *m,double Ts);
void Kalman_Step2Ame(Kalman_struc *Kalman,gsl_matrix *ymeas,gsl_matrix *u,int step_value);
void Kalman_Init_Dist(Kalman_struc *Kalman,Model *m,double Ts,double *Bd,double *Dd,int Ndi,int Ndo,double *Qkal,double *Rcal);
int Controllability(Kalman_struc *Kalman); ///returns rank of controllability matrix
int Observability(Kalman_struc *Kalman); ///returns rank of observability matrix
void findSVD(gsl_matrix *A,gsl_vector *S);
int matRank(gsl_matrix *R,double tol);
#endif
