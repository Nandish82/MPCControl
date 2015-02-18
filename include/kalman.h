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

void Kalman_Init(Kalman_struc *Kalman,gsl_matrix *Qcov, gsl_matrix *Rcov,Model *m);
void Kalman_Step(Kalman_struc *Kalman,gsl_matrix *ymeas,gsl_matrix *u);
void Kalman_Init_d(Kalman_struc *Kalman,Model *m);
void Kalman_Init_d_o(Kalman_struc *Kalman,Model *m);
void Kalman_Step2(Kalman_struc *Kalman,gsl_matrix *ymeas,gsl_matrix *u,int step_value);
void printKalmanData(Kalman_struc *m,int s,char *filename);
void Kalman_Init_d_servo(Kalman_struc *Kalman,Model *m);
void Kalman_Init_Ame_servo(Kalman_struc *Kalman,Model *m,double Ts);
void Kalman_Step2Ame(Kalman_struc *Kalman,gsl_matrix *ymeas,gsl_matrix *u,int step_value);
#endif
