
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>


void print2File(gsl_matrix *m,double *time_v,char *filename);

void print2scr(gsl_matrix *m); //print matrice to screen

// Takes a Matrix, inverts it and put the result in the A-1
void matrixInvert(gsl_matrix *A,gsl_matrix *Ainv);

gsl_matrix *MatInv2(gsl_matrix *m); //returns m^{-1}


void MatMul(gsl_matrix *m, gsl_matrix *n,gsl_matrix *p ); //multiplies two matrices  p<-m*n
gsl_matrix* MatTrans(gsl_matrix *n); //return n^T

gsl_matrix* MatMul2(gsl_matrix *m, gsl_matrix *n ); //returns m*n

gsl_matrix* MatAdd2(gsl_matrix *m, gsl_matrix *n ); //returns m+n

gsl_matrix* MatSub2(gsl_matrix *m, gsl_matrix *n );//returns n-m

void MatAdd(gsl_matrix *m, gsl_matrix *n,gsl_matrix *p); //p=m+n;

void printsizeMat(gsl_matrix *m,char *s);
void assign_MatP(gsl_matrix *m, double **val);
void assign_Mat(gsl_matrix *m,double *val);
void assignMatMat(gsl_matrix *m,gsl_matrix *n);

gsl_matrix* MatMulrec(gsl_matrix *m,int p);
void print2FileMat(gsl_matrix *m,FILE *filename);
int Diag(gsl_matrix *OUT,gsl_matrix *IN1,gsl_matrix *IN2);///takes two matrices and concatenates them diagonally
int createDiagonal(gsl_matrix *OUT,double *in);
