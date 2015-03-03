#include <math.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include "helper.h"
#include "model.h"
#include "kalman.h"
#include "simulate.h"
#include "mpccontrol.h"



//#include <mgl2/mgl_cf.h>
//#include <GL/freeglut.h>







//discretizes individual matrices.
void discretize(gsl_matrix *Ac,gsl_matrix *Bc,gsl_matrix *Cc,gsl_matrix *Ad,gsl_matrix *Bd,gsl_matrix *Cd,double Ts)
{
/*** Uses approximate discretization
Ad=(Id+Ac*Ts)
Bd=Ts*Bc
*/

int Ns=Ac->size1;
int Nu=Bc->size2;
int Ny=Cc->size1;

//We need 1 intermediate matrices Aint=Ac*Ts
gsl_matrix *Aint = gsl_matrix_alloc (Ns,Ns ); //allocate memory for matrix

//Set the values of identity matrix
gsl_matrix_set_identity (Ad); //Ad=Id
gsl_matrix_memcpy(Aint,Ac);//Aint=Ac
gsl_matrix_scale (Aint,Ts); //Aint=Aint*Ts
gsl_matrix_add(Ad,Aint);//Ad=Id+Aint*Ts

gsl_matrix_memcpy(Bd,Bc);//Bd=Bc
gsl_matrix_scale (Bd,Ts); //Bd=Bc*Ts

gsl_matrix_free(Aint);
}










// Linear Kalman Filter






void gnu_plot(char *c)
{
char d[100]="plot '";

//strcat(d,c);
//strcat(d,"' using 1:2 with lines\n");
FILE *pipe = popen("gnuplot -persist","w");
//fprintf(pipe, "set terminal wxt\n");
//fprintf(pipe, "set output 'test.png'\n");
//fprintf(pipe,"set multiplot\n");
//fprintf(pipe,d);
//fprintf(pipe,"plot 'test.dat' using 1:2 with lines\n");
//fprintf(pipe,"plot 'kalman.dat' using 1:2 with lines, 'test.dat' using 1:2 with lines\n");
//fprintf(pipe,"plot 'kalman.dat' using 1:7 with lines\n");
//fprintf(pipe,"plot 'kalman.dat' using 1:4 with lines, 'test.dat' using 1:4 with lines\n");
fprintf(pipe,c);
fclose(pipe);
}
//Read jacobian file and store it in struct model.








int main (void)
{

double a[]={0,1,-2,-3};
double b[]={0,1};
double c[]={1,0};
double d[]={0};

int Ns,Nu,Ny;


  structMPC mpc,*mpcptr;
  mpcptr=&mpc;
  mpcptr->A=gsl_matrix_alloc(2,2);
  mpcptr->B=gsl_matrix_alloc(2,1);
  mpcptr->C=gsl_matrix_alloc(1,2);
  mpcptr->D=gsl_matrix_alloc(1,1);

  assign_Mat(mpcptr->A,a);
  assign_Mat(mpcptr->B,b);
  assign_Mat(mpcptr->C,c);
  assign_Mat(mpcptr->D,d);

  mpcptr->predHor=10;
  mpcptr->contHor=5;

  Ns=mpcptr->A->size1;
  Nu=mpcptr->B->size2;
  Ny=mpcptr->C->size1;



 mpcptr->Q=gsl_matrix_alloc(Ny,Ny);
 gsl_matrix_set_identity(mpcptr->Q);
 mpcptr->P=gsl_matrix_alloc(Ny,Ny);
 gsl_matrix_memcpy(mpcptr->P,mpcptr->Q);
 mpcptr->R=gsl_matrix_alloc(Nu,Nu);
 gsl_matrix_set_identity(mpcptr->R);


 printf("MPC Prediction type%d %",mpcptr->predtype,mpcptr->type);

  mpcptr->predtype=OUTPUT;
  printf("MPC Prediction type%d %d",mpcptr->predtype,mpcptr->type);
  MPCpredmat(mpcptr,OUTPUT);


return 0;
}


/** Setting of parameters

//Model parameter does not need anything

//Kalman Filter needs model+disturbance model +Q+P


*/

