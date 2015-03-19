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
/**********DECLARATION OF VARIABLES************************/
double a[]={-1.2822,0,0.98,0,0,0,1,0,-5.4293,0,-1.8366,0,-128.2,128.2,0,0};
double b[]={-0.3,0,-17,0};
double c[]={0,1,0,0,0,0,0,1,-128.2,128.2,0,0};
double d[]={0,0,0};
double x0[]={0,0,0,0};

///Specific constraints
double EL_ANGLE_MIN=-0.262;
double EL_ANGLE_MAX=0.262;

double EL_SLEW_MAX=0.524;
double EL_SLEW_MIN=-0.524;

double PITCH_ANGLE_MAX=0.349;
double PITCH_ANGLE_MIN=-0.349;

double ALTITUDE_MIN=-1000000;
double ALTITUDE_MAX=1000000;

double ALTITUDE_RATE_MAX=1000000;
double ALTITUDE_RATE_MIN=-1000000;


///constraints
double lbu[]={EL_SLEW_MIN};
double ubu[]={EL_SLEW_MAX};
double lbx[]={-100000,PITCH_ANGLE_MIN,-100000,-10000,EL_ANGLE_MIN}  ;
double ubx[]={-100000,PITCH_ANGLE_MAX,-100000,-10000,EL_ANGLE_MAX};
double lby[]={PITCH_ANGLE_MIN,ALTITUDE_MIN,ALTITUDE_RATE_MIN};
double uby[]={PITCH_ANGLE_MAX,ALTITUDE_MAX,ALTITUDE_RATE_MAX};



int Ns=4,Nu=1,Ny=3,Np,Nc;
int i,j; ///counters

///model
Model model,*modelptr,modeld;
structMPC mpc,*mpcptr;



///weight matrices
gsl_matrix *Qx=gsl_matrix_alloc(4,4);  ///state prediciton
gsl_matrix_set_identity(Qx);
gsl_matrix *Qy=gsl_matrix_alloc(3,3);  ///output prediciton
gsl_matrix_set_identity(Qy);
gsl_matrix *R=gsl_matrix_alloc(1,1);   ///input weight
gsl_matrix_set_identity(R);
gsl_matrix *Rrate=gsl_matrix_alloc(1,1); ///rate input weight.
gsl_matrix_set_identity(Rrate);

/**
***
***          | NORMAL      | DELTA
---------------------------------------------
      STATES | Q=[Qx]  R=R | Q=[Qx R] R=Rrate|
---------------------------------------------
      OUTPUT | Q=[Qy]  R=R | Q=[Qy] R=Rrate|
---------------------------------------------

/**************ASSIGNMENT OF VARIABLES*************/
Np=5; ///predition horizon
Nc=3; ///control horizon

mpcptr=&mpc;
modelptr=&model;
LoadDoubles(modelptr,a,b,c,d,x0,Ns,Nu,Ny);
discretize_model(modelptr,&modeld,0.5);

printf("Model\n");
print2scr(modelptr->A);
print2scr(modelptr->B);
print2scr(modelptr->C);
print2scr(modelptr->D);


///AssignMPCweights(mpcptr,Q,R,Rrate); This has to be set.
double qx[]={1,1,1,1};
double qy[]={1,1,1};
double r[]={1};
double rrate[]={1};
//
createDiagonal(Qx,qx);
createDiagonal(Qy,qy);
createDiagonal(R,r);
createDiagonal(Rrate,rrate);

///this has to be done in a function. I am exceptionally doing it here
mpcptr->Q=gsl_matrix_alloc(Qy->size1,Qy->size2);
mpcptr->P=gsl_matrix_alloc(Qy->size1,Qy->size2);
mpcptr->R=gsl_matrix_alloc(Rrate->size1,Rrate->size2);

gsl_matrix_memcpy(mpcptr->Q,Qy);
gsl_matrix_memcpy(mpcptr->P,Qy);
gsl_matrix_memcpy(mpcptr->R,Rrate);

print2scr(mpcptr->Q);

print2scr(mpcptr->P);

print2scr(mpcptr->R);

///********************************************************************

InitMPCType(mpcptr,&modeld,DELTA,OUTPUT); ///sets model and type of formulation and type of prediction
/// FUNCTION TO ASSIGN WEIGHTS MISSING
MPCpredmat(mpcptr,Np,Nc);
InitMPCconstraints(mpcptr,lbu,ubu,lby,uby);


///states=[xxx pitch angle xxx altitude]
///outputs=[pitch angle altitude altitude rate]
///we want to control the altitude so
double Cref[]={0,0,0,1};
InitSteadyState(mpcptr,Cref,1);
printf("Steady State Matrix\n");
print2scr(mpcptr->SteadyState);


printf("Constant Constraints Matrices\n");
for(i=0;i<Nc*Nu;i++)
    printf("%f<=u<=%f\n",mpcptr->lb[i],mpcptr->ub[i]);


    for(i=0;i<Np*Ny;i++)
    printf("%f<=Su.x<=%f\n",mpcptr->lbA[i],mpcptr->ubA[i]);


printf("nC:%d nV:%d",mpcptr->nVar,mpcptr->nCon);
print2scr(mpcptr->C);

//gsl_matrix *Qdelta=gsl_matrix_alloc(Qx->size1+R->size1,Qx->size1+R->size1);
//
/////assuming moodel with output prediciton type and delta formulation
//mpcptr->Q=gsl_matrix_alloc(mpcptr->C->size1,mpcptr->C->size1); ///change if state predicition is used
//mpcptr->R=gsl_matrix_alloc(mpcptr->B->size2,mpcptr->B->size2);
//mpcptr->P=gsl_matrix_alloc(mpcptr->C->size1,mpcptr->C->size1);
//gsl_matrix_set_identity(mpcptr->Q);
//gsl_matrix_set_identity(mpcptr->P);
//gsl_matrix_set_identity(mpcptr->R);
//print2scr(mpcptr->Q);
//print2scr(mpcptr->R);
//print2scr(mpcptr->R);
//printf(" I am here");
//MPCpredmat(mpcptr,Np,Nc);
//
//InitMPCconstraints(mpcptr,lbu,ubu,lby,uby);
//
//gsl_matrix *Cref=gsl_matrix_alloc(1,2);
//gsl_matrix_set_zero(Cref);
//
//gsl_matrix_set(Cref,0,0,1);
//
//
//InitSteadyState(mpcptr,Cref);
//
//
//printf("Steady State Matrix\n");
//print2scr(mpcptr->SteadyState);






//
//  structMPC mpc,*mpcptr;
//  mpcptr=&mpc;
//  mpcptr->A=gsl_matrix_alloc(2,2);
//  mpcptr->B=gsl_matrix_alloc(2,1);
//  mpcptr->C=gsl_matrix_alloc(1,2);
//  mpcptr->D=gsl_matrix_alloc(1,1);
//
//  assign_Mat(mpcptr->A,a);
//  assign_Mat(mpcptr->B,b);
//  assign_Mat(mpcptr->C,c);
//  assign_Mat(mpcptr->D,d);
//
//  mpcptr->predHor=Np;
//  mpcptr->contHor=Nc;
//
//  Ns=mpcptr->A->size1;
//  Nu=mpcptr->B->size2;
//  Ny=mpcptr->C->size1;
//
//
//
// mpcptr->Q=gsl_matrix_alloc(Ny,Ny);
// gsl_matrix_set_identity(mpcptr->Q);
// mpcptr->P=gsl_matrix_alloc(Ny,Ny);
// gsl_matrix_memcpy(mpcptr->P,mpcptr->Q);
// mpcptr->R=gsl_matrix_alloc(Nu,Nu);
// gsl_matrix_set_identity(mpcptr->R);
//
//mpcptr->predtype=OUTPUT;
//MPCpredmat(mpcptr);
//
//  printf("Su:\n");
//  print2scr(mpcptr->Su);
//
//double lbu[]={-1};
//double ubu[]={1};
//
//double lbx[]={-2,-3};
//double ubx[]={2,3};
//
//InitMPCconstraints(mpcptr,lbu,ubu,lbx,ubx);
//
//for(i=0;i++;i<Np*Ns)
//{
//    printf("\n");
//    for(j=0;j++;j<Nc*Nu)
//        printf("%5.0f ",mpcptr->suval[i*Nc*Nu+j]);
//}
//
////printf("Nu:%d Ns:%d Ny:%d Np:%d Nc:%d\n",Nu,Ns,Ny,Np,Nc);
//for(i=0;i<Nc*Nu;i++)
//    printf("%f<=u<=%f\n",mpcptr->lb[i],mpcptr->ub[i]);
//
//
//    for(i=0;i<Np*Ny;i++)
//    printf("%f<=Su.x<=%f\n",mpcptr->lbA[i],mpcptr->ubA[i]);

return 0;
}


/** Setting of parameters

//Model parameter does not need anything

//Kalman Filter needs model+disturbance model +Q+P


*/

