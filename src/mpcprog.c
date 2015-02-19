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


FILE *fmpc;
Kalman_struc k,*k_p;
double Ts=0.01;
/*Simulation Time*/
double tsim=20.0;
long Npoints=tsim/Ts+1;
gsl_matrix *u=gsl_matrix_alloc(1,Npoints);
gsl_matrix_set_all(u,0); //all values of input to 1.
int i,j;
///Servo Simulation
Model servoC,servoD, *ptr_servoC,*ptr_servoD;
Simdata msim_servo,*msim_ptr_servo;
msim_ptr_servo=&msim_servo;
ptr_servoC=&servoC;
ptr_servoD=&servoD;
//Initialise model with file servojac
ReadJac(ptr_servoC,"ServoMech_Model_.jac0");
discretize_model(ptr_servoC,ptr_servoD,Ts);

gsl_matrix *Q=gsl_matrix_alloc(6,6);
gsl_matrix *R=gsl_matrix_alloc(2,2);
gsl_matrix_set_identity(Q);
gsl_matrix_set_identity(R);
//
 //Kalman_Init(Kalman_struc *Kalman,gsl_matrix *Qcov, gsl_matrix *Rcov,Model *m)
k_p=&k;
//print2scr(kalman_ptr->model->A);
//print2scr(kalman_ptr->Q);
gsl_matrix *X0=gsl_matrix_alloc(ptr_servoD->A->size1,1);
gsl_matrix_set_all(X0,0);
gsl_matrix *udata=gsl_matrix_alloc(ptr_servoD->B->size2,1);
//InitSimModel(ptr_servoD,X0,Ts,tsim);
gsl_matrix *ydis=gsl_matrix_alloc(ptr_servoD->C->size1,1);
gsl_matrix_set(ydis,0,0,50);
gsl_matrix_set(ydis,1,0,0.00);
double lb1[]={-1000};
double ub2[]={1000};
printf("I am here....\n");
structMPC mpckalman,*mpckalmanptr;
mpckalmanptr=&mpckalman;
double qmpcval2[]={1,1,1,1,1};
double pmpcval2[]={1,1,1,1,1};

double qmpcval[25];
double pmpcval[25];
for(i=0;i<25;i++)
{
   qmpcval[i]=0;
   pmpcval[i]=0;
}

for(i=0;i<5;i++)
{
    qmpcval[i*6]=qmpcval2[i];
    pmpcval[i*6]=pmpcval2[i];
}



gsl_matrix *Qmpckalman=gsl_matrix_alloc(ptr_servoD->A->size1,ptr_servoD->A->size1);
gsl_matrix_set_all(Qmpckalman,0);
gsl_matrix_set_identity(Qmpckalman);
assign_Mat(Qmpckalman,qmpcval);
gsl_matrix *Pmpckalman=gsl_matrix_alloc(ptr_servoD->A->size1,ptr_servoD->A->size1);
gsl_matrix_set_identity(Pmpckalman);
gsl_matrix_set_all(Pmpckalman,0);
assign_Mat(Pmpckalman,pmpcval);
double rmpcval[]={0.01};
gsl_matrix *Rmpckalman=gsl_matrix_alloc(ptr_servoD->B->size2,ptr_servoD->B->size2);
gsl_matrix_set_identity(Rmpckalman);
assign_Mat(Rmpckalman,rmpcval);
int cHorizon=10;
double lbmpckalman[]={-1000};
double ubmpckalman[]={1000};
double lbAmpc[]={-1000,-1000,-1000,-100,-1000};
double ubAmpc[]={1000,1000,1000,1000,1000};

/**disturbance matrix Td*/
double Bdval[]={0,0,0,-0.9549*Ts,0};
gsl_matrix *Bd=gsl_matrix_alloc(5,1);
assign_Mat(Bd,Bdval);
/**end*/

InitSimModel(ptr_servoD,X0,Ts,tsim,Bd);
Kalman_Init_d_servo(k_p,ptr_servoD);
InitMPCType(mpckalmanptr,ptr_servoD,NORMAL);
InitMPC(mpckalmanptr,cHorizon,Qmpckalman,Pmpckalman,Rmpckalman,lbmpckalman,ubmpckalman,lbAmpc,ubAmpc);
//InitStateConstraints(mpckalmanptr,lbAmpc,ubAmpc);
double qkalmanval[]={1,1,1,1,1,100};
gsl_matrix *qkalcov=gsl_matrix_alloc(6,6);
gsl_matrix_set_all(qkalcov,0);
for(i=0;i<6;i++)
    gsl_matrix_set(qkalcov,i,i,qkalmanval[i]);
gsl_matrix_memcpy(k_p->Q,qkalcov);


double x0mpc[]={0,0,0,0,0};
assign_Mat(X0,x0mpc);
//MPC_Step(mpckalmanptr,X0,u);




gsl_matrix *umpc=gsl_matrix_alloc(ptr_servoD->B->size2,1);
gsl_matrix_set_all(umpc,0);

double refr[]={10};
double input_dist[]={0};
double output_dist[]={0};

//gsl_matrix *Bd=gsl_matrix_alloc(mpckalmanptr->A->size1,1);
//gsl_matrix_set_all(Bd,0);
gsl_matrix *Cref=gsl_matrix_alloc(1,5);
double Crefval[]={0,0,0,0,1};
assign_Mat(Cref,Crefval);
double *ssval=malloc(6*sizeof(double));
ssval=MPCcalcSS(mpckalmanptr,refr,input_dist,output_dist,Bd,Cref);
gsl_matrix *mpcxdata=gsl_matrix_alloc(5,1);
assign_Mat(mpcxdata,x0mpc);

/**building disturbance matrix*/
gsl_matrix *torquedis=gsl_matrix_alloc(Npoints,1);
gsl_matrix_set_all(torquedis,0);

double torquedisval[Npoints];
  /**Adding disturbance*/
    for(i=0;i<Npoints;i++)
    if(i*Ts<10)
    {
        gsl_matrix_set(torquedis,i,0,0);
        torquedisval[i]=0;

    }

    else
    {
        gsl_matrix_set(torquedis,i,0,20.0);
        torquedisval[i]=20;
    }

double *mpcxdataval=malloc(ptr_servoD->A->size2*sizeof(double));

gsl_matrix_set_all(udata,50);

print2scr(Qmpckalman);
print2scr(Pmpckalman);
//Npoints=3;
int loop=1; //just to make the loop function for debuggin purposes
if(loop==1)
{

/**LOOP*/
for(j=1;j<Npoints;j++)
{
    if(j*Ts<3)
        refr[0]=0;
    else
        refr[0]=10;

    for(i=0;i<5;i++)
    gsl_matrix_set(mpcxdata,i,0,gsl_matrix_get(k_p->xdata,i,0));
    input_dist[0]=gsl_matrix_get(k_p->xdata,5,0);
    ssval=MPCcalcSS(mpckalmanptr,refr,input_dist,output_dist,Bd,Cref);
    /**recalculating the matrices to include the steady states*/
    for(i=0;i<ptr_servoD->A->size2;i++)
    {
        gsl_matrix_set(mpcxdata,i,0,gsl_matrix_get(ptr_servoD->currxdata,i,0)-ssval[i]);
        mpcxdataval[i]=gsl_matrix_get(ptr_servoD->currxdata,i,0)-ssval[i];
    }

    MPC_Step(mpckalmanptr,mpcxdata,umpc);
    for(i=0;i<ptr_servoD->B->size2;i++)
        gsl_matrix_set(udata,i,0,gsl_matrix_get(umpc,i,0)+ssval[i+ptr_servoD->A->size1]);

   ModelStep(ptr_servoD,j,udata,torquedisval[j]);
   //add disturbance to output
   //ptr_servoD->currydata=MatAdd2(ydis,ptr_servoD->currydata);
Kalman_Step2(k_p,ptr_servoD->currydata,udata,j);
if(j==1)
{
fmpc=fopen("MPCcodeBout2.txt","w");
fprintf(fmpc,"MPC--H\n");
print2FileMat(mpckalmanptr->H,fmpc);
fprintf(fmpc,"MPC--F\n");
print2FileMat(mpckalmanptr->F,fmpc);
fprintf(fmpc,"MPC--G\n");
print2FileMat(mpckalmanptr->G,fmpc);
fprintf(fmpc,"Kalman--A\n");
print2FileMat(k_p->A,fmpc);
fprintf(fmpc,"Kalman--B\n");
print2FileMat(k_p->B,fmpc);
fprintf(fmpc,"Kalman--C\n");
print2FileMat(k_p->C,fmpc);
fprintf(fmpc,"Kalman--D\n");
print2FileMat(k_p->D,fmpc);
fprintf(fmpc,"Kalman--Q\n");
print2FileMat(k_p->Q,fmpc);
fprintf(fmpc,"Kalman--R\n");
print2FileMat(k_p->R,fmpc);
fclose(fmpc);
}
}
printModeldata(ptr_servoD,0,"statedata.dat");
printKalmanData(k_p,0,"kalmandata2.dat");
gnu_plot("plot 'statedata.dat' using 1:6 with lines\n");
gnu_plot("plot 'kalmandata2.dat' using 1:7 with lines\n");
gnu_plot("plot 'statedata.dat' using 1:2 with lines\n");
gnu_plot("plot 'statedata.dat' using 1:5 with lines\n");
}



fmpc=fopen("MPCcodeBout.txt","w");
fprintf(fmpc,"Matrix H:\n");
print2FileMat(mpckalmanptr->H,fmpc);
fprintf(fmpc,"Matrix F:\n");
print2FileMat(mpckalmanptr->F,fmpc);
fprintf(fmpc,"Matrix G:\n");
print2FileMat(mpckalmanptr->G,fmpc);
fprintf(fmpc,"Matrix Q:\n");
print2FileMat(mpckalmanptr->Q,fmpc);
fprintf(fmpc,"Matrix R:\n");
print2FileMat(mpckalmanptr->R,fmpc);
fprintf(fmpc,"Prediction Horizon:%d\n",mpckalmanptr->contHor);
fprintf(fmpc,"Constrainsts: lb: %f ub: %f\n",mpckalmanptr->lb[5],mpckalmanptr->ub[5]);
for(i=0;i<5;i++)
    fprintf(fmpc,"xss[%d]:%f\n",i,mpckalmanptr->xss[i]);
for(i=0;i<1;i++)
    fprintf(fmpc,"uss[%d]:%f\n",i,mpckalmanptr->uss[i]);
fclose(fmpc);


structMPC mpc,*mpcptr;
mpcptr=&mpc;
InitMPCType(mpcptr,ptr_servoD,DELTA);
InitSteadyState(mpcptr,Cref);
print2scr(ptr_servoD->A);
print2scr(ptr_servoD->B);


print2scr(mpcptr->A);
print2scr(mpcptr->B);
print2scr(mpcptr->C);
print2scr(mpcptr->D);
print2scr(mpcptr->SteadyState);




return 0;
}


/** Setting of parameters

//Model parameter does not need anything

//Kalman Filter needs model+disturbance model +Q+P


*/

