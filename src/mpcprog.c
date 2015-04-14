#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_linalg.h>
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

void wait()
{
    fflush(stdin);
   char c;
    scanf("%d",&c);
}




void gnu_plot(char *c,char *e)
{
char d[100]="plot '";

strcat(d,c);
strcat(d,"' using ");
strcat(d,e);
strcat(d," with lines\n");
FILE *pipe = popen("gnuplot -persist","w");
//fprintf(pipe, "set terminal wxt\n");
//fprintf(pipe, "set output 'test.png'\n");
//fprintf(pipe,"set multiplot\n");
fprintf(pipe,d);
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

char *model_file="airplane.jac0";

FILE *fmpc;

///Specific constraints
double EL_ANGLE_MIN=-0.262;
double EL_ANGLE_MAX=0.262;

double EL_SLEW_MAX=0.524;
double EL_SLEW_MIN=-0.524;

double PITCH_ANGLE_MAX=0.349;
double PITCH_ANGLE_MIN=-0.349;

double ALTITUDE_MIN=-100000;
double ALTITUDE_MAX= 100000 ;

double ALTITUDE_RATE_MAX= 100000;
double ALTITUDE_RATE_MIN=-100000;


///constraints
double lbdelta[]={EL_SLEW_MIN};
double ubdelta[]={EL_SLEW_MAX};
double lbu[]={EL_ANGLE_MIN};
double ubu[]={EL_ANGLE_MAX};
double lbx[]={-100000,PITCH_ANGLE_MIN,-100000,-10000}  ;
double ubx[]={-100000,PITCH_ANGLE_MAX,-100000,-10000};
double lby[]={PITCH_ANGLE_MIN,ALTITUDE_MIN,ALTITUDE_RATE_MIN};
double uby[]={PITCH_ANGLE_MAX,ALTITUDE_MAX,ALTITUDE_RATE_MAX};



int Ns=4,Nu=1,Ny=3,Np,Nc;
int i,j; ///counters

///model
Model model,*modelptr,modeld,*modeldptr;
structMPC mpc,*mpcptr;



/**
***
***          | NORMAL      | DELTA
---------------------------------------------
      STATES | Q=[Qx]  R=R | Q=[Qx R] R=Rrate|
---------------------------------------------
      OUTPUT | Q=[Qy]  R=R | Q=[Qy R] R=Rrate|
---------------------------------------------

/**************ASSIGNMENT OF VARIABLES*************/
Np=10; ///predition horizon
Nc=3; ///control horizon

mpcptr=&mpc;
modelptr=&model;
modeldptr=&modeld;
ReadJac(modelptr,model_file);
//LoadDoubles(modelptr,a,b,c,d,x0,Ns,Nu,Ny);
discretize_model(modelptr,modeldptr,0.5);
//DiscrModel(modelptr,modeldptr,0.5);
printf("Model\n");
printf("A matrix\n");
print2scr(modelptr->A);
print2scr(modeldptr->A);
printf("B matrix\n");
print2scr(modelptr->B);
print2scr(modelptr->C);
print2scr(modelptr->D);


///AssignMPCweights(mpcptr,Q,R,Rrate); This has to be set.
double qx[]={1,1,1,1};
double qy[]={1,1,1};
double r[]={1};
double rrate[]={1};
printf("I am here at 190");
//
/*createDiagonal(Qx,qx);
createDiagonal(Qy,qy);
createDiagonal(R,r);
createDiagonal(Rrate,rrate);
printf("Qx");
print2scr(Qx);

///this has to be done in a function. I am exceptionally doing it here
mpcptr->Q=gsl_matrix_alloc(Qy->size1,Qy->size2);
mpcptr->P=gsl_matrix_alloc(Qy->size1,Qy->size2);
mpcptr->R=gsl_matrix_alloc(Rrate->size1,Rrate->size2);

gsl_matrix_memcpy(mpcptr->Q,Qy);
gsl_matrix_memcpy(mpcptr->P,Qy);
gsl_matrix_memcpy(mpcptr->R,Rrate);

print2scr(mpcptr->Q);

print2scr(mpcptr->P);

print2scr(mpcptr->R);*/


///********************************************************************

InitMPCType(mpcptr,modeldptr,DELTA,OUTPUT); ///sets model and type of formulation and type of prediction
AssignMPCWeights(mpcptr,qy,r,rrate);
MPCpredmat(mpcptr,Np,Nc);
InitMPCconstraints(mpcptr,lbu,ubu,lby,uby,lbdelta,ubdelta);
print2FileMPC(mpcptr,"test1.txt");

///states=[xxx pitch angle xxx altitude]
///outputs=[pitch angle altitude altitude rate]
///we want to control the altitude so
double Cref[]={0,0,0,1};
InitSteadyState(mpcptr,Cref,1);
printf("Steady State Matrix\n");
print2scr(mpcptr->SteadyState);
print2scr(mpcptr->C);


printf("Constant Constraints Matrices\n");
for(i=0;i<Nc*Nu;i++)
    printf("%d:%f<=u<=%f\n",i,mpcptr->lb[i],mpcptr->ub[i]);


    for(i=0;i<Np*mpcptr->C->size1;i++)
    printf("%d%f<=Su.x<=%f\n",i,mpcptr->lbA[i],mpcptr->ubA[i]);


//printf("nC:%d nV:%d",mpcptr->nVar,mpcptr->nCons);
//print2scr(mpcptr->C);
//print2scr(mpcptr->H);
//print2scr(mpcptr->F);
//print2scr(mpcptr->G);

/*****SIMULATION**********/
double Ts=0.5;
double tsim=20;
gsl_matrix *Bd=gsl_matrix_alloc(modeldptr->A->size1,1);
gsl_matrix_set_zero(Bd);
print2scr(Bd);
InitSimModel(modeldptr,modeldptr->X0,Ts,tsim,Bd);

gsl_matrix *umat=gsl_matrix_alloc(Nu,1);
double u[1]={0};
printf("Steady\n");
assign_Mat(umat,u);


/****MPC Parameters to pass for stepping*/
double refr[]={400};
double inputdist[]={0.0};
double Bdx[]={0,0,0,0,0};
double outputdist[]={0,0,0,0};
double xref[]={0.0,0.0,0.0,40,0.0};

printf("Steady\n");
print2scr(mpcptr->SteadyState);
double deltaU[]={0};
double x[5]={0};
int k,l;
printf("Before Loop:Curr X data");
print2scr(modeldptr->currxdata);



//for(i=0;i<Nc*Nu;i++)
//    printf("%d:%f<=u<=%f\n",i,mpcptr->lb[i],mpcptr->ub[i]);
//
//
//    for(i=0;i<Np*mpcptr->C->size1;i++)
//    printf("%d%f<=Su.x<=%f\n",i,mpcptr->lbA[i],mpcptr->ubA[i]);

gsl_matrix *deltaUdata=gsl_matrix_alloc(1,modelptr->Ndatapoints);
for(i=1;i<modeldptr->Ndatapoints;i++)
{
    StepSteadyState(mpcptr,refr,inputdist,outputdist,Bdx,0);
    for(k=0;k<4;k++)
       x[k]=gsl_matrix_get(modeldptr->currxdata,k,0)-mpcptr->xss[k];
      x[4]=u[0]-mpcptr->uss[0];
    StepMPC(mpcptr,x,deltaU);
    gsl_matrix_set(deltaUdata,0,i,deltaU[0]);
    //wait();
    u[0]=u[0]+deltaU[0];

    assign_Mat(umat,u);
    ModelStep(modeldptr,i,umat,0.0);

}
print2scr(modeldptr->statedata);
print2scr(modeldptr->inputdata);
printModeldata(modeldptr,2,"model.txt");
gnu_plot("model.txt","1:2");
gnu_plot("model.txt","1:3");
gnu_plot("model.txt","1:4");
print2scr(mpcptr->Q);
for(i=0;i<120;i++)
    printf("suval[%d]:%f\n",i,mpcptr->suval[i]);

for(i=0;i<Nc*Nu;i++)
    printf("%d:%f<=uss<=%f\n",i,mpcptr->lbuss[i],mpcptr->ubuss[i]);


    for(i=0;i<Np*mpcptr->C->size1;i++)
    {
        printf("%d%f<=Su.xss<=%f\n",i,mpcptr->lbAxss[i],mpcptr->ubAxss[i]);
        printf("%d%f<=Su.xss<=%f\n",i,mpcptr->lbA[i],mpcptr->ubA[i]);
    }

print2scr(mpcptr->Su);



return 0;
}


/** Setting of parameters

//Model parameter does not need anything

//Kalman Filter needs model+disturbance model +Q+P


*/
/**
Prints the different matrices after initialisatioin
*/

