#include "simulate.h"

void InitSim(Model *m,Simdata *sim,double Ts, double tsim)
{
int rows=m->A->size1;
int cols=m->A->size2;
int output=m->C->size1;
//sim.statedata=gsl_matrix_alloc(3,3);

sim->sampling_time=Ts;
sim->simulation_time=tsim;
sim->Ndatapoints=tsim/Ts+1;
sim->state=gsl_matrix_alloc(sim->Ndatapoints,rows); ///be careful here the states are the columns and time is rows
sim->output=gsl_matrix_alloc(sim->Ndatapoints,output);
sim->initflag=1;
gsl_matrix_set_all(sim->state,0);
sim->time_vec=(double *)malloc(sim->Ndatapoints*sizeof(double));
for (cols=0;cols<sim->Ndatapoints;cols++)
sim->time_vec[cols]=0;
}
///////////////////simulation
void Simulate(Model *m, Simdata *msim, gsl_matrix *u,double Ts,double tsim)
{
gsl_matrix *tmpA=gsl_matrix_alloc(m->A->size1,1); // if 2 state present  this is  [x1 x2] tmp=A*x[k]
gsl_matrix *tmpB=gsl_matrix_alloc(m->A->size1,u->size1);
gsl_matrix *xdata=gsl_matrix_alloc(m->A->size1,1);
gsl_matrix *udata=gsl_matrix_alloc(m->B->size2,1);
gsl_matrix *udataDis=gsl_matrix_alloc(m->B->size2,1);
gsl_matrix *ydata=gsl_matrix_alloc(m->C->size1,1);
int i,j;

InitSim(m,msim,Ts,tsim);
///Kalman Filter Parameters
 gsl_matrix *Q=gsl_matrix_alloc(m->A->size1,m->A->size1);
gsl_matrix *R=gsl_matrix_alloc(m->C->size1,m->C->size1);
gsl_matrix_set_identity(Q);
gsl_matrix_set_identity(R);

Kalman_struc kman,*kmanptr;
kmanptr=&kman;
//Kalman_Init(kmanptr,Q,R,m);
Kalman_Init_d(kmanptr,m);
printf("Init_d initialised");
print2scr(kmanptr->A);
gsl_matrix *kalmandata=gsl_matrix_alloc(msim->state->size1,kmanptr->A->size2);
printf("I am here");

//Kalman_Init(Kalman_struc *Kalman,gsl_matrix *Qcov, gsl_matrix *Rcov,Model *m)
//Kalman_Init_d(kmanptr,m);//Kalman_Init(Kalman_struc *Kalman,gsl_matrix *Qcov, gsl_matrix *Rcov,Model *m)
gsl_matrix_set_all(kmanptr->P,0);
gsl_matrix_set_all(kmanptr->xdata,0);





printsizeMat(tmpA,"tmpA");
printsizeMat(tmpB,"tmpB");
printsizeMat(u,"u");
printsizeMat(xdata,"xdata");
printsizeMat(udata,"udata");
printsizeMat(m->A,"Matrix A");
printsizeMat(m->B,"Matrix B");
printf("No of Datapoints:%d",msim->Ndatapoints);

print2scr(m->A);
print2scr(m->B);
print2scr(m->C);
print2scr(m->D);

//Intial value of xdata
for(i=0;i<m->A->size1;i++)
	gsl_matrix_set(xdata,i,0,gsl_matrix_get(m->X0,i,0));

	printf("initial values xdata");
	print2scr(xdata);
/*initial value for time	*/
msim->time_vec[0]=0*Ts;
for(i=0;i<m->A->size1;i++)
{
gsl_matrix_set(msim->state,0,i,gsl_matrix_get(xdata,i,0)); //initial value of state for kalman filter
}

for(i=0;i<kmanptr->xdata->size1;i++)
{
gsl_matrix_set(kalmandata,0,i,gsl_matrix_get(kmanptr->xdata,i,0)); //initial value of state for kalman filter
}
//initialise output[0]
print2scr(m->C);
print2scr(xdata);
print2scr(ydata);
MatMul(m->C,xdata,ydata);
for(i=0;i<m->C->size1;i++)
gsl_matrix_set(msim->output,0,i,gsl_matrix_get(ydata,i,0));
/*end inital values in matrices*/


//simulation loop
 for(j=1;j<msim->Ndatapoints;j++)
{
//take values from u matrix n puts it in udata vector
for(i=0;i<m->B->size2;i++)
{
  gsl_matrix_set(udata,i,0,gsl_matrix_get(u,i,j));
  gsl_matrix_set(udataDis,i,0,gsl_matrix_get(u,i,j));
}





///xdot=A*x+B*u
Kalman_Step(kmanptr,ydata,udata);//Kalman_Step(Kalman_struc *Kalman,gsl_matrix *ymeas,gsl_matrix *u)

MatMul(m->A,xdata,tmpA);
MatMul(m->B,udataDis,tmpB);
MatAdd(tmpA,tmpB,xdata);
MatMul(m->C,xdata,ydata);

msim->time_vec[j]=j*Ts;
for(i=0;i<m->A->size1;i++)
{
gsl_matrix_set(msim->state,j,i,gsl_matrix_get(xdata,i,0));
gsl_matrix_set(kalmandata,j,i,gsl_matrix_get(kmanptr->xdata,i,0));
}

for(i=0;i<kmanptr->A->size1;i++)
{
gsl_matrix_set(kalmandata,j,i,gsl_matrix_get(kmanptr->xdata,i,0));
}

for(i=0;i<m->C->size1;i++)
{
gsl_matrix_set(msim->output,j,i,gsl_matrix_get(ydata,i,0));

}

}

print2File(kalmandata,msim->time_vec,"kalman.dat");

gsl_matrix_free(tmpA);
gsl_matrix_free(tmpB);
gsl_matrix_free(xdata);
gsl_matrix_free(udata);
}



