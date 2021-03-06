#include "kalman.h"

void Kalman_Init(Kalman_struc *Kalman,Model *m,gsl_matrix *Qcov, gsl_matrix *Rcov)
{
int Ns,Ny,Nu;
long Npoints; //number os states, output, inputsKalman=m;
Ns=m->A->size1;
Nu=m->B->size2;
Ny=m->C->size1;

Kalman->A=gsl_matrix_alloc(Ns,Ns);
Kalman->B=gsl_matrix_alloc(Ns,Nu);
Kalman->C=gsl_matrix_alloc(Ny,Ns);
Kalman->D=gsl_matrix_alloc(Ny,Nu);

gsl_matrix_memcpy(Kalman->A,m->A);
gsl_matrix_memcpy(Kalman->B,m->B);
gsl_matrix_memcpy(Kalman->C,m->C);
gsl_matrix_memcpy(Kalman->D,m->D);

Kalman->Ts=m->Ts;

Kalman->Q=gsl_matrix_alloc(Ns,Ns);
gsl_matrix_set_identity(Kalman->Q);
Kalman->R=gsl_matrix_alloc(Ny,Ny);
gsl_matrix_set_identity(Kalman->R);
Kalman->P=gsl_matrix_alloc(Ns,Ns);
gsl_matrix_set_all(Kalman->P,0);
Kalman->K=gsl_matrix_alloc(Ns,Ny);
gsl_matrix_set_all(Kalman->K,0);
Kalman->xdata=gsl_matrix_alloc(Ns,1);

gsl_matrix_memcpy(Kalman->xdata,m->X0);

gsl_matrix_memcpy(Kalman->Q,Qcov);
gsl_matrix_memcpy(Kalman->R,Rcov);
//gsl_matrix_set_identity(Kalman->P);
}

void Kalman_Step(Kalman_struc *Kalman,gsl_matrix *ymeas,gsl_matrix *u)
{

//initialisation
gsl_matrix *xprev=gsl_matrix_alloc(Kalman->xdata->size1,Kalman->xdata->size2);
gsl_matrix *tmpA=gsl_matrix_alloc(Kalman->xdata->size1,1);
gsl_matrix *tmpB=gsl_matrix_alloc(Kalman->xdata->size1,2);
gsl_matrix *tmpC=gsl_matrix_alloc(Kalman->A->size1,Kalman->A->size1);
gsl_matrix *tmpD=gsl_matrix_alloc(Kalman->C->size1,Kalman->C->size1);
gsl_matrix *idenNs=gsl_matrix_alloc(Kalman->A->size1,Kalman->A->size1);

//create NsxNs identity matrix
gsl_matrix_set_identity(idenNs);
gsl_matrix_memcpy(xprev,Kalman->xdata);
/*Time update Equations*/
//x_[k]=A*x[k-1]+B*u[k]  x_ denotes a priori estimate //tmpA,tmpB
//P_[k]=A*P[k-1]*A^T+Q[k-1] //tmpQ
MatAdd(MatMul2(Kalman->A,Kalman->xdata),MatMul2(Kalman->B,u),xprev);
MatAdd(MatMul2(MatMul2(Kalman->A,Kalman->P),MatTrans(Kalman->A)),Kalman->Q,Kalman->P);


/* printf("Kalman Model A");
print2scr(Kalman->A);
printf("u");
print2scr(u);
print2scr(Kalman->B); */
/*Calculation of Kalman Gain*/
//K[k]=P_[k]*C[k]^T*inv(C[k]*P_[k]*C[k]^T+R[k])
//printsizeMat(tmpD);
tmpD=MatMul2(MatMul2(Kalman->C,Kalman->P),MatTrans(Kalman->C));
MatAdd(tmpD,Kalman->R,tmpD);
Kalman->K=MatMul2(MatMul2(Kalman->P,MatTrans(Kalman->C)),MatInv2(tmpD));
//Kalman->K=MatMul2(MatMul2(Kalman->P,MatTrans(Kalman->C)),MatInv2(tmpD));
/*Measurement Update Equations*/
//x[k]=x_[k]+K[k]*(y[k]-C*x_[k])
//P[k]=(I-K[k]*C[k])*P_[k]*(I-K[k]*C[k])^T+K[k]*R[k]*K[k];
Kalman->xdata=MatAdd2(xprev,MatMul2(Kalman->K,MatSub2(ymeas,MatMul2(Kalman->C,xprev))));
//(I-k[k]*C[k])=tmpC
tmpC=MatSub2(MatMul2(Kalman->K,Kalman->C),idenNs);
Kalman->P=MatMul2(tmpC,MatMul2(Kalman->P,MatTrans(tmpC)));
Kalman->P=MatAdd2(Kalman->P,MatMul2(Kalman->K,MatMul2(Kalman->R,MatTrans(Kalman->K))));

gsl_matrix_free(xprev);
gsl_matrix_free(tmpA);
gsl_matrix_free(tmpB);
gsl_matrix_free(tmpC);
gsl_matrix_free(tmpD);
gsl_matrix_free(idenNs);

}

void Kalman_Init_d_o(Kalman_struc *Kalman,Model *m)
{
//Assuming input disturbance
int Ns,Ny,Nu,Nsd,Nyd,Nud; //number os states, output, inputsKalman=m;
int i,j,k;//counters
gsl_matrix *Ad;
gsl_matrix *Cd;
//Disturbance model
//xd[k+1]=Ad.xd[k] being constant disturbance Ad=eye
//yd=Cd.xd[k]
////output disturbance..so number of disturbance output=number of output
//constant disturbance on each channel
Nyd=m->C->size1;
////Nsd--independent we can have as many states as we want to model disturbance
Nsd=Nyd; //let say its a constant disturbance acting on each input
Nud=0; //for now
Ns=m->A->size1;
Ns=Ns+Nsd;
Nu=m->B->size2;
Ny=m->C->size1;

//
Ad=gsl_matrix_alloc(Nsd,Nsd);
Cd=gsl_matrix_alloc(Nyd,Nyd);
//
////Augmented Model becomes
////Au=[A 0;0 Ad] Bu=[B;0] Cu=[C Cd]
//

//
//
Kalman->A=gsl_matrix_alloc(Ns,Ns);
Kalman->B=gsl_matrix_alloc(Ns,Nu);
Kalman->C=gsl_matrix_alloc(Ny,Ns);
Kalman->D=gsl_matrix_alloc(Ny,Nu);


Kalman->statedata=gsl_matrix_alloc(Ns,m->Ndatapoints);



//
////set Ad and Cd matrix
gsl_matrix_set_identity(Ad);
gsl_matrix_set_identity(Cd);
//

gsl_matrix_set_all(Kalman->A,0);
gsl_matrix_set_all(Kalman->C,0);
k=0;
for(i=0;i<Ns;i++)
for(j=0;j<Ns;j++)
{
   if(((i<(Ns-Nsd))&&(j<(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(m->A,i,j));
    }
    else if(((i<(Ns-Nsd))&&(j>=(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,0);

    }
    else if (((i>=(Ns-Nsd))&&(j<(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,0);
    }
    else
    {
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(Ad,i-Ns+Nsd,j-Ns+Nsd));
  }

}
//bmatrix
for(i=0;i<Ns;i++)
for(j=0;j<Nu;j++)
{
    if(i<(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->B,i,j,gsl_matrix_get(m->B,i,j));
    }
    else if(i>=(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->B,i,j,0);

    }
}
print2scr(Kalman->B);
//cmatrix
for(i=0;i<Ny;i++)
for(j=0;j<Ns;j++)
{
    if(j<(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->C,i,j,gsl_matrix_get(m->C,i,j));
    }
    else
    {
        gsl_matrix_set(Kalman->C,i,j,gsl_matrix_get(Cd,i,j-Ns+Nsd));

    }
}
Kalman->Q=gsl_matrix_alloc(Ns,Ns);
Kalman->R=gsl_matrix_alloc(Ny,Ny);
gsl_matrix_set_identity(Kalman->Q);
gsl_matrix_set_identity(Kalman->R);
//gsl_matrix_set(Kalman->Q,5,5,10);
//gsl_matrix_set(Kalman->Q,6,6,100);
printf("kalman->Q");
print2scr(Kalman->Q);
printsizeMat(Kalman->C,"f");
//Kalman->Q=gsl_matrix_alloc(Ns,Ns);
//Kalman->R=gsl_matrix_alloc(Ny,Ny);
Kalman->P=gsl_matrix_alloc(Ns,Ns);
Kalman->K=gsl_matrix_alloc(Ns,Ny);
Kalman->xdata=gsl_matrix_alloc(Ns,1);
gsl_matrix_set_all(Kalman->xdata,0);

//Kalman->R=Rcov;

for(i=0;i<Kalman->statedata->size1;i++)
    gsl_matrix_set(Kalman->statedata,i,0,gsl_matrix_get(Kalman->xdata,i,0));

    Kalman->Ts=m->Ts;

}

void Kalman_Init_d(Kalman_struc *Kalman,Model *m)
{
//Assuming input disturbance
int Ns,Ny,Nu,Nsd,Nyd,Nud; //number os states, output, inputsKalman=m;
int i,j,k;//counters

gsl_matrix *Ad;
gsl_matrix *Cd;
gsl_matrix *tmpBCd;

//Disturbance model
//xd[k+1]=Ad.xd[k] being constant disturbance Ad=eye
//yd=Cd.xd[k]
////input disturbance..so number of disturbance output=number of input
Nyd=m->B->size2;
////Nsd--independent we can have as many states as we want to model disturbance
Nsd=Nyd; //let say its a constant disturbance acting on each input
Nud=0; //for now
Ns=m->A->size1;
Ns=Ns+Nsd;
Nu=m->B->size2;
Ny=m->C->size1;
Kalman->statedata=gsl_matrix_alloc(Ns,m->Ndatapoints);
//
Ad=gsl_matrix_alloc(Nsd,Nsd);
Cd=gsl_matrix_alloc(Nyd,Nyd);
//
////Augmented Model becomes
////Au=[A B*Cd;0 Ad] Bu=[B;0] Cu=[C 0]
tmpBCd=gsl_matrix_alloc(m->B->size1,Cd->size2);
//

//
//
Kalman->A=gsl_matrix_alloc(Ns,Ns);
Kalman->B=gsl_matrix_alloc(Ns,Nu);
Kalman->C=gsl_matrix_alloc(Ny,Ns);
Kalman->D=gsl_matrix_alloc(Ny,Nu);



//
////set Ad and Cd matrix
gsl_matrix_set_identity(Ad);
gsl_matrix_set_identity(Cd);
//
////Construction of augmented matrix
tmpBCd=MatMul2(m->B,Cd);
gsl_matrix_set_identity(tmpBCd);
gsl_matrix_set_all(Kalman->A,0);
gsl_matrix_set_all(Kalman->C,0);
k=0;
//print2scr(m->A->size1);
print2scr(tmpBCd);
print2scr(Kalman->C);
for(i=0;i<Ns;i++)
for(j=0;j<Ns;j++)
{
   if(((i<(Ns-Nsd))&&(j<(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(m->A,i,j));
    }
    else if(((i<(Ns-Nsd))&&(j>=(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(tmpBCd,i,j-Ns+Nsd));

    }
    else if (((i>=(Ns-Nsd))&&(j<(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,0);
    }
    else
    {
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(Ad,i-Ns+Nsd,j-Ns+Nsd));
  }

}

//bmatrix
for(i=0;i<Ns;i++)
for(j=0;j<Nu;j++)
{
    if(i<(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->B,i,j,gsl_matrix_get(m->B,i,j));
    }
    else if(i>=(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->B,i,j,0);

    }
}
print2scr(Kalman->B);
//cmatrix
for(i=0;i<Ny;i++)
for(j=0;j<Ns;j++)
{
    if(j<(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->C,i,j,gsl_matrix_get(m->C,i,j));
    }
    else
    {
        gsl_matrix_set(Kalman->C,i,j,0);

    }
}
Kalman->Q=gsl_matrix_alloc(Ns,Ns);
Kalman->R=gsl_matrix_alloc(Ny,Ny);
gsl_matrix_set_identity(Kalman->Q);
gsl_matrix_set_identity(Kalman->R);
printf("kalman->Q");
print2scr(Kalman->Q);
printsizeMat(Kalman->C,"f");
//Kalman->Q=gsl_matrix_alloc(Ns,Ns);
//Kalman->R=gsl_matrix_alloc(Ny,Ny);
Kalman->P=gsl_matrix_alloc(Ns,Ns);
Kalman->K=gsl_matrix_alloc(Ns,Ny);
Kalman->xdata=gsl_matrix_alloc(Ns,1);
gsl_matrix_set_all(Kalman->xdata,0);

for(i=0;i<Kalman->statedata->size1;i++)
    gsl_matrix_set(Kalman->statedata,i,0,gsl_matrix_get(Kalman->xdata,i,0));

    Kalman->Ts=m->Ts;

}

void Kalman_Step2(Kalman_struc *Kalman,gsl_matrix *ymeas,gsl_matrix *u)
{
int i,j,Ns,Nu,Ny;

Ns=Kalman->A->size1;
Nu=Kalman->B->size2;
Ny=Kalman->C->size1;
//initialisation
gsl_matrix *xprev=gsl_matrix_alloc(Kalman->xdata->size1,Kalman->xdata->size2);
gsl_matrix *tmpA=gsl_matrix_alloc(Kalman->xdata->size1,1);
gsl_matrix *tmpB=gsl_matrix_alloc(Kalman->xdata->size1,2);
gsl_matrix *tmpC=gsl_matrix_alloc(Kalman->A->size1,Kalman->A->size1);
gsl_matrix *tmpD=gsl_matrix_alloc(Kalman->C->size1,Kalman->C->size1);
gsl_matrix *idenNs=gsl_matrix_alloc(Kalman->A->size1,Kalman->A->size1);

double xprevdouble[Ns];

//create NsxNs identity matrix
gsl_matrix_set_identity(idenNs);
gsl_matrix_memcpy(xprev,Kalman->xdata);

/*Time update Equations*/
//x_[k]=A*x[k-1]+B*u[k]  x_ denotes a priori estimate //tmpA,tmpB
//P_[k]=A*P[k-1]*A^T+Q[k-1] //tmpQ

///element by element multiplication
for(i=0;i<Ns;i++)
{
    xprevdouble[i]=0.0;
    for(j=0;j<Ns;j++)
        xprevdouble[i]=xprevdouble[i]+gsl_matrix_get(Kalman->A,i,j)*gsl_matrix_get(Kalman->xdata,j,0);
}

for(i=0;i<Ns;i++)
{
    for(j=0;j<Nu;j++)
        xprevdouble[i]=xprevdouble[i]+gsl_matrix_get(Kalman->B,i,j)*gsl_matrix_get(u,j,0);
}
assign_Mat(xprev,xprevdouble);

//MatAdd(MatMul2(Kalman->A,Kalman->xdata),MatMul2(Kalman->B,u),xprev);
MatAdd(MatMul2(MatMul2(Kalman->A,Kalman->P),MatTrans(Kalman->A)),Kalman->Q,Kalman->P);

/* printf("Kalman Model A");
print2scr(Kalman->A);
printf("u");
print2scr(u);
print2scr(Kalman->B); */
/*Calculation of Kalman Gain*/
//K[k]=P_[k]*C[k]^T*inv(C[k]*P_[k]*C[k]^T+R[k])
//printsizeMat(tmpD);
tmpD=MatMul2(MatMul2(Kalman->C,Kalman->P),MatTrans(Kalman->C));
MatAdd(tmpD,Kalman->R,tmpD);
Kalman->K=MatMul2(MatMul2(Kalman->P,MatTrans(Kalman->C)),MatInv2(tmpD));

//Kalman->K=MatMul2(MatMul2(Kalman->P,MatTrans(Kalman->C)),MatInv2(tmpD));
/*Measurement Update Equations*/
//x[k]=x_[k]+K[k]*(y[k]-C*x_[k])
//P[k]=(I-K[k]*C[k])*P_[k]*(I-K[k]*C[k])^T+K[k]*R[k]*K[k];
Kalman->xdata=MatAdd2(xprev,MatMul2(Kalman->K,MatSub2(ymeas,MatMul2(Kalman->C,xprev))));
//(I-k[k]*C[k])=tmpC
tmpC=MatSub2(MatMul2(Kalman->K,Kalman->C),idenNs);
Kalman->P=MatMul2(tmpC,MatMul2(Kalman->P,MatTrans(tmpC)));
Kalman->P=MatAdd2(Kalman->P,MatMul2(Kalman->K,MatMul2(Kalman->R,MatTrans(Kalman->K))));

//for(i=0;i<Kalman->statedata->size1;i++)
//gsl_matrix_set(Kalman->statedata,i,step_value,gsl_matrix_get(Kalman->xdata,i,0));

gsl_matrix_free(xprev);
gsl_matrix_free(tmpA);
gsl_matrix_free(tmpB);
gsl_matrix_free(tmpC);
gsl_matrix_free(tmpD);
gsl_matrix_free(idenNs);

}

void Kalman_Step2Ame(Kalman_struc *Kalman,gsl_matrix *ymeas,gsl_matrix *u,int step_value)
{
int i,j;
//initialisation
gsl_matrix *xprev=gsl_matrix_alloc(Kalman->xdata->size1,Kalman->xdata->size2);
gsl_matrix *tmpA=gsl_matrix_alloc(Kalman->xdata->size1,1);
gsl_matrix *tmpB=gsl_matrix_alloc(Kalman->xdata->size1,2);
gsl_matrix *tmpC=gsl_matrix_alloc(Kalman->A->size1,Kalman->A->size1);
gsl_matrix *tmpD=gsl_matrix_alloc(Kalman->C->size1,Kalman->C->size1);
gsl_matrix *idenNs=gsl_matrix_alloc(Kalman->A->size1,Kalman->A->size1);

//create NsxNs identity matrix
gsl_matrix_set_identity(idenNs);
gsl_matrix_memcpy(xprev,Kalman->xdata);
/*Time update Equations*/
//x_[k]=A*x[k-1]+B*u[k]  x_ denotes a priori estimate //tmpA,tmpB
//P_[k]=A*P[k-1]*A^T+Q[k-1] //tmpQ
MatAdd(MatMul2(Kalman->A,Kalman->xdata),MatMul2(Kalman->B,u),xprev);
MatAdd(MatMul2(MatMul2(Kalman->A,Kalman->P),MatTrans(Kalman->A)),Kalman->Q,Kalman->P);
/* printf("Kalman Model A");
print2scr(Kalman->A);
printf("u");
print2scr(u);
print2scr(Kalman->B); */
/*Calculation of Kalman Gain*/
//K[k]=P_[k]*C[k]^T*inv(C[k]*P_[k]*C[k]^T+R[k])
//printsizeMat(tmpD);
tmpD=MatMul2(MatMul2(Kalman->C,Kalman->P),MatTrans(Kalman->C));
MatAdd(tmpD,Kalman->R,tmpD);
Kalman->K=MatMul2(MatMul2(Kalman->P,MatTrans(Kalman->C)),MatInv2(tmpD));
//Kalman->K=MatMul2(MatMul2(Kalman->P,MatTrans(Kalman->C)),MatInv2(tmpD));
/*Measurement Update Equations*/
//x[k]=x_[k]+K[k]*(y[k]-C*x_[k])
//P[k]=(I-K[k]*C[k])*P_[k]*(I-K[k]*C[k])^T+K[k]*R[k]*K[k];
Kalman->xdata=MatAdd2(xprev,MatMul2(Kalman->K,MatSub2(ymeas,MatMul2(Kalman->C,xprev))));
//(I-k[k]*C[k])=tmpC
tmpC=MatSub2(MatMul2(Kalman->K,Kalman->C),idenNs);
Kalman->P=MatMul2(tmpC,MatMul2(Kalman->P,MatTrans(tmpC)));
Kalman->P=MatAdd2(Kalman->P,MatMul2(Kalman->K,MatMul2(Kalman->R,MatTrans(Kalman->K))));



gsl_matrix_free(xprev);
gsl_matrix_free(tmpA);
gsl_matrix_free(tmpB);
gsl_matrix_free(tmpC);
gsl_matrix_free(tmpD);
gsl_matrix_free(idenNs);


}
void printKalmanData(Kalman_struc *m,int s,char *filename)
{
    //s==1 prints states //s==0 prints output
   FILE *fp;
    int rows;
    int cols;
    int i,j; //counters for loops
    rows=m->statedata->size1;
    cols=m->statedata->size2;
    fp=fopen(filename,"w");
    for(j=0;j<cols;j++)
    {
        if(j!=0)
        fprintf(fp,"\n");

        fprintf(fp,"%f ",j*m->Ts);
        for(i=0;i<rows;i++)
        {
            fprintf(fp,"%f ",gsl_matrix_get(m->statedata,i,j));
        }
    }
    fclose(fp);


}

void Kalman_Init_d_servo(Kalman_struc *Kalman,Model *m)
{
//Assuming input disturbance
int Ns,Ny,Nu,Nsd,Nyd,Nud; //number os states, output, inputsKalman=m;
int i,j,k;//counters

gsl_matrix *Ad;
gsl_matrix *Cd;
gsl_matrix *tmpBCd;
double tempBCDfromfile[]={0,0,0,-0.9549,0};
//Disturbance model
//xd[k+1]=Ad.xd[k] being constant disturbance Ad=eye
//yd=Cd.xd[k]
////input disturbance..so number of disturbance output=number of input
Nyd=m->B->size2;
////Nsd--independent we can have as many states as we want to model disturbance
Nsd=Nyd; //let say its a constant disturbance acting on each input
Nud=0; //for now
Ns=m->A->size1;
Ns=Ns+Nsd;
Nu=m->B->size2;
Ny=m->C->size1;
Kalman->statedata=gsl_matrix_alloc(Ns,m->Ndatapoints);
//
Ad=gsl_matrix_alloc(Nsd,Nsd);
Cd=gsl_matrix_alloc(Nyd,Nyd);
//
////Augmented Model becomes
////Au=[A B*Cd;0 Ad] Bu=[B;0] Cu=[C 0]
tmpBCd=gsl_matrix_alloc(m->B->size1,Cd->size2);
//

//
//
Kalman->A=gsl_matrix_alloc(Ns,Ns);
Kalman->B=gsl_matrix_alloc(Ns,Nu);
Kalman->C=gsl_matrix_alloc(Ny,Ns);
Kalman->D=gsl_matrix_alloc(Ny,Nu);



//
////set Ad and Cd matrix
gsl_matrix_set_identity(Ad);
gsl_matrix_set_identity(Cd);
//
////Construction of augmented matrix
//tmpBCd=MatMul2(m->B,Cd);
for(i=0;i<tmpBCd->size1;i++)
    gsl_matrix_set(tmpBCd,i,0,m->Ts*tempBCDfromfile[i]);
//gsl_matrix_set_identity(tmpBCd);
gsl_matrix_set_all(Kalman->A,0);
gsl_matrix_set_all(Kalman->C,0);
k=0;
//print2scr(m->A->size1);
print2scr(tmpBCd);
print2scr(Kalman->C);
for(i=0;i<Ns;i++)
for(j=0;j<Ns;j++)
{
   if(((i<(Ns-Nsd))&&(j<(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(m->A,i,j));
    }
    else if(((i<(Ns-Nsd))&&(j>=(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(tmpBCd,i,j-Ns+Nsd));

    }
    else if (((i>=(Ns-Nsd))&&(j<(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,0);
    }
    else
    {
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(Ad,i-Ns+Nsd,j-Ns+Nsd));
  }

}

//bmatrix
for(i=0;i<Ns;i++)
for(j=0;j<Nu;j++)
{
    if(i<(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->B,i,j,gsl_matrix_get(m->B,i,j));
    }
    else if(i>=(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->B,i,j,0);

    }
}
print2scr(Kalman->B);
//cmatrix
for(i=0;i<Ny;i++)
for(j=0;j<Ns;j++)
{
    if(j<(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->C,i,j,gsl_matrix_get(m->C,i,j));
    }
    else
    {
        gsl_matrix_set(Kalman->C,i,j,0);

    }
}
Kalman->Q=gsl_matrix_alloc(Ns,Ns);
Kalman->R=gsl_matrix_alloc(Ny,Ny);
gsl_matrix_set_identity(Kalman->Q);
gsl_matrix_set_identity(Kalman->R);

//Kalman->Q=gsl_matrix_alloc(Ns,Ns);
//Kalman->R=gsl_matrix_alloc(Ny,Ny);
Kalman->P=gsl_matrix_alloc(Ns,Ns);
Kalman->K=gsl_matrix_alloc(Ns,Ny);
Kalman->xdata=gsl_matrix_alloc(Ns,1);
gsl_matrix_set_all(Kalman->xdata,0);

for(i=0;i<Kalman->statedata->size1;i++)
    gsl_matrix_set(Kalman->statedata,i,0,gsl_matrix_get(Kalman->xdata,i,0));

    Kalman->Ts=m->Ts;

}
void Kalman_Init_Ame_servo(Kalman_struc *Kalman,Model *m,double Ts)
{
//Assuming input disturbance
int Ns,Ny,Nu,Nsd,Nyd,Nud; //number os states, output, inputsKalman=m;
int i,j,k;//counters
gsl_matrix *Ad;
gsl_matrix *Cd;
gsl_matrix *tmpBCd;
double tempBCDfromfile[]={0,0,0,-0.9549,0};
//Disturbance model
//xd[k+1]=Ad.xd[k] being constant disturbance Ad=eye
//yd=Cd.xd[k]
////input disturbance..so number of disturbance output=number of input
Kalman->Ts=Ts;
Nyd=m->B->size2;
////Nsd--independent we can have as many states as we want to model disturbance
Nsd=Nyd; //let say its a constant disturbance acting on each input
Nud=0; //for now
Ns=m->A->size1;
Ns=Ns+Nsd;
Nu=m->B->size2;
Ny=m->C->size1;

//
Ad=gsl_matrix_alloc(Nsd,Nsd);
Cd=gsl_matrix_alloc(Nyd,Nyd);
//
////Augmented Model becomes
////Au=[A B*Cd;0 Ad] Bu=[B;0] Cu=[C 0]
tmpBCd=gsl_matrix_alloc(m->B->size1,Cd->size2);
//

//
//
Kalman->A=gsl_matrix_alloc(Ns,Ns);
Kalman->B=gsl_matrix_alloc(Ns,Nu);
Kalman->C=gsl_matrix_alloc(Ny,Ns);
Kalman->D=gsl_matrix_alloc(Ny,Nu);



//
////set Ad and Cd matrix
gsl_matrix_set_identity(Ad);
gsl_matrix_set_identity(Cd);
//
////Construction of augmented matrix
//tmpBCd=MatMul2(m->B,Cd);
for(i=0;i<tmpBCd->size1;i++)
    gsl_matrix_set(tmpBCd,i,0,Kalman->Ts*tempBCDfromfile[i]);
//gsl_matrix_set_identity(tmpBCd);
gsl_matrix_set_all(Kalman->A,0);
gsl_matrix_set_all(Kalman->C,0);
k=0;
//print2scr(m->A->size1);

for(i=0;i<Ns;i++)
for(j=0;j<Ns;j++)
{
   if(((i<(Ns-Nsd))&&(j<(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(m->A,i,j));
    }
    else if(((i<(Ns-Nsd))&&(j>=(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(tmpBCd,i,j-Ns+Nsd));

    }
    else if (((i>=(Ns-Nsd))&&(j<(Ns-Nsd))))
    {
        gsl_matrix_set(Kalman->A,i,j,0);
    }
    else
    {
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(Ad,i-Ns+Nsd,j-Ns+Nsd));
  }

}

//bmatrix
for(i=0;i<Ns;i++)
for(j=0;j<Nu;j++)
{
    if(i<(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->B,i,j,gsl_matrix_get(m->B,i,j));
    }
    else if(i>=(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->B,i,j,0);

    }
}

//cmatrix
for(i=0;i<Ny;i++)
for(j=0;j<Ns;j++)
{
    if(j<(Ns-Nsd))
    {
        gsl_matrix_set(Kalman->C,i,j,gsl_matrix_get(m->C,i,j));
    }
    else
    {
        gsl_matrix_set(Kalman->C,i,j,0);

    }
}
Kalman->Q=gsl_matrix_alloc(Ns,Ns);
Kalman->R=gsl_matrix_alloc(Ny,Ny);
gsl_matrix_set_identity(Kalman->Q);
gsl_matrix_set_identity(Kalman->R);

//Kalman->Q=gsl_matrix_alloc(Ns,Ns);
//Kalman->R=gsl_matrix_alloc(Ny,Ny);
Kalman->P=gsl_matrix_alloc(Ns,Ns);
Kalman->K=gsl_matrix_alloc(Ns,Ny);
Kalman->xdata=gsl_matrix_alloc(Ns,1);
gsl_matrix_set_all(Kalman->xdata,0);

gsl_matrix_free(Ad);
gsl_matrix_free(Cd);
gsl_matrix_free(tmpBCd);

}


void Kalman_Init_Dist(Kalman_struc *Kalman,Model *m,double Ts,double *Bd,double *Dd,int Ndi,int Ndo,double *Qkal,double *Rkal)
{
//Assuming input disturbance
int Ns,Ny,Nu,Nsd,Nyd,Nud,Ntotal; //number os states, output, inputsKalman=m;
int i,j,k;//counters
gsl_matrix *Ad;
gsl_matrix *Cd;
gsl_matrix *tmpBCd;
Kalman->Ts=Ts;

Ns=m->A->size1;
Nu=m->B->size2;
Ny=m->C->size1;
/// We consider system of the general form
/// x[k+1]=Ax[k]+Bu[k]+Bd.w
///Disturbance model
///xd[k+1]=Ad.xd[k] being constant disturbance Ad=eye
///w=Cd.xd[k]
///We are assuming a constant disturbance, so Ad=Cd=I (Ndi x Ndi)
///Nsd--independent we can have as many states as we want to model disturbance

Ntotal=Ns+Ndi+Ndo; ///Total number of states


///allocate matrices
Kalman->A=gsl_matrix_alloc(Ntotal,Ntotal);
Kalman->B=gsl_matrix_alloc(Ntotal,Nu);
Kalman->C=gsl_matrix_alloc(Ny,Ntotal);
Kalman->D=gsl_matrix_alloc(Ny,Nu);

///populate A matrix
///    [A Bd 0]
/// A= [0  Adi 0]
///    [0  0 Ado]
/// Adi and Ado (as of now are taken to be identity Matrices as disturbance is assumed constant.
/// Adi=I, Ado=I
/// little trick---after the row [ A Bd 0], we will just assign a diagonal matrice.
/// i.e for i<Ns we will assign a diagonal matrix

gsl_matrix_set_all(Kalman->A,0); ///set all entries to zero

///FIRST ROW ASSIGNMENT
for(i=0;i<Ns;i++)
    for(j=0;j<Ns+Ndi;j++)
{
    if(j<Ns)
        gsl_matrix_set(Kalman->A,i,j,gsl_matrix_get(m->A,i,j));
    if(j>=Ns)
        gsl_matrix_set(Kalman->A,i,j,Bd[i+(j-Ns)*Ns]);
}

///DIAGONAL ELEMENTS ASSIGNMENT; ///Here is where the Ad matrix should be set in case disturbance model is not a constant one
for(i=Ns;i<Ntotal;i++)
    gsl_matrix_set(Kalman->A,i,i,1.0);

/// B Matrix
gsl_matrix_set_all(Kalman->B,0);

for(i=0;i<Ns;i++)
    for(j=0;j<Nu;j++)
        gsl_matrix_set(Kalman->B,i,j,gsl_matrix_get(m->B,i,j));

/// C matrix
/// C=[C 0 Dd]
/// In theory we should have a value of Dd but here we take Dd=I
gsl_matrix_set_all(Kalman->C,0);

for(i=0;i<Ny;i++)
    for(j=0;j<Ns;j++)
        gsl_matrix_set(Kalman->C,i,j,gsl_matrix_get(m->C,i,j));

///Setting output disturbance
     for(i=0;i<10;i++)
        printf("Dd[%d]=%f\n",i,Dd[i]);
    for(i=0;i<Ny;i++)
        for(j=Ns+Ndi;j<Ntotal;j++)
        {
            gsl_matrix_set(Kalman->C,i,j,Dd[i+(j-Ns-Ndi)*Ny]);
            printf("Dd[%d]=%f\n",i+(j-Ns-Ndi)*Ny,Dd[i+(j-Ns-Ndi)*Ny]);

        }


/// D matrix
gsl_matrix_memcpy(Kalman->D,m->D);


/// Q,R,P and K matrices

Kalman->Q=gsl_matrix_alloc(Ntotal,Ntotal);
gsl_matrix_set_identity(Kalman->Q);
Kalman->R=gsl_matrix_alloc(Ny,Ny);
gsl_matrix_set_identity(Kalman->R);
Kalman->P=gsl_matrix_alloc(Ntotal,Ntotal);
gsl_matrix_set_all(Kalman->P,0);
Kalman->K=gsl_matrix_alloc(Ntotal,Ny);
gsl_matrix_set_all(Kalman->K,0);

assign_Diag(Kalman->Q,Qkal);
assign_Diag(Kalman->R,Rkal);





Kalman->xdata=gsl_matrix_alloc(Ntotal,1);
gsl_matrix_set_all(Kalman->xdata,0);
    for(i=0;i<Ns;i++)
    {
            gsl_matrix_set(Kalman->xdata,i,0,gsl_matrix_get(m->X0,i,0));
    }


Kalman->initflag=1;





}
int Controllability(Kalman_struc *Kalman)
{
  int Ns,Nu,Ny,i,j,k;
    gsl_matrix *Ctrb;
    gsl_matrix *tempmat;
    int rankmat;
    double tempval;

    Ns=Kalman->A->size1;
    Nu=Kalman->B->size2;

    Ctrb=gsl_matrix_alloc(Ns,Ns*Nu);
    tempmat=gsl_matrix_alloc(Ns,Nu);


    for(k=0;k<Ns;k++)
    {
        tempmat=MatMul2(MatMulrec(Kalman->A,k),Kalman->B);
        for(i=0;i<Ns;i++)
        {
            for(j=0;j<Nu;j++)
            {
                gsl_matrix_set(Ctrb,i,j+k*Nu,gsl_matrix_get(tempmat,i,j));
            }
        }
    }

      printf("Controllability matrix");
    print2scr(Ctrb);
    rankmat=matRank(Ctrb,1.0e-6);
    printf("Rank of Matrix is:%d\n",rankmat);

    gsl_matrix_free(tempmat);
    gsl_matrix_free(Ctrb);

    return(rankmat);
}
int Observability(Kalman_struc *Kalman)
{
    int Ns,Nu,Ny,i,j,k;
    gsl_matrix *Obsv;
    gsl_matrix *tempmat;
    int rankmat;
    double tempval;

    Ns=Kalman->A->size1;
    Ny=Kalman->C->size1;

   Obsv=gsl_matrix_alloc(Ns*Ny,Ns);
    tempmat=gsl_matrix_alloc(Ny,Ns);



    for(k=0;k<Ns;k++)
    {
        tempmat=MatMul2(Kalman->C,MatMulrec(Kalman->A,k));

        for(i=0;i<Ny;i++)
        {
            for(j=0;j<Ns;j++)
            {
                gsl_matrix_set(Obsv,i+k*Ny,j,gsl_matrix_get(tempmat,i,j));
            }
        }
    }

    printf("Observability matrix");
    print2scr(Obsv);
    rankmat=matRank(Obsv,1.0e-6);
    printf("Rank of Matrix is:%d\n",rankmat);
    gsl_matrix_free(tempmat);
    gsl_matrix_free(Obsv);

    return rankmat;

}

///findSVD does the SVD decomposition of the matrice A and return the diagonal elements in a vector
///if A is rectangular M x N matrice. S is an N x 1 matrice.

void findSVD(gsl_matrix *A,gsl_vector *S)
{
 int N,M,i;
 M=A->size1;
 N=A->size2;
 int error;

 gsl_matrix *U=gsl_matrix_alloc(M,N);
 gsl_matrix_memcpy(U,A);
gsl_vector *work=gsl_vector_alloc(N);
 gsl_matrix *V=gsl_matrix_alloc(N,N);

 error=gsl_linalg_SV_decomp(U,V,S,work);


/* printf("SVD decomposition\n");
 printf("U:");
 print2scr(U);
 printf("V:");
 print2scr(V);

printf("Singular Values S:\n");
 for(i=0;i<N;i++)
    printf("S[%d]:%f\n",i,gsl_vector_get(S,i));*/

gsl_matrix_free(U);
gsl_matrix_free(V);
gsl_vector_free(work);

}
int matRank(gsl_matrix *R,double tol)
{
    gsl_vector *S=gsl_vector_alloc(R->size2);
    int counter=0,i;

    findSVD(R,S);

    for(i=0;i<S->size;i++)
    {
        if (gsl_vector_get(S,i)>tol)
            counter=counter+1;
    }

    return counter;
}


