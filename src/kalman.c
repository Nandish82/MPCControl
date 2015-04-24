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


void Kalman_Init_delta(Kalman_struc *Kalman,Model *m,double Ts,double *Bd,int Ndi,int Ndo)
{
//Assuming input disturbance
int Ns,Ny,Nu,Nsd,Nyd,Nud,Ntotal; //number os states, output, inputsKalman=m;
int i,j,k;//counters
gsl_matrix *Ad;
gsl_matrix *Cd;
gsl_matrix *tmpBCd;
Kalman->Ts=Ts;
/// We consider system of the general form
/// x[k+1]=Ax[k]+Bu[k]+Bd.w
///Disturbance model
///xd[k+1]=Ad.xd[k] being constant disturbance Ad=eye
///w=Cd.xd[k]
///We are assuming a constant disturbance, so Ad=Cd=I (Ndi x Ndi)
Nsd=Ndi;
///Nsd--independent we can have as many states as we want to model disturbance
Nyd=Nsd; ///let say its a constant disturbance acting on each input
Nud=0; ///for now
Ns=m->A->size1;
Ns=Ns+Nsd; ///augmented state vector
Nu=m->B->size2;
Ny=m->C->size1;


double tempBCDfromfile[]={0,0,0,0,0};
//
Ad=gsl_matrix_alloc(Nsd,Nsd);
Cd=gsl_matrix_alloc(Nyd,Nyd);
//
/// Augmented Model becomes
/// Au=[A B*Cd;0 Ad] Bu=[B;0] Cu=[C 0]
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
void Controllability(Kalman_struc *Kalman)
{
    #ifdef AMESIM
        ameprintf(stderr,"");
    #endif // AMESIM
}




