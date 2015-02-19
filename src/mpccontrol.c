#include "mpccontrol.h"

/** \brief
 *  This function takes a pointer to an mpc structure, a pointer to a model and and MPCType.
 *  The MPCType determines whether the formulation is normal or delta formulation. It creates
    appropriate A,B,C,D matrices which it stores in the MPC struct.
 * \param mpcptr Pointer to a MPC structure
 * \param modelptr Pointer to an already discretized model
 * \return
 *
 */


void InitMPCType(structMPC *mpcptr,Model *modelptr,MPCType type)
{
    int Ns,Nu,Ny; ///number of states, inputs and outputs
    int i,j;

    /**Allocating the sizes */
    Ns=modelptr->A->size1;
    Nu=modelptr->B->size2;
    Ny=modelptr->C->size1;

    printf("Ns:%d,Nu:%d,Ny:%d",Ns,Nu,Ny);
    if(type==NORMAL)
    {
        mpcptr->A=gsl_matrix_alloc(Ns,Ns);
        mpcptr->B=gsl_matrix_alloc(Ns,Nu);
        mpcptr->C=gsl_matrix_alloc(Ny,Ns);
        mpcptr->D=gsl_matrix_alloc(Ny,Nu);

        gsl_matrix_memcpy(mpcptr->A,modelptr->A);
        gsl_matrix_memcpy(mpcptr->B,modelptr->B);
        gsl_matrix_memcpy(mpcptr->C,modelptr->C);
        gsl_matrix_memcpy(mpcptr->D,modelptr->D);
        mpcptr->type=NORMAL;
    }
    else if(type==DELTA)
    {
        mpcptr->type=DELTA;
        mpcptr->A=gsl_matrix_alloc(Ns+Nu,Ns+Nu);
        mpcptr->B=gsl_matrix_alloc(Ns+Nu,Nu);
        mpcptr->C=gsl_matrix_alloc(Ny,Ns+Nu);
        mpcptr->D=gsl_matrix_alloc(Ny,Nu);

        gsl_matrix_set_all(mpcptr->A,0);
        gsl_matrix_set_all(mpcptr->B,0);
        gsl_matrix_set_all(mpcptr->C,0);
        gsl_matrix_set_all(mpcptr->D,0);

        ///A_delta
        for(i=0;i<(Ns+Nu);i++)
            for(j=0;j<(Ns+Nu);j++)
        {
            if((i<Ns)&&(j<Ns))
                gsl_matrix_set(mpcptr->A,i,j,gsl_matrix_get(modelptr->A,i,j));
            else if((i<Ns)&&(j>=Ns))
                gsl_matrix_set(mpcptr->A,i,j,gsl_matrix_get(modelptr->B,i,j-Ns));
             else if((i>=Ns)&&(j<Ns))
                gsl_matrix_set(mpcptr->A,i,j,0);
             else if((i>=Ns)&&(j>=Ns))
                if(i==j)
                    gsl_matrix_set(mpcptr->A,i,j,1); ///careful here its an identity matrix
        }
        ///B_delta
         for(i=0;i<(Ns+Nu);i++)
            for(j=0;j<(Nu);j++)
        {
            if((i<Ns)&&(j<Nu))
                gsl_matrix_set(mpcptr->B,i,j,gsl_matrix_get(modelptr->B,i,j));
            else if((i>=Ns)&&(j<Nu))
                    if((i-Ns)==j)
                        gsl_matrix_set(mpcptr->B,i,j,1);
        }
        ///Cdelta
         for(i=0;i<(Ny);i++)
            for(j=0;j<(Ns+Nu);j++)
        {
            if((i<Ny)&&(j<Ns))
                gsl_matrix_set(mpcptr->C,i,j,gsl_matrix_get(modelptr->C,i,j));
            else if((i<Ny)&&(j>=Ns))
                        gsl_matrix_set(mpcptr->C,i,j,gsl_matrix_get(modelptr->D,i,j-Ns));
        }
        ///D_delta
        // gsl_matrix_memcpy(mpcptr->D,modelptr->D);

    }
}
/** \brief
        Computes the steady state matrix i.e
        [A-I B;
         Cref 0]
         This is done so that in the calculation
         algorithm less time is consumed
 *
 * \param takes pointer to initialised mpc structure and stores matrix in mpcptr->Steady State
 * \param
 * \return
 *
 */

void InitSteadyState(structMPC *mpcptr,gsl_matrix *Cref)
{
    int Ns,Nu,Ny,i,j;
    if(mpcptr->type==DELTA)
    {
        Nu=mpcptr->B->size2;
        Ns=mpcptr->A->size1-Nu;
        Ny=Cref->size1; ///Ny should be equal to the reference being tracked. To have a square matrix Ny[cref]=Nu

        mpcptr->SteadyState=gsl_matrix_alloc(Ns+Ny,Ns+Nu);
        for(i=0;i<Ns+Ny;i++)
            for(j=0;j<Ns+Nu;j++)
                if((i<Ns)&&(j<Ns))
                    gsl_matrix_set(mpcptr->SteadyState,i,j,gsl_matrix_get(mpcptr->A,i,j)-1);
                else if((i<Ns)&&(j>=Ns))
                    gsl_matrix_set(mpcptr->SteadyState,i,j,gsl_matrix_get(mpcptr->B,i,j-Ns));
                else if((i>=Ns)&&(j<Ns))
                    gsl_matrix_set(mpcptr->SteadyState,i,j,gsl_matrix_get(Cref,i-Ns,j));
                 else if((i>=Ns)&&(j>=Ns))
                    gsl_matrix_set(mpcptr->SteadyState,i,j,0);
    }
    else
    {
        Nu=mpcptr->B->size2;
        Ns=mpcptr->A->size1;
        Ny=Cref->size1;

        mpcptr->SteadyState=gsl_matrix_alloc(Ns+Ny,Ns+Nu);
        for(i=0;i<Ns+Ny;i++)
            for(j=0;j<Ns+Nu;j++)
                if((i<Ns)&&(j<Ns))
                {
                    if(i==j)
                        gsl_matrix_set(mpcptr->SteadyState,i,j,gsl_matrix_get(mpcptr->A,i,j)-1);
                    else
                     gsl_matrix_set(mpcptr->SteadyState,i,j,gsl_matrix_get(mpcptr->A,i,j));
                }

                else if((i<Ns)&&(j>=Ns))
                    gsl_matrix_set(mpcptr->SteadyState,i,j,gsl_matrix_get(mpcptr->B,i,j-Ns));
                else if((i>=Ns)&&(j<Ns))
                    gsl_matrix_set(mpcptr->SteadyState,i,j,gsl_matrix_get(Cref,i-Ns,j));
                 else if((i>=Ns)&&(j>=Ns))
                    gsl_matrix_set(mpcptr->SteadyState,i,j,0);
    }
}

void InitMPC(structMPC *mpcptr,int cHorizon,gsl_matrix *Q,gsl_matrix *P,gsl_matrix *R,double *lb,double *ub,double *lbA,double *ubA)
{

/**
This function constructs the H,F,G matrix given the prediction horizon and a model.
Constraints will be added later. These matrices are stored in a MPC structure defined
in the corresponding header file (mpccontrol.h)
The H matrix is given by
H=CM^T*Q*CM+R
The F matrix is given by
F=C^T*Q*M
and the G matrix is given by
G=M^T*Q*M+Q

M=[A A^2 A^3 ... A^(N-1)]^T

CM=[B 0 0 0 ...;
    AB B 0 0...;
    A^2B AB B 0;]

Q is the weighting on the states. The last (Ns,Nu) elements are the terminal weight (P).
R is the weighting on the input.
*/

mpcptr->Ts=mpcptr->Ts;
mpcptr->contHor=cHorizon;
int Ns,Nu,Ny,i,j,k,l;

Ns=mpcptr->A->size1;
Nu=mpcptr->B->size2;
Ny=mpcptr->C->size1;

/**
Assigning Q and R matrices*/
mpcptr->Q=gsl_matrix_alloc(Ns,Ns);
assignMatMat(mpcptr->Q,Q);
mpcptr->R=gsl_matrix_alloc(Nu,Nu);
assignMatMat(mpcptr->R,R);

/**
Construction of CM matrix
*/

mpcptr->CM=gsl_matrix_alloc(Ns*cHorizon,Nu*cHorizon);
gsl_matrix_set_all(mpcptr->CM,0);
gsl_matrix *tempCM=gsl_matrix_alloc(Ns,Nu);
printf("I am here befor CM");
for(i=0;i<cHorizon;i++)
    for(j=0;j<cHorizon;j++)
{
      if(i>=j)
    {
        tempCM=MatMul2(MatMulrec(mpcptr->A,i-j),mpcptr->B);//computes A^(n-1)*B
        for(k=0;k<Ns;k++)
            for(l=0;l<Nu;l++)
            gsl_matrix_set(mpcptr->CM,i*Ns+k,j*Nu+l,gsl_matrix_get(tempCM,k,l));
    }
}

/** End of Construction of CM Matrix*/
printf("I am here after CM");
printf("mpcptr->CM");
print2scr(mpcptr->CM);
/**Construction of the Q matrix
We will set the diagonal elements to 1.
A better way is to get the elements from AMESIM
*/
gsl_matrix *Qbar=gsl_matrix_alloc(cHorizon*Ns,cHorizon*Ns);
gsl_matrix_set_all(Qbar,0);
for(i=0;i<cHorizon;i++)
    for(j=0;j<cHorizon;j++)
{

      if(i==j)
    {

        for(k=0;k<Ns;k++)
            for(l=0;l<Ns;l++)
            {
                if(i<cHorizon-1)
                gsl_matrix_set(Qbar,i*Ns+k,j*Ns+l,gsl_matrix_get(Q,k,l));
                if(i==cHorizon-1)
                gsl_matrix_set(Qbar,i*Ns+k,j*Ns+l,gsl_matrix_get(P,k,l));

            }

    }
}
printf("I am here after Q");
/**Constructioin of the R matrix
Again we set everything to 1
*/
gsl_matrix *Rbar=gsl_matrix_alloc(cHorizon*Nu,cHorizon*Nu);
gsl_matrix_set_all(Rbar,0);
for(i=0;i<cHorizon;i++)
{
     for(k=0;k<Nu;k++)
            for(l=0;l<Nu;l++)
                gsl_matrix_set(Rbar,i*Nu+k,i*Nu+l,gsl_matrix_get(R,k,l));
}

/**Compute H matrix
H=CM^T*Q*CM+R */
mpcptr->H=gsl_matrix_alloc(cHorizon*Nu,cHorizon*Nu);
mpcptr->H=MatAdd2(MatMul2(MatTrans(mpcptr->CM),MatMul2(Qbar,mpcptr->CM)),Rbar);

/**Construction of M matrix */
mpcptr->M=gsl_matrix_alloc(Ns*cHorizon,Ns);
gsl_matrix *tempM=gsl_matrix_alloc(Ns,Ns);
for(i=0;i<cHorizon;i++)
{
    tempM=MatMulrec(mpcptr->A,i+1);
    for(k=0;k<Ns;k++)
            for(l=0;l<Ns;l++)
                gsl_matrix_set(mpcptr->M,i*Ns+k,l,gsl_matrix_get(tempM,k,l));

}
//print2scr(M);
/**Construction of F Matrix
F=C^T*Q*M */
mpcptr->F=gsl_matrix_alloc(cHorizon*Nu,Ns);
mpcptr->F=MatMul2(MatTrans(mpcptr->CM),MatMul2(Qbar,mpcptr->M));
print2scr(mpcptr->F);

/**Construction of G matrix
G=M^T*Qbar*M+Q*/
mpcptr->G=gsl_matrix_alloc(Ns,Ns);
mpcptr->G=MatAdd2(MatMul2(MatTrans(mpcptr->M),MatMul2(Qbar,mpcptr->M)),Q);
//print2scr(mpcptr->G);

/** Setting up the QP problem to use with qpOASES
qpOASES require the number of variables, nV and the number of constraints
the number of variables is the number of inputs Nu
the number of constraints: to be completed*/
mpcptr->hval=malloc(mpcptr->H->size1*mpcptr->H->size2*sizeof(double));
mpcptr->fval=malloc(mpcptr->F->size1*mpcptr->F->size2*sizeof(double));
for(i=0;i<mpcptr->H->size1;i++)
    for(j=0;j<mpcptr->H->size2;j++)
       mpcptr->hval[j+i*mpcptr->H->size2]=gsl_matrix_get(mpcptr->H,i,j);

for(i=0;i<mpcptr->F->size1;i++)
    for(j=0;j<mpcptr->F->size2;j++)
       mpcptr->fval[(i*mpcptr->F->size2)+j]=gsl_matrix_get(mpcptr->F,i,j);




/*************************New Constraints Matrices*////
/// lb<u<ub//

mpcptr->lb=malloc(cHorizon*Nu*sizeof(double));
mpcptr->ub=malloc(cHorizon*Nu*sizeof(double));

for(i=0;i<cHorizon;i++)
    for(k=0;k<Nu;k++)
{
    mpcptr->lb[k+i*Nu]=lb[k];
    mpcptr->ub[k+i*Nu]=ub[k];
}
for(i=0;i<cHorizon*Nu;i++)
    printf("MPC Constraints %f %f\n",mpcptr->lb[i],mpcptr->ub[i]);



    /**STATE CONSTRAINTS*/

     /** lbA and ubA are the lower and upper constraints respectively
    remember
    x[k+i|k]=C_iu[k]+M_i.x[k]
    hence if lbA<=x[k+i]<=ubA
    we have
    lbA<=C_i.u[k]+M_i.x[k]
    lbA-M_i.x[k]<=C_i[uk
    Similarly
    C_i.u[k]<=ubA-M_i.x[k]

    C-->we already computed it above
    ubA is the constraints of the states but we need ubA for all horizon*/

mpcptr->CMval=malloc(cHorizon*cHorizon*Ns*Nu*sizeof(double));
    for(i=0;i<cHorizon*Ns;i++)
        for(j=0;j<cHorizon*Nu;j++)
            mpcptr->CMval[i*(cHorizon*Nu)+j]=gsl_matrix_get(mpcptr->CM,i,j);

print2scr(mpcptr->CM);

/** lbA and ubA are the lower and upper constraints respectively
    remember
    x[k+i|k]=C_iu[k]+M_i.x[k]
    hence if lbA<=x[k+i]<=ubA
    we have
    lbA<=C_i.u[k]+M_i.x[k]
    lbA-M_i.x[k]<=C_i[uk
    Similarly
    C_i.u[k]<=ubA-M_i.x[k]
    M.x[k] should be calculated at every step*/

    mpcptr->ubA=malloc(Ns*cHorizon*sizeof(double));
    mpcptr->lbA=malloc(Ns*cHorizon*sizeof(double));

    for(i=0;i<cHorizon;i++)
        for(j=0;j<Ns;j++)
    {
        mpcptr->ubA[j+i*Ns]=ubA[j];
        mpcptr->lbA[j+i*Ns]=lbA[j];
    }

/** xss and uss assigned 0*/

mpcptr->xss=malloc(Ns*sizeof(double));
mpcptr->uss=malloc(Nu*sizeof(*(mpcptr->uss)));

gsl_matrix_free(tempCM);
gsl_matrix_free(Rbar);
gsl_matrix_free(tempM);

}

/**calculation of steady state values at every time step */
double* MPCcalcSS(structMPC *mpcptr, double *refr,double *input_dist, double *output_dist,gsl_matrix *Bd,gsl_matrix *Cref)
{
    /***Steady state value calculation
    x_ss=Ax_ss+Bu_ss+Bd.x_dss
    y_ss=r=Cref.x_ss
    [x_ss;u_ss]=inv([(A-I) B;Cref 0])*[-Bd*x_d;r-w0]*/

    int Nr=Cref->size1;
    int Ns=mpcptr->A->size1;
    int Nu=mpcptr->B->size2;
    int sizedistinput=Bd->size2;
    int sizeref=Cref->size1;

    double sum=0;
    int i,j,k,l;


    gsl_matrix *BdX=gsl_matrix_alloc(Bd->size1,1); /**-Bd*X*/
    gsl_matrix_set_all(BdX,0);
    gsl_matrix *Refdist=gsl_matrix_alloc(Cref->size1,1);/** to store r-w0*/
    gsl_matrix_set_all(Refdist,0);

    gsl_matrix *RefBdX=gsl_matrix_alloc((Bd->size1)+(Cref->size1),1);
    gsl_matrix_set_all(RefBdX,0);

    gsl_matrix *SSmat=gsl_matrix_alloc(Ns+Nu,1);
    gsl_matrix_set_all(SSmat,0);



    ///Bd*input_dist
    for(i=0;i<BdX->size1;i++)
    {
        sum=0;
        for(k=0;k<Bd->size2;k++)
            sum=sum+gsl_matrix_get(Bd,i,k)*input_dist[k]*-1; /**minus sign included here */
        gsl_matrix_set(BdX,i,0,sum);
    }


    ///RFBDx=[Bdx*inputdist;Ref-outdis]
    for(i=0;i<sizeref;i++)
        gsl_matrix_set(Refdist,i,0,refr[i]-output_dist[i]);

    for(i=0;i<RefBdX->size1;i++)
    {
        if(i<BdX->size1)
        {
            gsl_matrix_set(RefBdX,i,0,gsl_matrix_get(BdX,i,0));



        }

        else
        {
              gsl_matrix_set(RefBdX,i,0,gsl_matrix_get(Refdist,i-BdX->size1,0));
        }

    }



//    gsl_matrix *tempA=gsl_matrix_alloc(Ns,Ns);
//    gsl_matrix *IdA=gsl_matrix_alloc(Ns,Ns);
//    gsl_matrix *invMat=gsl_matrix_alloc(Ns+Nr,Ns+Nu);
//
//    gsl_matrix_set_all(tempA,0);
//    gsl_matrix_set_identity(IdA);
//    gsl_matrix_set_all(invMat,0);
//
//    tempA=MatSub2(mpcptr->A,IdA);
//
//    for(i=0;i<Ns+Nr;i++)
//        for(j=0;j<Ns+Nu;j++)
//    {
//        if((i<(Ns)) && (j<(Ns)))
//            gsl_matrix_set(invMat,i,j,gsl_matrix_get(tempA,i,j));
//        else if((i<(Ns)) && (j>=(Ns)))
//            gsl_matrix_set(invMat,i,j,gsl_matrix_get(mpcptr->B,i,j-Ns));
//        else if((i>=(Ns)) && (j<(Ns)))
//            gsl_matrix_set(invMat,i,j,gsl_matrix_get(Cref,i-Ns,j));
//        else
//            gsl_matrix_set(invMat,i,j,0);
//
//
//    }
//
//     mpcptr->type=NORMAL;
     InitSteadyState(mpcptr,Cref);

     //SSmat=MatMul2(MatInv2(invMat),RefBdX);
     SSmat=MatMul2(MatInv2(mpcptr->SteadyState),RefBdX);
     double *result=malloc((Ns+Nu)*sizeof(double));


        for(i=0;i<Nu+Ns;i++)
        {
            result[i]=gsl_matrix_get(SSmat,i,0);
            printf("Result[%d]=%f\n",i,result[i]);
        }

 /** Since we have uss and xss as operating point, constraints should be adjusted at these operating points
 -lb<u+uss<ub*/
for(i=0;i<Ns;i++)
    mpcptr->xss[i]=result[i];
for(i=0;i<Nu;i++)
    mpcptr->uss[i]=result[Ns+i];









    //gsl_matrix_free(tempA);
    //gsl_matrix_free(IdA);
    gsl_matrix_free(BdX);
    gsl_matrix_free(Refdist);
    gsl_matrix_free(RefBdX);
    gsl_matrix_free(SSmat);
    //gsl_matrix_free(invMat);
     return result;

}

/**
    @brief Calculates the next control move. It solves a quadratic optimization problem.
    @param[in] *mpcptr Pointer to the structure holding the model predictive model and associated H,F,G matrix.
    @param[in] *xdata pointer to matrix holding the current value (x[i]-xss).
    @param[out] *u pointer to matrix holding the next N (prediction horizon) control.

*/
void MPC_Step(structMPC *mpcptr,gsl_matrix *xdata,gsl_matrix *u)
{


int nV,nC;
int i,j;
int Nu=mpcptr->B->size2;
int Ns=mpcptr->B->size1;
int cHorizon=mpcptr->contHor;

nV=cHorizon*Nu;
nC=cHorizon*Ns;


/**Creating the F*x[k] matrix*/

double *Fx=malloc(cHorizon*Nu*sizeof(double));
gsl_matrix *Fxmat=gsl_matrix_alloc(cHorizon*Nu,1);
Fxmat=MatMul2(mpcptr->F,xdata);
for(i=0;i<cHorizon*Nu;i++)
    Fx[i]=gsl_matrix_get(Fxmat,i,0);
printf("Before Mx");
/**Creating the M*x[k] matrix*/
    gsl_matrix *MX=gsl_matrix_alloc(cHorizon*Ns,1);
    double *MXval=malloc(cHorizon*Ns*sizeof(double));
    mpcptr->ubAMX=malloc(cHorizon*Ns*sizeof(double));
    mpcptr->lbAMX=malloc(cHorizon*Ns*sizeof(double));

    MX=MatMul2(mpcptr->M,xdata);

    for(i=0;i<cHorizon*Ns;i++)
        MXval[i]=gsl_matrix_get(MX,i,0);

    /**calculation of ubAMX and lbAMX*/

    for(i=0;i<cHorizon*Ns;i++)
    {
        mpcptr->ubAMX[i]=mpcptr->ubA[i]-MXval[i];
        mpcptr->lbAMX[i]=mpcptr->lbA[i]-MXval[i];
    }



/***adjustment of constraints so that lb<u+uss<ub*/
double *lbuss=malloc(cHorizon*Nu*sizeof(double));
double *ubuss=malloc(cHorizon*Nu*sizeof(double));
for(i=0;i<cHorizon;i++)
    for(j=0;j<Nu;j++)
{
   lbuss[j+i*Nu]=mpcptr->lb[j+i*Nu]-mpcptr->uss[j];
   ubuss[j+i*Nu]=mpcptr->ub[j+i*Nu]-mpcptr->uss[j];
   //mpcptr->lb[j+i*Nu]=mpcptr->lb[j+i*Nu]-mpcptr->uss[j];
   //mpcptr->ub[j+i*Nu]=mpcptr->ub[j+i*Nu]-mpcptr->uss[j];
}

/***adjustment of constraints so that lb<x+xss<ub*/
     for(i=0;i<cHorizon;i++)
        for(j=0;j<Ns;j++)
    {
       mpcptr->ubAMX[j+i*Ns]=mpcptr->ubAMX[j+i*Ns]-mpcptr->xss[j];
        mpcptr->lbAMX[j+i*Ns]=mpcptr->lbAMX[j+i*Ns]-mpcptr->xss[j];
    }


double *xopt=malloc(cHorizon*Nu*sizeof(double));
double *sol=malloc(cHorizon*Nu*sizeof(double));
int nWSR;
int q;
double *cputime;

cputime=NULL;
nWSR=40;



/*myTrajectoryRelatedOptimization (nV,nC,
                                   H, g, A,
                                   lb, ub,
                                   lbA, ubA,
                                   xopt,sol,
                                   &nWSR,cputime
                                  );*/

myTrajectoryRelatedOptimization (nV,nC,
                                   mpcptr->hval, Fx,mpcptr->CMval,
                                   lbuss, ubuss,
                                   mpcptr->lbAMX, mpcptr->ubAMX,
                                   xopt,sol,
                                   &nWSR,cputime
                                  );

printf("Initialised Correctly\n");
for(i=0;i<cHorizon*Nu;i++)
{
    printf("s[i]=%f\n",sol[i]);
    printf("xopt[i]=%f\n",xopt[i]);
}

for(i=0;i<Nu;i++)
    gsl_matrix_set(u,i,0,sol[i]);

    free(Fx);
    free(xopt);
    free(sol);
    free(lbuss);
    free(ubuss);
    free(MXval);
    free(cputime);
    gsl_matrix_free(Fxmat);
}


