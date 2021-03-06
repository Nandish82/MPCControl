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
/**
This function initialises the prediction type as well as the formulation type
and creates the corresponding matrices.
*/

void InitMPCType(structMPC *mpcptr,Model *modelptr,MPCType type,MPCPredictionType predtype)
{
    int Ns,Nu,Ny; ///number of states, inputs and outputs
    int i,j;

    /**Allocating the sizes */
    Ns=modelptr->A->size1;
    Nu=modelptr->B->size2;
    Ny=modelptr->C->size1;

    mpcptr->type==type;
    mpcptr->predtype=predtype;

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
        mpcptr->C=gsl_matrix_alloc(Ny+Nu,Ns+Nu);
        mpcptr->D=gsl_matrix_alloc(Ny+Nu,Nu);

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
         for(i=0;i<(Ny+Nu);i++)
            for(j=0;j<(Ns+Nu);j++)
        {
            if((i<Ny)&&(j<Ns))
                gsl_matrix_set(mpcptr->C,i,j,gsl_matrix_get(modelptr->C,i,j));
            else if((i<Ny)&&(j>=Ns))
                        gsl_matrix_set(mpcptr->C,i,j,gsl_matrix_get(modelptr->D,i,j-Ns));
            else if((i>=Ny)&&(j>=Ns))
                if((i-Ny)==(j-Ns))
                        gsl_matrix_set(mpcptr->C,i,j,1); ///added input as output also to facilitate the use of constraints.

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

int InitSteadyState(structMPC *mpcptr,double *Cref,int Ntr)
{
    int Ns,Nu,Ny,i,j;

    ///Ntr is number of tracked outputs.

    if(!((mpcptr->type==DELTA)||(mpcptr->type==NORMAL)))
    {
        printf("Error: MPC type has not been initialised");
        return 2;
    }
    if(mpcptr->type==DELTA)
    {
        Nu=mpcptr->B->size2;
        Ns=mpcptr->A->size1-Nu;
        Ny=mpcptr->C->size1-Nu;
        printf("Delta Ns:%d Nu:%d Ny:%d\n",Ns,Nu,Ny);
    }
    else
    {
        Nu=mpcptr->B->size2;
        Ns=mpcptr->A->size1;
        Ny=mpcptr->C->size1;
        printf("Ns:%d Nu:%d Ny:%d\n",Ns,Nu,Ny);
    }
    if(Ntr!=Nu)
    {
        printf("Number of tracked outputs is greater than number of inputs. Infeasible");
        return 2;
    }

        mpcptr->SteadyState=gsl_matrix_alloc(Ns+Ntr,Ns+Nu);
        for(i=0;i<Ns+Ntr;i++)
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
                    gsl_matrix_set(mpcptr->SteadyState,i,j,Cref[j+(i-Ns)*Ntr]);
                 else if((i>=Ns)&&(j>=Ns))
                    gsl_matrix_set(mpcptr->SteadyState,i,j,0);

    mpcptr->uss=malloc((Nu)*sizeof(double));
    mpcptr->xss=malloc((Ns)*sizeof(double));
    mpcptr->yss=malloc((Ny)*sizeof(double));

  return 1;
}

void InitMPC(structMPC *mpcptr,Model *md,int cHorizon,gsl_matrix *Q,gsl_matrix *P,gsl_matrix *R,double *lb,double *ub,double *lbA,double *ubA)
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

mpcptr->Ts=md->Ts;
mpcptr->predHor=cHorizon;
int Ns,Nu,Ny,i,j,k,l;

Ns=md->A->size1;
Nu=md->B->size2;
Ny=md->C->size1;

mpcptr->A=gsl_matrix_alloc(Ns,Ns);
mpcptr->B=gsl_matrix_alloc(Ns,Nu);
mpcptr->C=gsl_matrix_alloc(Ny,Ns);
mpcptr->D=gsl_matrix_alloc(Ny,Nu);

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
//     InitSteadyState(mpcptr,Cref);

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
int cHorizon=mpcptr->predHor;

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
int MPCpredmat(structMPC *mpcptr,int Np,int Nc)
{

int Ns,Nu,Ny,i,j,k,l,Nsy;

gsl_matrix *C;

if(!((mpcptr->predtype==STATE) || (mpcptr->predtype==OUTPUT)))
{
    printf("MPC has not yet been initialised");
    return 2;
}

mpcptr->predHor=Np;
mpcptr->contHor=Nc; ///

Ns=mpcptr->A->size1;
Nu=mpcptr->B->size2;
Ny=mpcptr->C->size1;


///Nsy determines whether state or output matrices should be formed
if(mpcptr->predtype==STATE)
{
    ///We set the C matrix to identity so that we can use same equation for both state and output predictions.
    C=gsl_matrix_alloc(Ns,Ns);
    gsl_matrix_set_identity(C);
    Nsy=C->size1;
        if((mpcptr->Q->size1!=Ns) && (mpcptr->P->size1!=Ns))
        {
            printf("[Error] Check Size of matrix Q and/or P\n");
            return 2;
        }
}
else
{
    C=gsl_matrix_alloc(Ny,Ns);
    gsl_matrix_memcpy(C,mpcptr->C);
    Nsy=C->size1;

    if((mpcptr->Q->size1!=Ny) && (mpcptr->P->size1!=Ns))
        {
            printf("[Error] Check Size of matrix Q and/or P\n");
            return 2;
        }
}


/**

Np is the prediction Horizon
Nc is the control Horizon

    //matrix sizes
    // Sux=(Np x Ns)*(Nc x Nu)
    // Suy=(Np x Ny) *(Nc*Nu)
    // Sxx=(Np*Ns) * (Ns)
    // Sxy=(Np*Ny) *(Ns) //care should be taken here as we made Ns=Ny...the column of this matrix is the number of states
    // Qtildex=(Np*Ns) * (Np*Ns)
    // Qtildey=(Np*Ny) * (Np*Ny)
    // Rtildex=(Nc*Nu) * (Nc*Nu)
    // Rtildex=Rtildey
    // H=Su'*Qtildexy*Su+Rt1ilde (Nc x Nu)*(Nc x Nu)
    // F=Su'*Qtildexy*Sx
    // G=


*/

///Declare temporary matrices, they all have to be freed afterwards
mpcptr->Su=gsl_matrix_alloc(Np*Nsy,Nc*Nu);
gsl_matrix_set_all(mpcptr->Su,0);
printf("\n Here \n");
mpcptr->Sx=gsl_matrix_alloc(Np*Nsy,Ns);
gsl_matrix_set_all(mpcptr->Sx,0);
gsl_matrix *Qtilde=gsl_matrix_alloc(Np*Nsy,Np*Nsy);
gsl_matrix_set_zero(Qtilde);
gsl_matrix *Rtilde=gsl_matrix_alloc(Nc*Nu,Nc*Nu);
gsl_matrix_set_zero(Rtilde);
///intermediate matrices to be freed also
gsl_matrix *tempCAB=gsl_matrix_alloc(Nsy,Nu);
gsl_matrix *tempCA=gsl_matrix_alloc(Nsy,Ns);
gsl_matrix_set_zero(tempCAB);
print2scr(tempCA);
///Populating Su[Np x Nsy,Nc x Nu] matrix
for(i=0;i<Np;i++)
{
    for(k=0;k<Nsy;k++)
        for(l=0;l<Ns;l++)

        {
            tempCA=MatMul2(C,MatMulrec(mpcptr->A,i+1));
             gsl_matrix_set(mpcptr->Sx,i*Nsy+k,l,gsl_matrix_get(tempCA,k,l));///C*(A)^(n-1)
        }
        gsl_matrix_set_zero(tempCAB);
    for(j=0;j<Nc;j++)
    {
        if((i>=j)&& (i<Nc))
           {
               tempCAB=MatMul2(C,MatMul2(MatMulrec(mpcptr->A,i-j),mpcptr->B));
               for(k=0;k<Nsy;k++)
                for(l=0;l<Nu;l++)
                gsl_matrix_set(mpcptr->Su,i*Nsy+k,j*Nu+l,gsl_matrix_get(tempCAB,k,l));
           }
        gsl_matrix_set_zero(tempCAB);
        if((i>=j)&&i>=Nc)
        {
            if(j==Nc-1)
            {
                for(k=0;k<=(i-Nc+1);k++)
                {
                    tempCAB=MatAdd2(tempCAB,MatMul2(C,MatMul2(MatMulrec(mpcptr->A,k),mpcptr->B)));
                }
                 for(k=0;k<Nsy;k++)
                    for(l=0;l<Nu;l++)
                        gsl_matrix_set(mpcptr->Su,i*Nsy+k,j*Nu+l,gsl_matrix_get(tempCAB,k,l));
                    gsl_matrix_set_zero(tempCAB);

            }
            else
            {
                tempCAB=MatMul2(C,MatMul2(MatMulrec(mpcptr->A,i-j),mpcptr->B));
                for(k=0;k<C->size1;k++)
                    for(l=0;l<Nu;l++)
                        gsl_matrix_set(mpcptr->Su,i*Nsy+k,j*Nu+l,gsl_matrix_get(tempCAB,k,l));

            }
        }
    }




} ///end for loop

///Qtilde
for(i=0;i<Np;i++)
{
    if(i<Np-1)
    {
        for(k=0;k<Nsy;k++)
            for(l=0;l<Nsy;l++)
            gsl_matrix_set(Qtilde,i*Nsy+k,i*Nsy+l,gsl_matrix_get(mpcptr->Q,k,l));
    }
    else
    {
        for(k=0;k<Nsy;k++)
            for(l=0;l<Nsy;l++)
            gsl_matrix_set(Qtilde,i*Nsy+k,i*Nsy+l,gsl_matrix_get(mpcptr->P,k,l));
    }
}///end for loop

///Rtilde
for(i=0;i<Nc;i++)
{

        for(k=0;k<Nu;k++)
            for(l=0;l<Nu;l++)
            gsl_matrix_set(Rtilde,i*Nu+k,i*Nu+l,gsl_matrix_get(mpcptr->R,k,l));
}



mpcptr->H=gsl_matrix_alloc(Nc*Nu,Nc*Nu);
mpcptr->H=MatAdd2(MatMul2(MatTrans(mpcptr->Su),MatMul2(Qtilde,mpcptr->Su)),Rtilde);
print2scr(mpcptr->H);

mpcptr->F=gsl_matrix_alloc(Nc*Nu,Ns);
mpcptr->F=MatMul2(MatTrans(mpcptr->Su),MatMul2(Qtilde,mpcptr->Sx));
//print2scr(mpcptr->F);

mpcptr->G=gsl_matrix_alloc(Ns,Ns);
mpcptr->G=MatAdd2(MatMul2(MatTrans(C),MatMul2(mpcptr->Q,C)),MatMul2(MatTrans(mpcptr->Sx),MatMul2(Qtilde,mpcptr->Sx)));
/***Casting into qpoases form*/
mpcptr->hval=malloc(Nc*Nu*Nc*Nu*sizeof(double));
for(i=0;i<Nc*Nu;i++)
    for(j=0;j<Nc*Nu;j++)
    mpcptr->hval[j+i*(Nc*Nu)]=gsl_matrix_get(mpcptr->H,i,j);

mpcptr->fval=malloc(Nc*Nu*Ns*sizeof(double));
for(i=0;i<Nc*Nu;i++)
    for(j=0;j<Ns;j++)
    mpcptr->fval[j+i*(Ns)]=gsl_matrix_get(mpcptr->F,i,j);

mpcptr->gval=malloc(Ns*Ns*sizeof(double));
for(i=0;i<Ns;i++)
    for(j=0;j<Ns;j++)
    mpcptr->gval[j+i*(Ns)]=gsl_matrix_get(mpcptr->G,i,j);
//print2scr(mpcptr->G);
/***free matrices**/
gsl_matrix_free(Qtilde);
gsl_matrix_free(Rtilde);
gsl_matrix_free(tempCAB);
gsl_matrix_free(tempCA);

return 0;

}

/** \brief This function creates the matrices Axcon,bxcon,bucon.
    In qpoases, the constraints has to be in the form
    lbx<Ax<ubx
    lbu<x<ubu
    Since we have to optimise for the variable u or delta u in MPC control
    our constraints become as follows
    For input constraints
    lbu<u+uss<ubu, size of u= Nc*Nu
    For state/output contraints
    lbx<(x+xss)<ubx
    lbx-xss<x<ubx-xss
    lbx-xss<Su.u+Sx.x<ubx-uss
    lbs-xss-Sx.x<Su.u<ubx-xss-Sx.x..............1

    As we see all the elements are already created in MPCpredmat and MPCsteadystate.
    qpoases requires
    number of variables to be optimised which is [Nc*Nu]
    number of constraints is [Np*Ns]----rows of matrice Su

    lbx can be either state or output constraint. This function casts everything into double row form

    There will be one function InitMPCconstraints where the constraints will be initialised and one StepMPCconstraints
    where the constraints will be calculated and feeded to the qpoases algorithm.

    Su is also called CM variable sometimes and Sx is called Mm variable sometimes.

 *
 * \param
 * \param
 * \return
 *
 */

int InitMPCconstraints(structMPC *mpcptr,double *lbu,double *ubu, double *lbxy, double *ubxy,double *lbdelta,double *ubdelta)
{

    ///size of lbu = Nu x 1
    ///size of lbx =(Ns or Ny) x 1

    int i,j,k,l; ///counters
    int Ns,Nu,Ny,Np,Nc,Nsy;

    double *lb1;
    double *ub1;

    double *lbA1;
    double *ubA1;

    Ns=mpcptr->A->size1;
    Nu=mpcptr->B->size2;
    Ny=mpcptr->C->size1;

    Nc=mpcptr->contHor;
    Np=mpcptr->predHor;


/***Creates constraints according to type and predtype***/
/// for Delta lbA=[constraints states/output constraints input]
///           lb=[constraints rate(lbdelta)]



    if(mpcptr->predtype==OUTPUT)
    {
        Nsy=Ny;
    }
    else if(mpcptr->predtype==STATE)
    {
        Nsy=Ns;

    }
    else
    {
        printf("\n Prediction type has not been initialised");
        return 2; ///unsucessful
    }

        lb1=malloc(Nu*sizeof(double));
        ub1=malloc(Nu*sizeof(double));
        lbA1=malloc(Nsy*sizeof(double));
        ubA1=malloc(Nsy*sizeof(double));



    if(mpcptr->type==DELTA)
    {


        for(i=0;i<Nu;i++)
        {
            lb1[i]=lbdelta[i];
            ub1[i]=ubdelta[i];
        }

        for(i=0;i<Nsy;i++)
        {
            if(i<(Nsy-Nu))
            {
                lbA1[i]=lbxy[i];
                ubA1[i]=ubxy[i];
            }
            if(i>=(Nsy-Nu))
            {
                lbA1[i]=lbu[i-Nsy+Nu];
                ubA1[i]=ubu[i-Nsy+Nu];
            }
        }
    }
    else
    {
         for(i=0;i<Nu;i++)
        {
            lb1[i]=lbu[i];
            ub1[i]=ubu[i];
        }

        for(i=0;i<Nsy;i++)
        {
            lbA1[i]=lbxy[i];
            ubA1[i]=ubxy[i];
        }



    }
     mpcptr->lb=malloc(Nc*Nu*sizeof(double));
     mpcptr->ub=malloc(Nc*Nu*sizeof(double));

     mpcptr->lbuss=malloc(Nc*Nu*sizeof(double));
     mpcptr->ubuss=malloc(Nc*Nu*sizeof(double));

     mpcptr->lbA=malloc(Np*Nsy*sizeof(double));
     mpcptr->ubA=malloc(Np*Nsy*sizeof(double));

     mpcptr->lbAxss=malloc(Np*Nsy*sizeof(double));
     mpcptr->ubAxss=malloc(Np*Nsy*sizeof(double));

     mpcptr->suval=malloc(Np*Nsy*Nc*Nu*sizeof(double));

    ///populating the matrices due to control and prediction horizon..for each of them the constraints has to be repeated
    for(i=0;i<Nc;i++)
        for(j=0;j<Nu;j++)
        {
            mpcptr->lb[i*Nu+j]=lb1[j];
            mpcptr->ub[i*Nu+j]=ub1[j];

            mpcptr->lbuss[i*Nu+j]=lb1[j];
            mpcptr->ubuss[i*Nu+j]=ub1[j];
        }
    for(i=0;i<Np;i++)
        for(j=0;j<Nsy;j++)
    {
        mpcptr->lbA[i*Nsy+j]=lbA1[j];
        mpcptr->ubA[i*Nsy+j]=ubA1[j];

        mpcptr->lbAxss[i*Nsy+j]=lbA1[j];
        mpcptr->ubAxss[i*Nsy+j]=ubA1[j];
    }

    ///suval is required for the qpoases problem. So here we convert the Su matrix to suval double.
    for(i=0;i<Np*Nsy;i++)
        for(j=0;j<Nc*Nu;j++)
        mpcptr->suval[(i*Nc*Nu+j)]=gsl_matrix_get(mpcptr->Su,i,j);

    ///Assigning nCon and nVar

    mpcptr->nCons=Nsy*Np; ///number of constraints
    mpcptr->nVar=Nc*Nu; ///number of variables to optimize

    free(lb1);
    free(ub1);
    free(lbA1);
    free(ubA1);
    return 0;

}
int StepMPCconstraints(structMPC *mpcptr,double *xdata)
{
    /**
    lbs-xss-Sx.x<Su.u<ubx-xss-Sx.x.
    i.e lbA=lbx-xss-Sx.x
        ubA=ubx-xss-Sx.x
        lb=lbu-uss
        ub=ubu-uss

        We assume that the function StepSteadyMPC has calculated the values of xss and uss and stored it in the
        variables mpcptr->xss and mpcptr->uss and has already initialised these latter arrays of doubles to their
        correct size.


    */


    int Ns,Nu,Ny,Np,Nc,Nsy,i,j,k;
    double *ssxy,*sstemp;

    Ns=mpcptr->A->size1;
    Nu=mpcptr->B->size2;
    Ny=mpcptr->C->size1;

    Np=mpcptr->predHor;
    Nc=mpcptr->contHor;

     if(mpcptr->predtype==OUTPUT)
    {
        Nsy=Ny;
        sstemp=&mpcptr->yss[0];
    }
    else if(mpcptr->predtype==STATE)
    {
        Nsy=Ns;
        sstemp=&mpcptr->xss[0];

    }
    else
    {
        printf("\n Prediction type has not been initialised");
        return 2; ///unsucessful
    }

    if(mpcptr->type==DELTA)
    {
        ssxy=malloc(Nsy*sizeof(double));
        for(i=0;i<Nsy;i++)
        {
             if(i<Nsy-Nu)
             {
                 ssxy[i]=sstemp[i];
                 printf("ssxy[%d]:%f\n",i,ssxy[i]);
             }
            else
            {
                ssxy[i]=mpcptr->uss[i-(Nsy-Nu)];
                 printf("ssxy[%d]:%f\n",i-Nsy+Nu,ssxy[i]);
            }


        }

    }
    else
    {
         ssxy=malloc(Nsy*sizeof(double));
          for(i=0;i<Nsy;i++)
                ssxy[i]=sstemp[i];
    }

    double temp=0.0;

    for(i=0;i<Np;i++)
    {
           for(j=0;j<Ns;j++)
            {
                temp=temp+xdata[j]*gsl_matrix_get(mpcptr->Sx,i,j);///Sx.
            }

            for(j=0;j<Nsy;j++)
            {
                 mpcptr->lbAxss[j+i*Nsy]=mpcptr->lbA[j+i*Nsy]-temp-ssxy[j];///lbx=lbx-Sx.x-xss
                 mpcptr->ubAxss[j+i*Nsy]=mpcptr->ubA[j+i*Nsy]-temp-ssxy[j];
            }
            printf("temp:%f\n",temp);
            temp=0;


    }
    if(mpcptr->type==NORMAL)
    {
        for(i=0;i<Nc*Nu;i++)
    {
       mpcptr->lbuss[i]=mpcptr->lb[i]-mpcptr->uss[i%Nu];///lb=lb-uss REMEMBER: lb is a Nc x Nu matrix whereas uss is an Nu, so values have to be repeated
                                                        ///every Nu times. Same for xss
        mpcptr->ubuss[i]=mpcptr->ub[i]-mpcptr->uss[i%Nu];
    }
    }
    else
    {
        for(i=0;i<Nc*Nu;i++)
    {
       mpcptr->lbuss[i]=mpcptr->lb[i];///lb=lb-uss REMEMBER: lb is a Nc x Nu matrix whereas uss is an Nu, so values have to be repeated
                                                        ///every Nu times. Same for xss
        mpcptr->ubuss[i]=mpcptr->ub[i];
    }
    }


    return 0;

}

// TODO (ncalcha#1#): Find a way to assign weight depending on whether we are using state predicitions or output predicitions. Furthermore, it also depends whetehter we are using Delta representation or normal representation. ...
//

/****
This function assigns weights
**/

int AssignMPCWeights(structMPC *mpcptr,double *qxy,double *rinput,double *rrate)
{
    int Ns,Ny,Nu,i,j,Nsy;


    Ns=mpcptr->A->size1;
    Nu=mpcptr->B->size2;
    Ny=mpcptr->C->size1;

    if (mpcptr->predtype==OUTPUT)
        Nsy=Ny;
    else
        Nsy=Ns;



    /**allocate matrices*/
    mpcptr->Q=gsl_matrix_alloc(Nsy,Nsy);
    gsl_matrix_set_zero(mpcptr->Q);
    mpcptr->R=gsl_matrix_alloc(Nu,Nu);
    gsl_matrix_set_zero(mpcptr->R);
    mpcptr->P=gsl_matrix_alloc(Nsy,Nsy);
    gsl_matrix_set_zero(mpcptr->P);

    for(i=0;i<Nsy;i++)
    {
        if(mpcptr->type==DELTA)
        {
            if(i<Nu)
                gsl_matrix_set(mpcptr->R,i,i,rrate[i]);


            if(i<(Nsy-Nu))
                gsl_matrix_set(mpcptr->Q,i,i,qxy[i]);
            else
                gsl_matrix_set(mpcptr->Q,i,i,rinput[i-Nsy]);
        }
        else if(mpcptr->type==NORMAL)
        {
            if(i<Nu)
            gsl_matrix_set(mpcptr->R,i,i,rinput[i]);

            gsl_matrix_set(mpcptr->Q,i,i,qxy[i]);

        }

    }

    gsl_matrix_memcpy(mpcptr->P,mpcptr->Q);






}


int StepSteadyState(structMPC *mpcptr,double *refr,double *inputdist,double *outputdist, double *Bd, int Nid)
{
   ///Nid is the number of input disturbances
    ///we dont need the number of output disturbance since we can obtain it from the number of tracked outputs.

    /**
    **  We are going to
    ** inv[SteadyStateMatrix] * [-Bd*inputdist]    =[xss]
                                [refr-outputdist]   [uss]
    ** Here we are using the normal matrix not the DELTA formulation
    */
   int Ns,Nu,Ny,Nss,i,j,k;
    double temp=0;

    Nss=mpcptr->SteadyState->size1; ///size of steady state matrix

    Nu=mpcptr->B->size2;
    ///we cannot have more tracked outputs than inputs
    ///else our steady state matrix is not square
    /// and cannot be inverted.
    if(mpcptr->type==DELTA)
    {
        Ny=mpcptr->C->size1-mpcptr->B->size2;
        Ns=mpcptr->A->size1-mpcptr->B->size2; ///original size of number of states
    }
    else
    {
         Ny=mpcptr->C->size1;
          Ns=mpcptr->A->size1;
    }


    gsl_matrix *RefBdx=gsl_matrix_alloc(Nss,1);
    gsl_matrix_set_zero(RefBdx);
    for(i=0;i<Ns+Nu;i++)
    {
        if(i<Ns)
    {

            for(k=0;k<Nid;k++)
            {
                temp=temp+Bd[i+k*Ns]*inputdist[k]; ///Bd*x

            }


         gsl_matrix_set(RefBdx,i,0,-temp);

         temp=0;
    }
        if(i>=Ns)
            gsl_matrix_set(RefBdx,i,0,refr[i-Ns]-outputdist[i-Ns]);
    }

    printf("RefBdx\n");
    print2scr(RefBdx);
    gsl_matrix *ssmat=gsl_matrix_alloc(Nss,1);
    ssmat=MatMul2(MatInv2(mpcptr->SteadyState),RefBdx);




    for(i=0;i<Ns;i++)
    {
         mpcptr->xss[i]=gsl_matrix_get(ssmat,i,0);
    }

    for(i=0;i<Nu;i++)
    {
        mpcptr->uss[i]=gsl_matrix_get(ssmat,i+Ns,0);
    }


    ///Populating mpcptr->yss mat yss=C*xss///
     temp=0;
    for(i=0;i<Ny;i++)
    {
        for(j=0;j<Ns;j++)
        {
             temp=temp+gsl_matrix_get(mpcptr->C,i,j)*mpcptr->xss[j];

        }
        mpcptr->yss[i]=temp;
        temp=0;


    }



    gsl_matrix_free(RefBdx);
    gsl_matrix_free(ssmat);
return 0;
}
int StepMPC(structMPC *mpcptr,double *x,double *u)
{
    int Ns,Nu,Ny,i,j,Nsy,Np,Nc;
    double *ssxy;
    double *fvalx;


///Printing of results
FILE *fp;
fp=fopen("step.txt","a+");

Nc=mpcptr->contHor;
Np=mpcptr->predHor;
    if(mpcptr->type==DELTA)
    {
        Nu=mpcptr->B->size2;
        Ns=mpcptr->A->size1;
        Ny=mpcptr->C->size1;

    }

    else
    {
        Nu=mpcptr->B->size2;
        Ns=mpcptr->A->size1;
        Ny=mpcptr->C->size1;
    }
        printf("Nu:%d Ns:%d Ny:%d\n",Nu,Ns,Ny);
    if(mpcptr->predtype==OUTPUT)
    {
        Nsy=Ny;
        //ssxy=&mpcptr->yss[0];
    }

        else
        {
            Nsy=Ns;
            //ssxy=&mpcptr->xss[0];
        }

fprintf(fp,"Values passed:\n");
for(i=0;i<Ns;i++)
    fprintf(fp,"x[%d]:%f\n",i,x[i]);


fvalx=malloc(mpcptr->contHor*Nu*sizeof(double));
double *xopt=malloc(mpcptr->contHor*Nu*sizeof(double));
double *sol=malloc(mpcptr->contHor*Nu*sizeof(double));
int nWSR;
int q;
double *cputime;

cputime=NULL;
nWSR=40;


for(i=0;i<mpcptr->contHor*Nu;i++)
{
     for(j=0;j<Ns;j++)
    {
        fvalx[i]=fvalx[i]+mpcptr->fval[j+i*Ns]*x[j];
    }
    fprintf(fp,"fvalx[%d]: %f\n",i,fvalx[i]);
}


/*myTrajectoryRelatedOptimization (nV,nC,
                                   H, g, A,
                                   lb, ub,
                                   lbA, ubA,
                                   xopt,sol,
                                   &nWSR,cputime*/

//for(i=0;i<Nc*Nu;i++)
//    printf("%d:%f<=uss<=%f\n",i,mpcptr->lbuss[i],mpcptr->ubuss[i]);
//
//
//for(i=0;i<Np*mpcptr->C->size1;i++)
//    printf("%d:%f<=Su.xss<=%f\n",i,mpcptr->lbA[i],mpcptr->ubA[i]);

//StepMPCconstraints(mpcptr,x);
///****************************************TO BE REMOVED LATER************************/


    double *MXval=malloc(Np*Nsy*sizeof(double));

    printf("Nsy:%d Ns:%d\n",Nsy,Ns);

    ///Calculate Sx*x needed for constraints
    for(i=0;i<Np*Nsy;i++)
    {
        MXval[i]=0;
        for(j=0;j<Ns;j++)
        {
            MXval[i]=MXval[i]+gsl_matrix_get(mpcptr->Sx,i,j)*x[j];

        }
        fprintf(fp,"MXval[%d]: %f\n",i,MXval[i]);
    }




    /**calculation of ubAMX and lbAMX*/

    for(i=0;i<Np*Nsy;i++)
    {
        mpcptr->ubAxss[i]=mpcptr->ubA[i]-MXval[i];
        mpcptr->lbAxss[i]=mpcptr->lbA[i]-MXval[i];
        fprintf(fp,"Mpc->lbAxss[%d]: %f\n",i,mpcptr->lbAxss[i]);
        fprintf(fp,"Mpc-<ubAxxl[%d]: %f\n",i,mpcptr->ubAxss[i]);

    }



/***adjustment of constraints so that lb<u+uss<ub*/
double *lbuss=malloc(Nc*Nu*sizeof(double));
double *ubuss=malloc(Nc*Nu*sizeof(double));
for(i=0;i<Nc;i++)
    for(j=0;j<Nu;j++)
{
   lbuss[j+i*Nu]=mpcptr->lb[j+i*Nu];//-mpcptr->uss[j];
   ubuss[j+i*Nu]=mpcptr->ub[j+i*Nu];//-mpcptr->uss[j];
   //mpcptr->lb[j+i*Nu]=mpcptr->lb[j+i*Nu]-mpcptr->uss[j];
   //mpcptr->ub[j+i*Nu]=mpcptr->ub[j+i*Nu]-mpcptr->uss[j];
}

/***adjustment of constraints so that lb<x+xss<ub*/
for(i=0;i<Np;i++)
        for(j=0;j<Nsy;j++)
    {
       mpcptr->ubAxss[j+i*Ns]=mpcptr->ubAxss[j+i*Ns];//-ssxy[j];
       mpcptr->lbAxss[j+i*Ns]=mpcptr->lbAxss[j+i*Ns];//-ssxy[j];
    }

///************************************************************************************
myTrajectoryRelatedOptimization (mpcptr->nVar,mpcptr->nCons,
                                   mpcptr->hval, fvalx,mpcptr->suval,
                                   mpcptr->lbuss, mpcptr->ubuss,
                                   mpcptr->lbAxss, mpcptr->ubAxss,
                                   xopt,sol,
                                   &nWSR,cputime
                                  );
for(i=0;i<Nu;i++)
{
    u[i]=sol[i];
    fprintf(fp,"u/delta[%d]:%f\n",i,sol[i]);
}
fclose(fp);
free(MXval);
free(fvalx);
free(lbuss);
free(ubuss);
free(xopt);
free(sol);

//    free(xmat);
}

void print2FileMPC(structMPC *mpcptr,char *filename)
{
    FILE *fp;
    fp=fopen(filename,"w");

    int Nsy,Nu,Ns;
    int i,j;
    int typeflag;
    int formflag;
    int predflag;
    char stringform[50]="Type is: ";

    Nu=mpcptr->B->size2;
    Ns=mpcptr->A->size1;
    if(mpcptr->predtype==OUTPUT)
    {
        predflag=1;
        strcat(stringform,"OUTPUT/");
        Nsy=mpcptr->C->size1;

    }

    else
    {
         predflag=0;
         strcat(stringform,"STATE/");
         Nsy=mpcptr->A->size1;

    }

    if(mpcptr->type==DELTA)
    {
         formflag=1;
         strcat(stringform,"DELTA");

    }
    else
    {
        formflag=0;
        strcat(stringform,"NORMAL");

    }

    /**Prints delta/normal and state/output formulation**/
    fprintf(fp,"%s\n",stringform);





    /***Prints A,B,C,D Matrices*/
    fprintf(fp,"Matrix A:\n");
    print2FileMat(mpcptr->A,fp);
    fprintf(fp,"Matrix B:\n");
    print2FileMat(mpcptr->B,fp);
    fprintf(fp,"Matrix C:\n");
    print2FileMat(mpcptr->C,fp);
    fprintf(fp,"Matrix D:\n");
    print2FileMat(mpcptr->D,fp);

    /**Prints weights**/
     fprintf(fp,"Matrix Q:\n");
    print2FileMat(mpcptr->Q,fp);
    fprintf(fp,"Matrix P:\n");
    print2FileMat(mpcptr->P,fp);
    fprintf(fp,"Matrix R:\n");
    print2FileMat(mpcptr->R,fp);

    /**Control and Prediciton Horizon*/
    fprintf(fp,"Prediction Horizon is: %d \nControl Horizon is: %d \n",mpcptr->predHor,mpcptr->contHor);

    /**Prediction Matrices**/
    fprintf(fp,"Matrix H:\n");
    print2FileMat(mpcptr->H,fp);
    fprintf(fp,"Matrix F:\n");
    print2FileMat(mpcptr->F,fp);
    fprintf(fp,"Matrix G:\n");
    print2FileMat(mpcptr->G,fp);


    /***Constraint Matrices*/
    fprintf(fp,"Constraints on the input/delta matrix\n");

    for(i=0;i<(mpcptr->contHor*(Nu));i++)
    fprintf(fp,"%f<=u<=%f\n",mpcptr->lb[i],mpcptr->ub[i]);


    for(i=0;i<(mpcptr->predHor*Nsy);i++)
    {
        if(i%Nsy==0)
            fprintf(fp, "Prediction Step %d\n",i/Nsy+1);
        fprintf(fp,"%f<=Su.x<=%f\n",mpcptr->lbA[i],mpcptr->ubA[i]);
    }

    fprintf(fp,"Steady State Matrix\n");
    print2FileMat(mpcptr->SteadyState,fp);
    fprintf(fp,"Su Matrix:\n");
    print2FileMat(mpcptr->Su,fp);
    fprintf(fp,"Sx Matrix:\n");
    print2FileMat(mpcptr->Sx,fp);


for(i=0;i<mpcptr->contHor*Nu;i++)
{

    for(j=0;j<Ns;j++)
    {
        fprintf(fp,"fval[%d]:%f:%f\n",j+i*Ns,mpcptr->fval[j+i*Ns],gsl_matrix_get(mpcptr->F,i,j));
    }
    fprintf(fp,"\n");
}

    fclose(fp);



}

void printSteadyState(structMPC *mpcptr)
{
    int Ns,Ny,Nu,i,j;
     Nu=mpcptr->B->size2;
     Ns=mpcptr->A->size1;
     Ny=mpcptr->C->size1;
    if(mpcptr->type==DELTA)
    {
       Ns=Ns-Nu;
       Ny=Ny-Nu;
    }

    for(i=0;i<Ns;i++)
        printf("xss[%d]=%f\n",i,mpcptr->xss[i]);
    for(i=0;i<Nu;i++)
        printf("uss[%d]=%f\n",i,mpcptr->uss[i]);
}
