#include "model.h"


void InitModel(Model *m,gsl_matrix *A,gsl_matrix *B,gsl_matrix *C,gsl_matrix *D,gsl_matrix *X0)
{
int Ns=A->size1;
int Nu=B->size2;
int Ny=C->size1;

m->A=gsl_matrix_alloc(Ns,Ns);
m->B=gsl_matrix_alloc(Ns,Nu);
m->C=gsl_matrix_alloc(Ny,Ns);
m->D=gsl_matrix_alloc(Ny,Nu);
m->X0=gsl_matrix_alloc(Ns,1);

printf("Adress of m->Bd before allocation:%p\n",m->Bd);
m->Bd=gsl_matrix_alloc(Ns,1);
printf("Adress of m->Bd after allocation:%p\n",m->Bd);
gsl_matrix_set_zero(m->Bd);
printf("In Init Model\n");
printf("m->Bd");
print2scr(m->Bd);


gsl_matrix_memcpy(m->A,A);
gsl_matrix_memcpy(m->B,B);
gsl_matrix_memcpy(m->C,C);
gsl_matrix_memcpy(m->D,D);
gsl_matrix_memcpy(m->X0,X0);

}
void InitSimModel(Model *m,gsl_matrix *X0,double Ts,double tsim,gsl_matrix *Bd)
{
    int i,j;
    long Ndatapoints=tsim/Ts+1; //number of datapoints to be simulated
    gsl_matrix *y0=gsl_matrix_alloc(m->C->size1,1);
    m->statedata=gsl_matrix_alloc(m->A->size1,Ndatapoints);
    m->outputdata=gsl_matrix_alloc(m->C->size1,Ndatapoints);
    m->currxdata=gsl_matrix_alloc(m->A->size1,1);
    m->currydata=gsl_matrix_alloc(m->C->size1,1);
    m->inputdata=gsl_matrix_alloc(m->B->size2,Ndatapoints);
    m->X0=X0;
    m->currxdata=X0;
    y0=MatMul2(m->C,X0);
    m->currydata=y0;
    gsl_matrix_set_all(m->statedata,0);
    gsl_matrix_set_all(m->outputdata,0);
    for(i=0;i<m->A->size1;i++)
        gsl_matrix_set(m->statedata,i,0,gsl_matrix_get(X0,i,0));
    for(i=0;i<m->C->size1;i++)
        gsl_matrix_set(m->outputdata,i,0,gsl_matrix_get(y0,i,0));

    m->Ts=Ts;
    m->Ndatapoints=Ndatapoints;
    printf("In InitSimModel:\n");
    printf("Bd");
    print2scr(Bd);
    printf("m->Bd");
    print2scr(m->Bd);
    gsl_matrix_memcpy(m->Bd,Bd);
    //m->Bd=Bd;

}

void ModelStep(Model *m,int step_time,gsl_matrix *u,double dist)
{

    double temp;

    int j;
    m->currxdata=MatAdd2(MatMul2(m->A,m->currxdata),MatMul2(m->B,u));
    for(j=0;j<m->A->size1;j++)
    {

         gsl_matrix_set(m->currxdata,j,0,gsl_matrix_get(m->currxdata,j,0)+(gsl_matrix_get(m->Bd,j,0)*dist));
    }

    m->currydata=MatMul2(m->C,m->currxdata);
    for(j=0;j<m->A->size1;j++)
        gsl_matrix_set(m->statedata,j,step_time,gsl_matrix_get(m->currxdata,j,0));
    for(j=0;j<m->C->size1;j++)
        gsl_matrix_set(m->outputdata,j,step_time,gsl_matrix_get(m->currydata,j,0));
    for(j=0;j<m->B->size2;j++)
        gsl_matrix_set(m->inputdata,j,step_time,gsl_matrix_get(u,j,0));

}

void printModeldata(Model *m,int s,char *filename)
{
    ///s==1 prints states //s==0 prints output s==2 prints output
    gsl_matrix *plotdata;
    if (s==0)
    {
        plotdata=gsl_matrix_alloc(m->statedata->size1,m->statedata->size2);
        plotdata=m->statedata;

    }
    else if(s==1)
    {
        plotdata=gsl_matrix_alloc(m->outputdata->size1,m->outputdata->size2);
        plotdata=m->outputdata;
    }
    else
    {
        plotdata=gsl_matrix_alloc(m->inputdata->size1,m->outputdata->size2);
        plotdata=m->inputdata;
    }
   FILE *fp;
    int rows;
    int cols;
    int i,j; //counters for loops
    rows=plotdata->size1;
    cols=plotdata->size2;
    fp=fopen(filename,"w");
    for(j=0;j<cols;j++)
    {
        if(j!=0)
        fprintf(fp,"\n");

        fprintf(fp,"%f ",j*m->Ts);
        for(i=0;i<rows;i++)
        {
            fprintf(fp,"%f ",gsl_matrix_get(plotdata,i,j));
        }
    }
    fclose(fp);
gsl_matrix_free(plotdata);
}

void ReadJac(Model *m,char *s) ///once you read a jac file it initialises model
{
int Ns,Nu,Ny; //number states, number inputs, number outputs
int i,j;
char buffer[50];
double t1,t2;
double **a;  //store matrix A
double **b; //store matrix B
double **c; //store matrix C
double **d; //store matrix D
double **x0; //initial conditions


//matrices A,B,C,D to store values.

gsl_matrix *A;
gsl_matrix *B;
gsl_matrix *C;
gsl_matrix *D;
gsl_matrix *X0;

FILE *jac;

jac=fopen(s,"r");
fscanf(jac,"%d%d%d",&Ns,&Nu,&Ny);
printf("\nValues read from file are %d %d %d",Ns,Nu,Ny);

a=malloc(Ns*sizeof(double *));
b=malloc(Ns*sizeof(double *));
c=malloc(Ny*sizeof(double *));
d=malloc(Ny*sizeof(double *));
x0=malloc(Ns*sizeof(double));

for(i=0;i<Ns;i++)
{
a[i]=malloc(Ns*sizeof(double));
b[i]=malloc(Nu*sizeof(double));
x0[i]=malloc(sizeof(double));
}

for(i=0;i<Ny;i++)
{
c[i]=malloc(Ns*sizeof(double));
d[i]=malloc(Nu*sizeof(double));
}
fgets(buffer,50,jac); //move to next line

//Read Matrix A
for(i=0;i<Ns;i++)
{
	printf("\n");
	for(j=0;j<Ns;j++)
	{
		fscanf(jac,"%lf",&a[i][j]);
		printf(" %lf",a[i][j]);
	}
}
//Read Matrix B
for(i=0;i<Ns;i++)
{
	printf("\n");
	for(j=0;j<Nu;j++)
	{
		fscanf(jac,"%lf",&b[i][j]);
		printf(" %lf",b[i][j]);
	}
}
//Read Matrix c
printf("\n C matrix");
for(i=0;i<Ny;i++)
{
	printf("\n");
	for(j=0;j<Ns;j++)
	{
		fscanf(jac,"%lf",&c[i][j]);
		printf(" %lf",c[i][j]);
	}
}
printf("\n D matrix");
//Read Matrix d
for(i=0;i<Ny;i++)
{
	printf("\n");
	for(j=0;j<Nu;j++)
	{
		fscanf(jac,"%lf",&d[i][j]);
		printf(" %lf",d[i][j]);
	}
}


//Read x0 initial conditions
printf("\n Initial Conditions");
for(i=0;i<Ns;i++)
{
fscanf(jac,"%lf",&x0[i][0]);
printf("\n %lf",x0[i][0]);
}
fclose(jac);
//allocation for matrices.
A=gsl_matrix_alloc(Ns,Ns);
B=gsl_matrix_alloc(Ns,Nu);
C=gsl_matrix_alloc(Ny,Ns);
D=gsl_matrix_alloc(Ny,Nu);
X0=gsl_matrix_alloc(Ns,1);
//copying values to matrices
assign_MatP(A,a);
assign_MatP(B,b);
assign_MatP(C,c);
assign_MatP(D,d);
assign_MatP(X0,x0);

printf("\n Matrix A \n");
print2scr(A);

printf("\n Matrix X0 \n");
print2scr(X0);
InitModel(m,A,B,C,D,X0);

gsl_matrix_free(A);
gsl_matrix_free(B);
gsl_matrix_free(C);
gsl_matrix_free(D);
gsl_matrix_free(X0);

for(i=0;i<Ns;i++)
{
        free(a[i]);
        free(b[i]);
        free(x0[i]);
}

free(a);
free(b);
free(x0);

for(i=0;i<Ny;i++)
{
     free(c[i]);
     free(d[i]);
}
free(c);
free(d);


}

//discretizes model
void discretize_model(Model *mc,Model *md,double Ts)
{
int Ns=mc->A->size1;
int Nu=mc->B->size2;
int Ny=mc->C->size1;

md->A=gsl_matrix_alloc(Ns,Ns);
md->B=gsl_matrix_alloc(Ns,Nu);
md->C=gsl_matrix_alloc(Ny,Ns);
md->D=gsl_matrix_alloc(Ny,Nu);
md->X0=gsl_matrix_alloc(Ns,1);
md->Bd=gsl_matrix_alloc(Ns,1);

gsl_matrix_memcpy(md->C,mc->C);
gsl_matrix_memcpy(md->D,mc->D);
gsl_matrix_memcpy(md->X0,mc->X0);


//We need 1 intermediate matrices Aint=Ac*Ts
gsl_matrix *Aint = gsl_matrix_alloc (Ns,Ns ); //allocate memory for matrix

//Set the values of identity matrix
gsl_matrix_set_identity (md->A); //Ad=Id
gsl_matrix_memcpy(Aint,mc->A);//Aint=Ac
gsl_matrix_scale (Aint,Ts); //Aint=Aint*Ts
gsl_matrix_add(md->A,Aint);//Ad=Id+Aint*Ts

gsl_matrix_memcpy(md->B,mc->B);//Bd=Bc
gsl_matrix_scale (md->B,Ts); //Bd=Bc*Ts

gsl_matrix_free(Aint);

}

void InitSimModelAme(Model *m,gsl_matrix *X0,double Ts,double tsim)
{
    int i,j;
    long Ndatapoints=tsim/Ts+1; //number of datapoints to be simulated
  //  gsl_matrix *y0=gsl_matrix_alloc(m->C->size1,1);
   // m->statedata=gsl_matrix_alloc(m->A->size1,Ndatapoints);
  //  m->outputdata=gsl_matrix_alloc(m->C->size1,Ndatapoints);
  //  m->currxdata=gsl_matrix_alloc(m->A->size1,1);
 //   m->currydata=gsl_matrix_alloc(m->C->size1,1);
    m->X0=X0;
   // m->currxdata=X0;
   // y0=MatMul2(m->C,X0);
   // m->currydata=y0;
    //gsl_matrix_set_all(m->statedata,0);
   // gsl_matrix_set_all(m->outputdata,0);
   // for(i=0;i<m->A->size1;i++)
   //     gsl_matrix_set(m->statedata,i,0,gsl_matrix_get(X0,i,0));
   // for(i=0;i<m->C->size1;i++)
        //gsl_matrix_set(m->outputdata,i,0,gsl_matrix_get(y0,i,0));

    m->Ts=Ts;
    m->Ndatapoints=Ndatapoints;

}
void LoadDoubles(Model *m,double *a,double *b, double *c, double *d, double *x0,int Ns,int Nu,int Ny)
{
m->A=gsl_matrix_alloc(Ns,Ns);
m->B=gsl_matrix_alloc(Ns,Nu);
m->C=gsl_matrix_alloc(Ny,Ns);
m->D=gsl_matrix_alloc(Ny,Nu);
m->X0=gsl_matrix_alloc(Ns,1);

assign_Mat(m->A,a);
assign_Mat(m->B,b);
assign_Mat(m->C,c);
assign_Mat(m->D,d);
assign_Mat(m->X0,x0);

}
void DiscrModel(Model *src,Model *dest, double Ts)
{
    gsl_matrix *expA;


    int Nu,Ns,Ny,i,j;

    Nu=src->B->size2;
    Ny=src->C->size1;
    Ns=src->A->size1;

    dest->A=gsl_matrix_alloc(Ns,Ns);
    dest->B=gsl_matrix_alloc(Ns,Nu);
    dest->C=gsl_matrix_alloc(Ny,Ns);
    dest->D=gsl_matrix_alloc(Ny,Nu);
    dest->X0=gsl_matrix_alloc(Ns,1);
    dest->Bd=gsl_matrix_alloc(Ns,1);


    expA=gsl_matrix_alloc(Ns+Nu,Ns+Nu);
    gsl_matrix_set_zero(expA);

    for(i=0;i<Ns;i++)
    {
        for(j=0;j<Ns+Nu;j++)
        {
             if(j<Ns)
            gsl_matrix_set(expA,i,j,gsl_matrix_get(src->A,i,j)*Ts);
            else
             gsl_matrix_set(expA,i,j,gsl_matrix_get(src->B,i,j-Ns)*Ts);

        }


    }

    printf("Before Discretization:\n");
    print2scr(expA);
    gsl_linalg_exponential_ss(expA,expA,1);
    printf("After Discretization:\n");
    print2scr(expA);
    for(i=0;i<Ns;i++)
    {
        for(j=0;j<Ns+Nu;j++)
        {
             if(j<Ns)
            gsl_matrix_set(dest->A,i,j,gsl_matrix_get(expA,i,j));
            else
             gsl_matrix_set(dest->B,i,j-Ns,gsl_matrix_get(expA,i,j));

        }


    }

    gsl_matrix_memcpy(dest->C,src->C);
    gsl_matrix_memcpy(dest->D,src->D);
    gsl_matrix_memcpy(dest->X0,src->X0);
    gsl_matrix_memcpy(dest->Bd,src->Bd);

    dest->Ts=Ts;

    gsl_matrix_free(expA);

}

/// discretizes input disturbance
/// takes in the model and continous form of Bd--Bdc
/// returns Bdk discrete form
///num Disturbances
void DiscrDisturbance(Model *cont,double* Bdc,double* Bdk,int numDisturbances,double Ts)
{




    int Ns,Nu,Ny,i,j;

    gsl_matrix *expA;

    Ns=cont->B->size1;
    Nu=cont->B->size2;

    expA=gsl_matrix_alloc(Ns+numDisturbances,Ns+numDisturbances);

    gsl_matrix_set_zero(expA);

    ///we need matrix of the form
    /// [A B]
    /// [0 0]

    ///TO BE DELELTED///
    ///VALUES PASSED//





    for (i=0;i<Ns;i++)
        for(j=0;j<Ns+numDisturbances;j++)
    {
        if(j<Ns)
            gsl_matrix_set(expA,i,j,gsl_matrix_get(cont->A,i,j)*Ts);
        else
            gsl_matrix_set(expA,i,j,Ts*Bdc[i+(j-Ns)*Ns]);
    }



    /// Discretize funtion
    gsl_linalg_exponential_ss(expA,expA,1);


    print2scr(expA);
    for(i=0;i<Ns;i++)
        for(j=Ns;j<Ns+numDisturbances;j++)
        {
             Bdk[i+(j-Ns)*Ns]=gsl_matrix_get(expA,i,j);
             printf("bdk[%d]=%f\n",i+(j-Ns)*Ns,Bdk[i+(j-Ns)*Ns]);
        }







    gsl_matrix_free(expA);

}

