#include "helper.h"

void print2File(gsl_matrix *m,double *time_v,char *filename)
{

    FILE *fp;
    int rows;
    int cols;
    int i,j; //counters for loops
    rows=m->size1;
    cols=m->size2;
    fp=fopen(filename,"w");
    for(i=0;i<rows;i++)
    {
        if(i!=0)
        fprintf(fp,"\n");

        fprintf(fp,"%f ",time_v[i]);
        for(j=0;j<cols;j++)
        {
            fprintf(fp,"%f ",gsl_matrix_get(m,i,j));
        }
    }
    fclose(fp);


}

void print2scr(gsl_matrix *m)
{
    int rows;
    int cols;
    int i,j; //counters for loops
    rows=m->size1;
    cols=m->size2;
	printf("\n");
    for(i=0;i<rows;i++)
    {
        if(i!=0)
        printf("\n");
        for(j=0;j<cols;j++)
        {
            printf("%f ",gsl_matrix_get(m,i,j));
        }
    }
	printf("\n");
}

// Takes a Matrix, inverts it and put the result in the A-1
void matrixInvert(gsl_matrix *A,gsl_matrix *Ainv)
{
int nRows=A->size1;
int nCols=A->size2;
int s; //signum for LU decomposition

// Define all the used matrices
	gsl_permutation *perm = gsl_permutation_alloc(nRows);

	// Make LU decomposition of matrix A
	gsl_linalg_LU_decomp(A,perm,&s);

	// Invert the matrix A
	gsl_linalg_LU_invert (A,perm,Ainv);


}

gsl_matrix *MatInv2(gsl_matrix *m) //returns m^{-1}
{
int nRows=m->size1;
int nCols=m->size2;
int s; //signum for LU decomposition
gsl_matrix *minv=gsl_matrix_alloc(nRows,nCols);

// Define all the used matrices
	gsl_permutation *perm = gsl_permutation_alloc(nRows);

	// Make LU decomposition of matrix A
	gsl_linalg_LU_decomp(m,perm,&s);

	// Invert the matrix A
	gsl_linalg_LU_invert (m,perm,minv);
	return minv;
}


void MatMul(gsl_matrix *m, gsl_matrix *n,gsl_matrix *p ) //multiplies two matrices  p<-m*n
{
//   Matrix multiplication:   gsl_blas_dgemm(Operator A , Operator B , alpha , matrix A , matrix B , beta , matrix C);
//                            Compute:   p = (alpha)mn + (beta)p
//                            Operator:   CblasNoTrans,CBlasTrans,CblasConjTrans
//                            alpha,beta:   scaling value
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m,n,0.0,p);
}

gsl_matrix* MatTrans(gsl_matrix *n) //return n^T
{
gsl_matrix *m=gsl_matrix_alloc(n->size2,n->size1);
gsl_matrix_transpose_memcpy(m,n);
return m;
}

gsl_matrix* MatMul2(gsl_matrix *m, gsl_matrix *n ) //returns m*n
{
int rows,cols;
gsl_matrix *p2;
rows=m->size2;
cols=n->size1;
if(rows!=cols)
printf("Cannot multiply");
p2=gsl_matrix_alloc(m->size1,n->size2);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m,n,0.0,p2);
return p2;
}

gsl_matrix* MatMulrec(gsl_matrix *m,int p) //returns m^n
{

gsl_matrix *temp=gsl_matrix_alloc(m->size1,m->size2);
gsl_matrix *p2=gsl_matrix_alloc(m->size1,m->size2);
if(p>0)
{
 gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m,MatMulrec(m,p-1),0.0,p2);
 return p2;
}
else
{
  gsl_matrix_set_identity(p2);
  return p2;
}

}

gsl_matrix* MatAdd2(gsl_matrix *m, gsl_matrix *n ) //returns m+n
{
gsl_matrix* add=gsl_matrix_alloc(m->size1,m->size2);
gsl_matrix_memcpy(add,m);//p=m
gsl_matrix_add(add,n);//p=p+n
return add;
}

gsl_matrix* MatSub2(gsl_matrix *m, gsl_matrix *n )//returns n-m
{
gsl_matrix* p=gsl_matrix_alloc(m->size1,m->size2);
gsl_matrix_memcpy(p,m);//p=m
gsl_matrix_sub(p,n);//p=n-p
return p;
}

void MatAdd(gsl_matrix *m, gsl_matrix *n,gsl_matrix *p) //p=m+n
{
gsl_matrix_memcpy(p,m);//p=m
gsl_matrix_add(p,n);//p=p+n
}

void printsizeMat(gsl_matrix *m,char *s)
{
printf("Size of %s:%dx%d\n",s,m->size1,m->size2);
}

void assign_Mat(gsl_matrix *m,double *val)
{
// Takes a double array and assigns the value to a matrix
    int rows=m->size1;
    int cols=m->size2;
    int i,j,k=0;
    for (i = 0,k=0; i <rows; i++)
        for (j = 0; j <cols; j++,k++)
            gsl_matrix_set (m, i, j, val[k]);
}

void assign_MatP(gsl_matrix *m, double **val)
{
//Takes a pointer to pointer (multidimensional array) and assigns it to a matrix
int rows=m->size1;
    int cols=m->size2;
    int i,j;
    for (i = 0; i <rows; i++)
        for (j = 0; j <cols; j++)
            gsl_matrix_set (m, i, j, val[i][j]);
}

void print2FileMat(gsl_matrix *m,FILE *filename)
{
    int rows=m->size1;
    int cols=m->size2;
    int i,j;
    for (i=0;i<rows;i++)
    {
        for(j=0;j<cols;j++)
            fprintf(filename,"%f ",gsl_matrix_get(m,i,j));
        fprintf(filename,"\n");
    }
}

void assignMatMat(gsl_matrix *m,gsl_matrix *n)
{
    int i,j;

    for(i=0;i<m->size1;i++)
        for(j=0;j<m->size2;j++)
            gsl_matrix_set(m,i,j,gsl_matrix_get(n,i,j));

}
