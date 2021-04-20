#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "mytypes.h"
#include "bone.h"
#ifndef NDEBUG
int symmetric(real **c, int dim);
#endif

extern int rdim, npart, nagent, **lump;
extern real **m, ***F, *sensitivity;
extern int *y;

/* log likelihood is multivariate normal */
real lnL(real *p)
{
  int i, j;
  real  lndet, mh, rval, u;
  static real **c=NULL;
  static real **c2=NULL;
  static real *x=NULL;
  static real *x2=NULL;

  if(c==NULL)
  {
    x = (real *) mustalloc(npart * sizeof(real));
    x2 = (real *) mustalloc(rdim * sizeof(real));
    c = (real **) alloc2d(npart, npart, sizeof(real));
    c2 = (real **) alloc2d(rdim, rdim, sizeof(real));
    if(c==NULL || c2==NULL)
      error("lnL: alloc2d");
  }

  /* get vector x = y - E{y} */
  for(i=0; i<npart; i++)
  {
#if 1
    /* new code */
    u = 0.0;
    for(j=0; j<nagent; j++)
      u += alpha[j] * agent[j]->m[i];
    u *= LAMBDA;
#else    
    /* old code */
    u = LAMBDA*(ALPHA*agent[0]->m[i] + (1.0-ALPHA)*agent[1]->m[i]);
#endif
#if ATTRITION    
    u *= exp(-BETA*bones->sensitivity[i]); 
#endif
    x[i] = y[i] - u;
  }

#if 0
  privec(y, npart, "y");
  prfvec(x, npart, "y - mu");
#endif

  /* get covariance matrix */
  for(i=0; i<npart; i++)
  {
    c[i][i] = LAMBDA*(ALPHA*agent[0]->F[i][i]
		      + (1.0-ALPHA)*agent[1]->F[i][i]);
#if ATTRITION    
    c[i][i] *= exp(-2.0*BETA*bones->sensitivity[i]);
#endif
    for(j=0; j<i; j++)
    {
      c[i][j] = LAMBDA*(ALPHA*agent[0]->F[i][j]
			+ (1.0-ALPHA)*agent[1]->F[i][j]);
#if ATTRITION      
      c[i][j] *= exp(-BETA*(bones->sensitivity[i]
			    + bones->sensitivity[j]));
#endif
      c[j][i] = c[i][j];
    }
  }
  assert(symmetric(c, npart));
#if 0
  prfmat(c, rdim, rdim, "C");
#endif

  /* reduce dimension */
  vTx(x, lump, x2, npart, rdim);
  xTax(lump, c, c2, npart, rdim);

  assert(symmetric(c2, rdim));
#if 0
  prfmat(c2, rdim, rdim, "reduced C");
#endif

  /* calculate cholesky decomposition */
  if(cholesky(c, rdim) != 0)
  {
    fflush(stdout);
    fprintf(stderr,"\nlnL: reduced C is not PD\n");
    exit(1);
    rval = -HUGE_VAL;
  }else
  {
    lndet = log_determinant(c, rdim);  /* determinant */
    mh = mahal(x, c, rdim);           /* mahalanobis distance */
    rval = -0.5*(lndet + mh);
  }
#if 0
  printf("\nlndet=%g mahal=%g", lndet, mh);
  printf("\nlnL(l=%g, a=%g) = %g", LAMBDA, ALPHA, rval);
#endif  

  return(rval);
}

#ifndef NDEBUG
/* test a square matrix of reals for symmetry */
int symmetric(real **c, int dim)
{
  int i, j;
  real x;
  for(i=0; i<dim; i++)
    for(j=0; j<i; j++)
    {
      if(fabs(c[i][j]) >= fabs(c[j][i]))
	x = fabs(c[i][j]);
      else
	x = fabs(c[j][i]);
      if(fabs(c[i][j]-c[j][i]) > x*0.0001)
      {
	printf("\nasymmetric: mat[%d][%d]=%.20g mat[%d][%d]=%.20g\n",
	       i, j, c[i][j], j, i, c[j][i]);
	return(0);    /* asymmetric */
      }
    }
  return(1);          /* symmetric */
}
#endif

/* make a lump matrix */
int **lumpmat(int a, int b, int dim)
{
  int row,col;
  int **mat;

  mat = (int**) alloc2d(dim, dim-1, sizeof(int));
  if(mat==NULL)
    return(mat);

  /* make a < b */ 
  if(b<a)
  {
    col=b;
    b=a;
    a=col;
  }
  assert(a<b);

  /* rows 0..(b-1) are an identity matrix */
  for(row=0; row<b; row++)
  {
    for(col=0; col<dim-1; col++)
      mat[row][col] = (row==col ? 1 : 0);
  }

  /* row b has 1 in column a */
  for(col=0; col<dim-1; col++)
    mat[b][col] = (col==a ? 1 : 0);

  /* rows b+1..(dim-1) are an offset identity matrix */
  for(row=b+1; row < dim; row++)
    for(col=0; col<dim-1; col++)
      mat[row][col] = (col == row-1 ? 1 : 0);

  return(mat);
}

/* make an identity matrix of ints */
int **id_mat(int dim)
{
  int row,col;
  int **mat;

  mat = (int**) alloc2d(dim, dim, sizeof(int));
  if(mat==NULL)
    return(mat);

  for(row=0; row<dim; row++)
    for(col=0; col<dim; col++)
      mat[row][col] = (row==col ? 1 : 0);

  return(mat);
}

/**************************************************************** 
Form b = x' a x

x : matrix of ndim X kdim ints
a : symmetric matrix of ndim X ndim reals
b : a kdim X kdim real matrix into which the answer will be put

****************************************************************/
void xTax(int **x, real **a, real **b, int ndim, int kdim)
{
  int i,j,k,l;

  assert(b!=NULL);

  for(i=0; i<kdim; i++)
    for(j=0; j<kdim; j++)
    {
      b[i][j] = 0.0;
      for(k=0; k<ndim; k++)
	for(l=0; l<ndim; l++)
	  b[i][j] += x[k][i] * a[k][l] * x[l][j];
    }
}

/**************************************************************** 
  Form b = v' * x

x : matrix of ndim X kdim ints
v : vector of ndim reals
b : resulting vector of kdim reals  
****************************************************************/  
void vTx(real *v, int **x, real *b, int ndim, int kdim)
{
  int i, j;

  for(i=0; i<kdim; i++)
  {
    b[i] = 0.0;
    for(j=0; j<ndim; j++)
      b[i] += v[j] * x[j][i];
  }
}

/* lump two skeletal parts */
void dolump(int a, int b)
{
  int i, j, k, **lmat, **lump2;

  lmat = lumpmat(a,b,rdim);
  if(lump==NULL)
  {
    lump = lmat;
    return;
  }else
  {
    /* allocate new lump matrix */
    lump2 = (int **) alloc2d(npart,rdim-1, sizeof(int));
    if(lump2==NULL)
    {
      fprintf(stderr,"\ndolump:alloc2d\n");
      exit(1);
    }
    /* matrix multiplication: lump2 = lump * lmat */
    for(i=0; i<npart; i++)
      for(j=0; j<rdim-1; j++)
      {
	lump2[i][j] = 0;
	for(k=0; k<rdim; k++)
	  lump2[i][j] += lump[i][k] * lmat[k][j];
      }
    free2d((void **) lump);
    free2d((void **) lmat);
    lump = lump2;
  }
  rdim -= 1;
}

/* calculate the mean F matrix */
void getFmean(real **Fmean, real ***F, int dim, real *alpha, int nagent)
{
  int i, j, k;
  real sumalpha;

  for(j=0; j<dim; j++)
    for(k=0; k<dim; k++)
    {
      Fmean[j][k] = sumalpha = 0.0;
      for(i=0; i<nagent-1; i++)
      {
	sumalpha += alpha[i];
	Fmean[j][k] += alpha[i]*F[i][j][k];
      }
      Fmean[j][k] += (1.0 - sumalpha)*F[nagent-1][j][k];
    }
}
