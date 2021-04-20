/****************************************************************
    chol: functions involving the Cholesky decomposition of a positive
    definite matrix.
    Copyright (C) 1999 Alan R. Rogers

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    Alan R. Rogers, Department of Anthropology, University of Utah,
    Salt Lake City, UT 84112. rogers@anthro.utah.edu
****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mytypes.h"
#include "chol.h"
/*****************************************************************
 * These routines compute and manipulate the Cholesky decomposition of
 * a positive definite matrix.
 *   The Cholesky factor is a lower triangular matrix, L, such that
 * L * L' = A, where the "'" denotes matrix transposition.
 *
 * Alan R. Rogers, Dept. of Anthropology, University of Utah, S.L.C.,
 * UT 84108.    February 23, 1998
****************************************************************/
#if 0
/********** Prototypes **********************/
int   cholesky(real **A, int n);
real log_determinant(real **L, int n);
real dotprod(real *x, real *y, int n);
void  matmultLLT(real **L, real **B, int n);
int   col_L_solve(real **L, int n, real *b);
int   row_L_solve(real **L, int n, real *b);
real mahal(real *x, real **L, int n);
#endif

/*****************************************************************
 * Find the left Cholesky factor, L, of a positive definite
 * matrix, A.  The code here was modified from the gaxpy version of
 * the cholesky algorithm on p. 143 of Matrix Computations (2nd ed),
 * by Golub and Van Loan. 
 * 
 * On entry
 *     A     is an n X n matrix of reals as allocated by alloc2d().
 *     n     is the dimension of A.
 * On return
 *     A     contains in its lower triangle the left (lower) Cholesky
 *           factor, L, of A.  L is a lower triangular matrix such
 *           that L L' = A.
 * The routine returns 0 if successful, or 1 if A is not positive 
 * definite.
 *                   Alan R. Rogers Feb 23, 1998
 *****************************************************************/
#ifdef PIVOTING
int   cholesky(real **A, int n, int *rowpvt, int *colpvt)
#else
int cholesky(real **A, int n)
#endif
{
  int i,j;
  real root;

  for(i=0; i<n; i++)
  {
    if(i>0)
    {
      for(j=i; j<n; j++)
	A[j][i] -= dotprod(A[j], A[i], i);
    }
    if(A[i][i] <= 0.0)
      return(1);
    root = sqrt(A[i][i]);
    for(j=i; j<n; j++)
      A[j][i] /= root;
  }
  return(0);
}

/**************************************************************** 
Return a vector of ints that sorts the diagonal entries of
A in decreasing order, except that all non-finite values are
placed at the end.  Return the number of finite values.  
****************************************************************/
int *srt_diag_ndx=NULL;
real *srt_diag_val=NULL;
int srt_diag_vec_size=0;
int srt_diag(real **A, int *ndx, int n)
{
  int i, nfinite;
  
  /* allocate memory if necessary */
  if(srt_diag_vec_size < n)
  {
    if(srt_diag_ndx == NULL)
      free(srt_diag_ndx);
    if(srt_diag_val == NULL)
      free(srt_diag_val);
    srt_diag_vec_size = n;
    srt_diag_ndx = (int *) malloc(n * sizeof(int));
    srt_diag_val = (real *) malloc(n * sizeof(real));
    if(srt_diag_ndx == NULL || srt_diag_val == NULL)
    {
      fprintf(stderr,"\nsrt_diag: memory\n");
      exit(1);
    }
  }
  /* initialize */
  for(i=0; i<n; i++)
  {
    srt_diag_ndx[i] = i;
    srt_diag_val[i] = A[i][i];
  }
  qsort(srt_diag_ndx, (unsigned) n, sizeof(int),
	(int (*)(const void*, const void*)) srt_diag_compar);
  for(i=0; i<n; i++)
    ndx[i] = srt_diag_ndx[i];
  for(nfinite=0;
      nfinite < n && FINITE(srt_diag_val[ndx[nfinite]]);
      nfinite++)
    ;
  return(nfinite);
}

/* comparison for srt_diag */
int srt_diag_compar(int *i, int *j)
{
  if(FINITE(srt_diag_val[*i]) && !FINITE(srt_diag_val[*j]))
     return(-1);
  if(FINITE(srt_diag_val[*j]) && !FINITE(srt_diag_val[*i]))
    return(1);
  if(srt_diag_val[*i] > srt_diag_val[*j])
    return(-1);
  if(srt_diag_val[*i] < srt_diag_val[*j])
    return(1);
  return(0);
}

/* invert a permutation vector */
void invert_permutation(int *permutation, int *inv, int n)
{
  int i;
  for(i=0; i<n; i++)
    inv[permutation[i]] = i;
}

/* multiply a lower triangular matrix L by its transpose. */
/* Answer is returned in matrix B.                        */
void matmultLLT(real **L, real **B, int n)
{
  int i, j, k;

  for(i=0; i<n; i++)
    for(j=0; j<=i; j++)
    {
      B[i][j] = 0.0;
      for(k=0; k<=j; k++)
	B[i][j] += L[i][k] * L[j][k];

      if(j < i)
	B[j][i] = B[i][j];
    }
}

/* Use cholesky decomposition to get log determinant */
/* L should already be in Cholesky form */
real log_determinant(real **L, int n)
{
  real det;
  int i;

  det = log(L[0][0]);
  for(i=1; i<n; i++)
    det += log(L[i][i]);
  return(2.0*det);
}

/* return inner product of two vectors */
real dotprod(real *x, real *y, int n)
{
  int i;
  real z=0.0;

  for(i=0; i<n; i++)
    z += x[i]*y[i];

  return(z);
}
/****************************************************************
 * Solve
 *                          b = L * x
 *
 * where b and x are column vectors, and L is lower triangular with
 * diagonal entries not necessarily equal to 1.0.  The diagonal entries
 * of L *ARE* referenced.  On return b contains the solution vector x;
 ****************************************************************/
int col_L_solve(real **L, int n, real *b)
{
    int i;
    register int j;
    register real x;
    
    for(i=0; i<n; i++)
    {
      if(L[i][i] <= 0.0)
	return(1);
      b[i] = x = b[i]/L[i][i];
      for(j=i+1; j<n; j++)
	b[j] -= x * L[j][i];
    }
    return(0);
}
/*****************************************************************
 * Solve 
 *                            b' = x' * L
 * where b' and x' are row vectors, and L is lower triangular with diagonal
 * entries not necessarily equal to 1.0.  The diagonal entries of L *ARE*
 * referenced.  On return b contains the solution vector x;
 *****************************************************************/
int row_L_solve(real **L, int n, real *b)
{
    int i;
    register int j;
    register real x;
    
    for(i=n-1; i>=0; i--)
    {
      if(L[i][i] <= 0.0)
	return(1);
      b[i] = x = b[i]/L[i][i];
      for(j=0; j<i; j++)
	b[j] -= x * L[i][j];
    }
    return(0);
}

/****************************************************************
Mahalanobis distance: x'*inv(A)*x.
Should be used after a call to cholesky.
ON INPUT:
x	a vector of n reals
L	an nXn matrix of reals, as allocated by alloc2d.  The lower
        triangle of a should contain the left Cholesky factor of A.
n	dimension of square matrix L.
ON RETURN:
Function returns x'*inv(a)*x if successful.  This value is always positive,
	and negative returns indicate an error.
L	is unchanged.
x	is destroyed.
ERROR RETURN:
Returned value is < 0 if Cholesky factor is singular.

ALGORITHM:

The Cholesky decomposition turns
                  -1
        q = x' * A  * x.
 
into
                  -1   -1         2
        q = x' * L'  * L * x = |z|
where
             -1
        z = L * x
 
To find z we need only solve the triangular system
 
        L * z = x.
*****************************************************************/
real mahal(real *x, real **L, int n)
{
  int i;
  real dist;

  /* solve triangular system */
  if(col_L_solve(L, n, x))
  {
    fprintf(stderr,"\nmahal: cholesky factor is singular");
    exit(1);
  }

  /* get squared norm of solution vector */
  dist = 0.0;
  for(i = 0; i<n; i++)
    dist += x[i]*x[i] ;
  
  return(dist);
}

/***************************************************************

On input L is the left cholesky factor of a n X n matrix of
reals, called "a".  In other words, L is a lower triangular matrix
such that L * L' = a.   

On return, L contains the inverse of the input matrix L, and
s contains the inverse of "a", calculated as inv(L)' * inv(L).

Sources: Numerical Recipes, from p 312 of the Linpack Users' Guide,
1979 edition, and the Linpack routine spodi.f. 
***************************************************************/
void cholinv(real **L, real **s, int n)
{
  int i, j, k;
  real sum;

  /* invert L in place */
  for(i=0; i<n; i++)
  {
    L[i][i] = 1.0/L[i][i];
    for(j=i+1; j<n; j++)
    {
      sum = 0.0;
      for(k=i; k<j; k++)
	sum -= L[j][k] * L[k][i];
      L[j][i] = sum/L[j][j];
    }
  }

  /* form inv(A) = inv(L)' * inv(L) */
  for(i=0; i<n; i++)
  {
    for(j=0; j<=i; j++)
    {
      s[i][j] = 0.0;
      for(k=i; k<n; k++)
	s[i][j] += L[k][i] * L[k][j];
    }
    for(j=i+1; j<n; j++)
    {
      s[i][j] = 0.0;
      for(k=j; k<n; k++)
	s[i][j] += L[k][i] * L[k][j];
    }
  }
}
