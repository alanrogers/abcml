/****************************************************************
    infomat: Calculate the information matrix
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
#include "mytypes.h"
#include "infomat.h"
/****************************************************************
Quadratic form of vector v, square matrix m, and vector u, all of
dimension dim: ans = v' * m * u
****************************************************************/
real quadratic_form(real *v, real **m, real *u, int dim)
{
  real ans=0.0;
  
  for(i=0; i<dim; i++)
    for(j=0; j<dim; j++)
      ans += v[i] * m[i][j] * u[j];

  return(ans);
}

/* 2nd partial of l w/ respect to kappa */
real l_kk(void)
{
  return(-0.5*(lDetC_kk + Q_kk));
}

/* 2nd mixed partial of l w/ respect to kappa and beta */
real l_kb(void)
{
  return(-0.5*(lDetC_kb + Q_kb));
}

/* 2nd mixed partial of l w/ respect to kappa and alpha[k] */
real l_ka(void)
{
  return(-0.5*(lDetC_ka + Q_ka));
}

/* 2nd partial of l w/ respect to beta */
real l_bb(void)
{
  return(-0.5*(lDetC_bb + Q_bb));
}

/* 2nd mixed partial of l w/ respect to beta and alpha[k] */
real l_ba(void)
{
  return(-0.5*(lDetC_ba + Q_ba));
}

/* 2nd mixed partial of l w/ respect to alpha_k and alpha_j */
real l_aa(void)
{
  return(-0.5*(lDetC_aa + Q_aa));
}

/* 2nd partial of log |C| w/ respect to kappa */
real lDetC_kk(void)
{
  return( -P/(kappa*kappa) );
}

/* 2nd mixed partial of log |C| w/ respect to kappa and beta */
real lDetC_kb(void)
{
  return( 0 );
}

/* 2nd mixed partial of log |C| w/ respect to kappa and alpha[k] */
real lDetC_ka(void)
{
  return( 0 );
}

/* 2nd partial of log |C| w/ respect to beta */
real lDetC_bb(void)
{
  return( 0 );
}

/* 2nd mixed partial of log |C| w/ respect to beta and alpha[k] */
real lDetC_ba(void)
{
  return( 0 );
}

/****************************************************************
 2nd mixed partial of log |C| w/ respect to alpha[k] and alpha[j]
 equals trace( A^{-1} F_j A^{-1} F_k ).  To solve it, first solve
 A X = F_j and A Y = F_k for X and Y.  Then form product X Y and
 take its trace.

 Is there a simpler algorithm?
****************************************************************/
real lDetC_aa(void)
{
}

/****************************************************************
 2nd partial of Q w/ respect to kappa is

(1/kappa^2) [ mu' Cinv mu + y' Cinv y ]
****************************************************************/
real Q_kk(void)
{
  real ans;

  ans = mahal(mu, C_chol, P) + mahal(y, C_chol, P);
  ans /= kappa*kappa;
  return(ans);
}

/****************************************************************
2nd mixed partial of Q w/ respect to kappa, beta.  
 ****************************************************************/
real Q_kb(void)
{
  int i, j;
  real ans, **mat;

  /* check diagonal entries of S */
  j=0;
  for(i=1; i<P; i++)
  {
    if(S[i][i] != S[i-1][i-1])
    {
      j=1;
      break;
    }
  }
  if(j==0)  /* early return: diagonal entries of S are equal */
    return(0.0);

  mat = alloc2d(P, P, sizeof(real));
  if(mat==NULL)
  {
    fprintf(stderr,"\nQ_kb: memory\n");
    exit(1);
  }

  /* form S*Cinv */
  for(i=0; i<P; i++)
  {
    for(j=0; j<P; j++)
      mat[i][j] = S[i][j];
    row_L_solve(C_chol, P, mat[i]);
  }

  /* form mat = S*Cinv - Cinv*S */
  for(i=0; i<P; i++)
  {
    mat[i][i] = 0;
    for(j=0; j<i; j++)
    {
      mat[i][j] -= mat[j][i];
      mat[j][i] = mat[i][j];
    }
  }

  /* set: ans = mu' * mat * mu - y' * mat * y */
  ans = quadratic_form(mu, mat, mu, P) - quadratic_form(y, mat, y, P);

  /* set: ans /= kappa */
  ans /= kappa;

  free2d(mat, P);
  return(ans);
}

/****************************************************************
2nd mixed partial of Q w/ respect to kappa, alpha_k:

2 * mu' * Cinv * exp(-beta*S) * m[k]
+ (y-mu)' Cinv exp(-beta*S) * F[k] exp(-beta*S) * Cinv (y+mu)
 ****************************************************************/
real Q_ka(void)
{
  real ans;
  
  /* mat1 = Cinv * exp(-beta*S)    */
  /* mat2 = mat1 * F[k]            */
  /* mat3 = mat2 * transpose(mat1) */
  /* v1 = y-mu                     */
  /* v2 = y+mu                     */
  /* ans = 2 * mu' * mat1 * m[k]   */
  /* ans += v1' * mat2 * v2        */

  return(ans);
}

/****************************************************************
Calculate information matrix.

On entry:
  imat   is an nparam X nparam matrix of reals, into which the 
         answer will be placed,
  C      is the covariance matrix, an PXP matrix of reals, 
  S      is a PXP diagonal matrix of reals, whose i'th diagonal 
         entry is s[i], the sensitivity of the i'th skeletal part,
  nparam is the number of parameters in the model,
  P      is the dimension of C and S
****************************************************************/
int infomat(real **imat, real **C, **S, int iparam, int P)
{
  int i, j, k;

  
}













