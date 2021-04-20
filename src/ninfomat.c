/****************************************************************
    ninfomat: Calculate the information matrix by numerical differentiation
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
#include "misc.h"
#include "chol.h"
#include "io.h"
#include "ninfomat.h"
/****************************************************************
This file handles calculation of the information matrix by numerical
differentiation.  
****************************************************************/
void dummy(double x);

/* the following externals are demanded by io.c but are not used */
int dim=-1;
int *live=NULL;
real *density=NULL;


/* This dummy function is really needed */
void dummy(double x){}

/****************************************************************
Return values:
  0  success
  1  information matrix was not PD
****************************************************************/
int ninfomat(real **C, real *p0, int dimension, real (*f)(real *pvec), int attrition)
{
  real *p;
  real *add;                  /* vector of additions */
  real *sub;                  /* vector of subtractions */
  real **infomat;             /* information matrix */
  real **srtmat;              /* sorted matrix */
  real eps, base_eps;         /* scale of error */
  real x, y;
  int i, j, first_alpha, rdimension;
  int *ndx, *indx;

  p = (real *) mustalloc( dimension * sizeof(real));
  add = (real *) mustalloc( dimension * sizeof(real));
  sub = (real *) mustalloc( dimension * sizeof(real));
  ndx = (int *) mustalloc( dimension * sizeof(int));
  indx = (int *) mustalloc( dimension * sizeof(int));
  infomat = (real **) alloc2d(dimension, dimension, sizeof(real));
  srtmat = (real **) alloc2d(dimension, dimension, sizeof(real));
  if(infomat == NULL || srtmat==NULL)
    error("ninfomat: memory");
  
  base_eps = 0.001;
  first_alpha = (attrition ? 2 : 1);  /* index of first alpha */

  /****************************************************************
  Set magnitude of deviations, making sure that p+add and p differ
  by exact machine numbers and ditto for p-sub and p.

  I check the inequality constraints below in a manner slightly 
  different than in lnL.  lnL returns -Inf unless sum(alpha[i]) = 1.
  Here, I only check that 0 <= alpha[i] <= 1.  This is necessary
  because when the last alpha is 0, any increase in any other alpha
  makes the sum exceed 1.  Such increases are necessary for 
  calculating numerical derivatives.
  ****************************************************************/
  for(i=0; i<dimension; i++)
  {
#if 1    
    if(p0[i] < base_eps)
    {
      eps = p0[i];
    }else
    {
      eps = base_eps;
    }
#else
    eps = base_eps;
#endif
    add[i] = eps;
    if(p0[i] > 1.0)
      add[i] *= p0[i];
    sub[i] = add[i];

    x = p0[i] + add[i];
    dummy(x);           /* circumvent optimizing compilers */
    add[i] = x - p0[i];

    x = p0[i] - sub[i];
    dummy(x);
    sub[i] = p0[i] - x;

    p[i] = p0[i];
  }

  for(i=0; i<dimension; i++)
  {
    p[i]  = p0[i] + add[i];
    /* make sure alpha[i] <= 1 */
    if(i >= first_alpha && (p[i]*(1.0+FTOL) > 1.0 || p0[i]*(1.0+FTOL) > 1.0))
      infomat[i][i] = MISSING;
    else
    {
      x = f(p);
      p[i] = p0[i] - sub[i];
      /* make sure p[i] >= 0.0 */
      if(p[i] < FTOL || p0[i] < FTOL)
	infomat[i][i] = MISSING;
      else
      {
	x += f(p);
	x -= 2.0 * f(p0);
	y = add[i] + sub[i];
	infomat[i][i] = -4.0 * x / (y*y);
      }
    }
    p[i] = p0[i];
    
    for(j=0; j<i; j++)
    {
      if(!FINITE(infomat[j][j]) || !FINITE(infomat[i][i]))
      {
	infomat[j][i] = infomat[i][j] = MISSING;
	continue;
      }
      /* f(x+h, y+h) */
      p[i] = p0[i] + add[i];
      p[j] = p0[j] + add[j];
      x = f(p);

      /* f(x-h, y+h) */
      p[i] = p0[i] - sub[i];
      y = f(p);

      /* f(x-h, y-h) */
      p[j] = p0[j] - sub[j];
      x += f(p);

      /* f(x+h, y-h) */
      p[i] = p0[i] + add[i];
      y += f(p);

      p[i] = p0[i];
      p[j] = p0[j];

      x -= y;
      y = add[i] + sub[i];
      y *= add[j] + sub[j];
      infomat[i][j] = -x / y;
      infomat[j][i] = infomat[i][j];
    }
  }

  /****************************************************************
  Find ndx vector that sorts diagonal entries.  rdimension is the
  number of finite diagonal entries.
  ****************************************************************/
  rdimension = srt_diag(infomat, ndx, dimension);
  invert_permutation(ndx, indx, dimension);

  /* Use ndx to sort infomat, producing srtmat */
  for(i=0; i<rdimension; i++)
    for(j=0; j<rdimension; j++)
      srtmat[i][j] = infomat[ndx[i]][ndx[j]];
#if 0
  prfmat(srtmat, rdimension, rdimension, "srtmat");
#endif

  /* cholesky factorization of srtmat */
  i = cholesky(srtmat, rdimension);

  if(i!= 0)                  /* error return */
  {
    printf("\n#Warning from ninfomat: Information matrix is not PD.");
    printf("  (cholesky returned %d)", i);
#if 1
    prfmat(infomat, dimension, dimension, "Information matrix");
    prfmat(srtmat, rdimension, rdimension, "Partially factorized matrix");
#endif
    free2d((void **) infomat, dimension);
    free2d((void **) srtmat, dimension);
    return(i);
  }

  /* inverse of sorted information matrix */
  cholinv(srtmat, infomat, rdimension);

  /* put NaNs into appropriate parts of infomat */
  for(i=0; i<rdimension; i++)
    for(j=rdimension; j<dimension; j++)
      infomat[i][j] = 0.0/0.0;
  for(i=rdimension; i<dimension; i++)
    for(j=0; j<dimension; j++)
      infomat[i][j] = 0.0/0.0;

  /* put rows and columns back into original order */
  for(i=0; i<dimension; i++)
    for(j=0; j<dimension; j++)
      C[i][j] = infomat[indx[i]][indx[j]];
  free(p);
  free(add);
  free(sub);
  free(ndx);
  free(indx);
  free2d((void **) infomat, dimension);
  free2d((void **) srtmat, dimension);
  return(0);                 /* normal return */
}
