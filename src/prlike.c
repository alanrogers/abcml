/****************************************************************
    prlike: Print a transect through the likelihood function
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
#include <assert.h>
#include <math.h>
#include "mytypes.h"
#include "prlike.h"
#include "misc.h"
#include "likelihood.h"
/****************************************************************
Print a transect through the likelihood function.

On input:
  k   is an integer that determines which parameter to vary
  p   is the base parameter vector
  len is the length of p
  l   is low value of parameter k
  h   is high value of parameter k
  n   is the number of points along the transect

The function prints (x,y) pairs where  
  x   contains the entry with index k of vectors:
      p - l*inc, ..., p - inc, p, p + inc, ... , p + (n-l-1)*inc 
  y   contains the likelihoods of these parameter vectors
****************************************************************/
void ltransect(int k, real *p, int len, real l, real h, int n)
{
  extern int attrition;
  int i, nn;
  real x, y, *u=NULL, inc;

  /* allocate array */
  u = (real *) mustalloc(n*sizeof(real));

  assert(k < len);
  assert(l <= p[k]);
  assert(h >= p[k]);
  assert(h > l);

  /* how many negative points ? */
  nn = n * (p[k] - l)/(h-l);

  /* what is the increment? */
  if((p[k] - l) > (h-p[k]))
    inc = (p[k] - l) / nn;
  else
    inc = (h - p[k]) / (n - nn - 1);

  /* initialize u */
  for(i=0; i<len; i++)
    u[i] = p[i];

  /* print label */
  switch(k)
  {
  case 0:
    printf("\n#kappa l=%g h=%g\n", l, h);
    break;
  case 1:
    if(attrition)
      printf("\n#beta l=%g h=%g\n", l, h);
    else
      printf("\n#alpha[0] l=%g h=%g\n", l, h);
    break;
  default:
    if(attrition)
      printf("\n#alpha[%d] l=%g h=%g\n", k-2, l, h);
    else
      printf("\n#alpha[%d] l=%g h=%g\n", k-1, l, h);
  }

  /* print likelihood transect */
  for(i=0; i<n; i++)
  {
    x = u[k] = p[k] + (i-nn)*inc; 
    y = lnL(u);
    if((i+1)%3 == 0)
      putchar('\n');
    printf(" %g %g", x, y);
  }
}

/**************************************************************** 
print likelihood transects for all dimensions.  p is the base
parameter vector and has length len.  The upper bounds of kappa and
beta are given by hkappa and hbeta.  n is the number of points that
will be printed along each transect.
****************************************************************/
void all_ltransects(real *p, int len, int n, real lkappa, real hkappa,
		    real lbeta, real hbeta)
{
  extern int attrition, nagent;
  int i, k=0;
  real sum;

  assert(attrition==0 || attrition==1);
  /* kappa */
  ltransect(k, p, len, lkappa, hkappa, n);
  k++;

  if(attrition)   /* beta */
  {
    ltransect(k, p, len, lbeta, hbeta, n);
    k++;
  }

  /* loop over agents */
  while(k < attrition+nagent)
  {
    sum=0.0;
    for(i=0; i<nagent-1; i++)
    {
      if(i==k)
	continue;
      sum += p[1+attrition+i];  /* sum of alphas excluding alpha[k] */
    }
    ltransect(k++, p, len, 0.0, 1.0-sum, n);
  }
}
