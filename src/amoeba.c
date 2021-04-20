/****************************************************************
    amoeba: minimization by downhill simplex method.

    This code was adapted from that in Numerical Recipes in C, 2nd
    edition (by WH Press, SA Teukolksy, and WT Vetterling., Cambridge
    University Press).
****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mytypes.h"
#include "amoeba.h"

#define NMAX 5000
#define Alpha 1.0
#define Beta 0.5
#define GAMMA 2.0
void get_psum(real *psum, real **p, int ndim, int mpts);
real amotry(real **p, real *y, real *psum, int ndim,
	     real (*funk)(real *psum), int ihi, int *nfunk,
	     real fac);

void get_psum(real *psum, real **p, int ndim, int mpts)
{
  int i, j;

  for (j=0; j<ndim; j++)
  {
    psum[j] = 0.0;
    for (i=0; i<mpts; i++)
      psum[j] += p[i][j];
  }  
}

/* returns number of function evaluations on success,  */
/* the negative of that number on failure 1 on failure */
int amoeba(real **p,real *y, int ndim, real ftol,
	    real (*funk)(real *psum))
{
  int i, j, ilo, ihi, inhi, mpts=ndim+1, nfunk=0;
  real ytry, ysave, rtol, *psum;

  psum= (real *) malloc( ndim * sizeof(real) );
  if(psum==NULL)
  {
    fprintf(stderr,"\namoeba: no memory\n");
    exit(1);
  }
  get_psum(psum, p, ndim, mpts);
  for (;;)
  {
    ilo=0;
    ihi = y[0]>y[1] ? (inhi=1, 0) : (inhi=0, 1);
    for (i=0; i<mpts; i++)
    {
      if (y[i] < y[ilo])
	ilo=i;
      if (y[i] > y[ihi])
      {
	inhi = ihi;
	ihi = i;
      }else if(y[i] > y[inhi])
	if(i != ihi)
	  inhi=i;
    }
    rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
    if (rtol < ftol)
      break;
    if (nfunk >= NMAX)
      return(-nfunk);
    ytry=amotry(p, y, psum, ndim, funk, ihi, &nfunk, -Alpha);
    if (ytry <= y[ilo])
      ytry=amotry(p, y, psum, ndim, funk, ihi, &nfunk, GAMMA);
    else if (ytry >= y[inhi])
    {
      ysave=y[ihi];
      ytry=amotry(p, y, psum, ndim, funk, ihi, &nfunk, Beta);
      if (ytry >= ysave)
      {
	for (i=0; i<mpts; i++)
	{
	  if (i != ilo)
	  {
	    for (j=0; j<ndim; j++)
	    {
	      psum[j]=0.5*(p[i][j]+p[ilo][j]);
	      p[i][j]=psum[j];
	    }
	    y[i]=(*funk)(psum);
	  }
	}
	nfunk += ndim;
	get_psum(psum,  p, ndim, mpts);
      }
    }
  }
  free(psum);
  return(nfunk);
}

real amotry(real **p, real *y, real *psum, int ndim,
	     real (*funk)(real *psum), int ihi, int *nfunk,
	     real fac)
{
  int j;
  real fac1, fac2, ytry, *ptry;

  ptry = (real *) malloc(ndim * sizeof(real));
  if(ptry==NULL)
  {
    fprintf(stderr,"\namotry: no memory\n");
    exit(1);
  }
  fac1=(1.0-fac)/ndim;
  fac2=fac1-fac;
  for (j=0; j<ndim; j++)
    ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;
  ytry=(*funk)(ptry);
  ++(*nfunk);
  if (ytry < y[ihi])
  {
    y[ihi]=ytry;
    for (j=0; j<ndim; j++)
    {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j];
    }
  }
  free(ptry);
  return(ytry);
}

#undef Alpha
#undef Beta
#undef GAMMA
#undef NMAX
