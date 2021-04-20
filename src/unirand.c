/****************************************************************
This is a rewrite of "newrand.c", which was a rewrite of "rand.c".  The
present version is based on the uni() function 
described by 

@Book{Kahaner:NMS-89,
  author = 	 "Kahaner, David and Cleve Moler and Stephen Nash",
  title = 	 "Numerical Methods and Software",
  publisher = 	 "Prentice Hall",
  year = 	 1989,
  address = 	 "Englewood Cliffs, NJ"
}

Alan R. Rogers  4/23/93
****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mytypes.h"
#include "unirand.h"
#define MAXPERM 1000
#if !defined( FLT_MANT_DIG )
#include <float.h>
#endif

/********* globals **************************/
static real u[17];
int uni_initialized = 0;

/****************************************************************
Get a number from the operating system to use as a random number seed.
****************************************************************/
int     getseed(void)
{
  return ((int) time((time_t *) 0));
}

/*************Initialize random number generator*************/
int    initrand(int seed)
{
  if (seed == 0)
    seed = getseed();
  return(ustart(seed));
}

/****************************************************************
Fill vec with random permutation of integers from 0 to n-1
****************************************************************/
int    *randperm(int *vec, int n)
{
  int     i;
  unsigned int rnd;
  int     seq[MAXPERM], *s;

  if (n > MAXPERM)
  {
    fprintf(stderr, "\nrandperm: can't make an %d-vector. Max is %d\n",
	    n, MAXPERM);
    exit(1);
  }
  for (i = 0; i < n; i++)
    seq[i] = i;
  s = seq;
  for (i = 0; i < n; i++)
  {
    rnd = randint((int) n - i);	/* get number in interval (0,n-i-1) */
#if DEBUG
    if (rnd < 0 || rnd >= n - i)
    {
      printf("\nerr: randint(%d-%d) returned %d", n, i, rnd);
      exit(1);
    }
#endif
    vec[i] = s[rnd];
    s[rnd] = s[0];		/* s[1,2,...] contains unused numbers */
    s++;
  }
  return (vec);
}

/*****************************************************************
Input: cum[] is an n-vector containing in its i'th element the
sum of the probabilities of categories 0,1,...,i.
Output: Returns i with the appropriate probability.
****************************************************************/
int     multinomial(float *cum, int n)
{
  int     i;
  real  x;

  x = uni();
  for (i = 0; cum[i] < x && i < n; i++)
    ;
  if (i == n)
  {
    fputs("\nError in multinomial: i==n", stdout);
    exit(1);
  }
  return (i);
}

/****************************************************************
uni(): Generates random numbers distributed uniformly on [0,1].  Translated
from function of same name in:

@Book{Kahaner:NMS-89,
  author = 	 "Kahaner, David and Cleve Moler and Stephen Nash",
  title = 	 "Numerical Methods and Software",
  publisher = 	 "Prentice Hall",
  year = 	 1989,
  address = 	 "Englewood Cliffs, NJ"
}

Alan R. Rogers  4/23/93
****************************************************************/
real uni(void)
{
  static real c = 362436.0/16777216.0;   /**  2^24=16777216 **/
  static real cd = 7654321.0/16777216.0;
  static real cm = 16777213.0/16777216.0;
  static int i = 16;  /* = 17 - 1: Lag 1 = 17 */
  static int j =  4;  /* = 5 - 1 : Lag 2 = 5  */
  real result;

  if(!uni_initialized)
    (void) ustart(getseed());

  /*** Basic generator is Fibonacci ***/
      result = u[i]-u[j];
      if(result < 0.0)
	result += 1.0;
      u[i] = result;
      if(--i < 0)
	i = 16;
      if(--j < 0)
	j = 16;

  /*** Second generator is congruential ***/
      c -= cd;
      if(c<0.0)
	c += cm;

  /*** Combination generator **/
      result -= c;
      if(result < 0.0)
	result += 1.0;
      return(result);
}

/****************************************************************
ustart(seed): Initialize random number generator with given seed.
Call this before first call to uni()
****************************************************************/

int ustart(int iseed)
{
  real s,t;
  int ii,jj;
  int i1,jj1,k1,l1,m1;
  static int k = 24;

  /****************************************************************
          Set up ...
          Convert ISEED to four smallish positive integers.
  ****************************************************************/	  
  i1 = (abs(iseed) % 177)+1;
  jj1 = (abs(iseed) % 167)+1;
  k1 = (abs(iseed) % 157)+1;
  l1 = (abs(iseed) % 147)+1;

  /* Generate random bit pattern in array based on given seed. */
  for(ii=0; ii<17; ii++)
  {
    s = 0.0;
    t = 0.5;
    /****************************************
    Do for each of the bits of mantissa of word 
    Loop  over K bits, where K is defaulted to 24 but can
    be changed by defining SETK to turn on the code that follows. 
    *****************************************/
#ifdef SETK
    if(sizeof(float) == sizeof(real))  /* we're in single precision */
      k = FLT_MANT_DIG;  /* # of digits in mantissa of a float */
    else                 /* we're in double precision */
      k = DBL_MANT_DIG;  /* # of digits in mantissa of a double */
#endif    

    for(jj=1; jj<=k; jj++)
    {
      m1 = (((i1*jj1) % 179)*k1) % 179;
      i1 = jj1;
      jj1 = k1;
      k1 = m1;
      l1 = (53*l1+1) % 169;
      if( (l1*m1) % 64 > 32)
	s=s+t;
      t *= 0.5;
    }
    u[ii] = s;
  }
  /****************************************************************
  set flag so that one can tell if uni has been initialized.
  ****************************************************************/
  uni_initialized = 1;
  return(iseed);
}
  
/****************************************************************
rnorm: generate standard normal deviates.  Based on unirand.c
****************************************************************/
real rnorm(void)
{
  real v1, v2, x1, y, s;
  static real x2;
  static int gotone=0;

  if(gotone)
  {
    gotone = 0;
    return(x2);
  }
  do{
    /* form uniform variates on [-1,1] */
    v1 = 2.0*uni() - 1.0;
    v2 = 2.0*uni() - 1.0;
    s = v1*v1 + v2*v2;
  }while(s > 1.0);
  /* we now have a random point (v1,v2) in the unit circle */
  y = sqrt(-2.0*log(s)/s);
  /* now x1 and x2 are uniform random variates */
  x1 = v1*y;
  x2 = v2*y;
  gotone = 1;
  return(x1);
}
