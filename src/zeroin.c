#include "zeroin.h"
#define abs(x)	((x)>=0 ? (x) : -(x))
/*****************************************************************
zeroin: find the zero of a function within a given interval.

This code was translated from the fortran in the book by Forsythe, 
Malcolm, and Moler.

A zero of the function  f(x)  is computed in the interval ax,bx .

input..

ax     left endpoint of initial interval 
bx     right endpoint of initial interval
f      function subprogram which evaluates f(x) for any x in the
       interval  ax,bx
tol    desired length of the interval of uncertainty of the final
       result ( >= 0.0)

output..

zeroin abcissa approximating a zero of  f  in the interval ax,bx

It is assumed  that   f(ax)   and   f(bx)   have  opposite  signs without
a  check.  zeroin  returns a zero  x  in the given interval ax,bx  to
within a tolerance  4*macheps*abs(x) + tol, where macheps is the
relative machine precision. this function subprogram is a slightly
modified  translation  of the algol 60 procedure  zero  given in
Richard Brent, Algorithms for minimization without derivatives, 
Prentice-Hall, Inc. (1973).
*****************************************************************/
double  dsign(double x, double y);

double  zeroin(double a, double b, double (*f)(double), double tol)
{
  double  c, d, e, eps, fa, fb, fc, tol1, xm, p, q, r, s;
  /*
   * compute eps, the relative machine precision
   */
  for (eps = 1.0; 1.0 + eps > 1.0; eps /= 2.0)
    ;
  /*
   * initialization
   */
  fa = (*f) (a);
  fb = (*f) (b);
  while (1)
  {
    /*
     * begin step
     */
    c = a;
    fc = fa;
    d = b - a;
    e = d;
    do
    {
      if (abs(fc) < abs(fb))
      {				/* swap */
	a = b;
	b = c;
	c = a;
	fa = fb;
	fb = fc;
	fc = fa;
      }
      /*
       * convergence test
       */
      tol1 = 2.0 * eps * abs(b) + 0.5 * tol;
      xm = (c - b) / 2.0;
      if (abs(xm) <= tol1 || fb == 0.0)
	return (b);
      /*
       * is bisection necessary
       */
      if (abs(e) >= tol1 && abs(fa) > abs(fb))
      {
	if (a == c)
	{			/* linear interpolation */
	  s = fb / fa;
	  p = 2.0 * xm * s;
	  q = 1.0 - s;
	} else
	{			/* inverse quadratic interpolation */
	  q = fa / fc;
	  r = fb / fc;
	  s = fb / fa;
	  p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
	  q = (q - 1.0) * (r - 1.0) * (s - 1.0);
	}
	/*
	 * adjust signs
	 */
	if (p > 0.0)
	  q = -q;
	p = abs(p);
	/*
	 * is interpolation acceptable
	 */
	if (2.0 * p < 3.0 * xm * q - abs(tol1 * q) && p < abs(0.5 * e * q))
	{
	  e = d;
	  d = p / q;
	  goto NEWB;
	}
      }
      /*
       * bisection
       */
      d = xm;
      e = d;
      /*
       * complete step
       */
  NEWB:a = b;
      fa = fb;
      if (abs(d) > tol1)
	b = b + d;
      else
	b = b + dsign(tol1, xm);
      fb = (*f) (b);
    } while (fb * (fc / abs(fc)) <= 0.0);
  }
  return (-99.99);  /* NOTREACHED */
}
/*
 * I really don't know what the fortran function dsign() does. This is a
 * guess. This function returns a value whose magnitude is the same as x,
 * and whose sign is the same as that of y.
 */
double  dsign(double x, double y)
{
  x = abs(x);
  if (y >= 0.0)
    return (x);
  return (-x);
}
