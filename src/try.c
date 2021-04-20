#include <stdio.h>
#include <math.h>
#if 1
#define FINITE(x) ((x) > -HUGE_VAL && (x) < HUGE_VAL)
#else
#define FINITE(x) ((x) == (x))
#endif
int main(void)
{
  float infinity, finite=3.0, nan;

  infinity = 1.0/0.0;
  nan = 0.0/0.0;

  if(FINITE(finite))
    printf("\n%g is finite", finite);
  else
    printf("\n%g is infinite", finite);
  if(FINITE(infinity))
    printf("\n%g is finite", infinity);
  else
    printf("\n%g is infinite", infinity);
  if(FINITE(nan))
    printf("\n%g is finite", nan);
  else
    printf("\n%g is infinite", nan);
  putchar('\n');
  return(0);
}
