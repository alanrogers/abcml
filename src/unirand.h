/** unirand.h: Prototypes for functions in unirand.c **/

/*Must somewhere define a type called "real" by e.g. "typedef double real"*/
/****************************************************************
randint(n) returns a random integer between 0 and n-1.  It is now
implemented as a macro, which should improve speed.  The float-int
conversion rounds down, thus giving an int uniformly distributed on
0,1,...,(n-1).  We will never get n because uni() is uniformly
distributed on [0,1), not [0,1].
****************************************************************/
#define randint(n)  ( (int) (uni() * (n)))
/*********Defined in unirand.c************/
int     getseed(void);
int     initrand(int seed);
int    *randperm(int *vec, int n);
int     multinomial(float *cum, int n);
real    uni(void);
int     ustart(int iseed);
real    rnorm(void);
