/****************************************************************
    wgtcfg: Assign probabilities to configurations based on a model
    that assumes that foragers prefer parts that are light but rich
    in meat.
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
#include <string.h>
#include <time.h>
#include <assert.h>
#include "mytypes.h"
#include "misc.h"
#include "io.h"
#include "unirand.h"
#include "header.h"
#define YES(x) ((x) == 0 ? "No" : "Yes")
#define ON(x) ((x) == 0 ? "Off" : "On")

/* function prototypes */
void     usage(const char *msg);
int      real_lookup(real value, real *vec, int length);
real     kwatts(real l, real w, real v, real g, real eta);
double   lnfact(int n);
double   poisson(int x, real mean);
double   gammln(double xx);


/* external variables */
int npart;                  /* number of parts */
int nagent=0;               /* number of agents */
int K=100;                  /* number of animals in each assemblage */
int ndataset=1;             /* number of simulated datasets */
int utility_only=1;         /* weight by utility only */
real *alpha;                /* vector of alpha values */
real beta=1.0;              /* intensity of attrition */
int multiply=0;             /* multiply weights by pr ? */
real inflate_size=1.0;      /* scale the animal up or down */
real tol = 1e-4;            /* governs number of terms in sum */
real lambda = 0.25;         /* mean number of extra carriers */
real max_load = 20.0;       /* max: 20 kg per carrier */
real hours = 2.0;    /* hours that load must be carried */

#define KEEPSIZE 100
/* log of factorial */
double lnfact(int n)
{
  static int first_call = 1;
  static double a[KEEPSIZE];
  int j;

  if(first_call)
  {
    first_call = 0;
    for(j=0; j<KEEPSIZE; j++)
      a[j] = 0.0;
  }
  assert(n >= 0);

  if(n <= 1)
    return(0.0);
  if (n >= KEEPSIZE)
    return(gammln(n+1.0));
  if(a[n] == 0.0)
    a[n] = gammln(n + 1.0);
  return a[n];
}

/* poisson probability function */
double poisson(int x, double mean)
{
  double lnprob;

  assert(x >= 0);
  assert(mean >= 0.0);
  switch(x)   /* for small x, do it the simple way */
  {
  case 0:
    return(exp(-mean));
  case 1:
    return(mean*exp(-mean));
  case 2:
    return(mean*mean*exp(-mean)/2.0);
  case 3:
    return(mean*mean*mean*exp(-mean)/6.0);
  }

  /* for larger x, work with logs to avoid large numbers */
  lnprob = x*log(mean) - mean - lnfact(x);
  return(exp(lnprob));
}

/* log gamma, this came from Numerical Recipes */
double   gammln(double xx)
{
  double  x, tmp, ser;
  static double cof[6] = {76.18009173, -86.50532033, 24.01409822,
  -1.231739516, 0.120858003e-2, -0.536382e-5};
  int     j;

  x = xx - 1.0;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.0;
  for (j = 0; j <= 5; j++)
  {
    x += 1.0;
    ser += cof[j] / x;
  }
  return -tmp + log(2.50662827465 * ser);
}

void usage(const char *msg)
{
  fflush(stdout);
  if(msg != NULL)
    fprintf(stderr,"\n\nError on command line: %s", msg);
  fputs("\n\nusage: wgtcfg xxx.wgt xxx.cfg [options]",
	stderr);
  fputs("\nwhere:",stderr);
  fputs("\n      .wgt file defines weights of skeletal parts", stderr);
  fputs("\n      .cfg file configures agent",
	stderr);
  fputs("\n options may include:", stderr);
  fprintf(stderr,"\n -c x Mean number of extra carriers. Def: %g",
	  lambda);
  fprintf(stderr,"\n -h x Hours that load must be carried. Def: %g",
	  hours);
  fprintf(stderr,"\n -i x Inflate size by factor x. Def: %g",
	  inflate_size);
  fprintf(stderr,"\n -l x Maximum load per carrier. Def: %g kg",
	  max_load);
  fprintf(stderr,"\n -m   Multiply weights by existing pr values?");
  fprintf(stderr,"  Def: %s", YES(multiply));
  fprintf(stderr,"\n -u   Weight by utility only. Def: %s",
	  YES(utility_only));
  putc('\n', stderr);
  exit(1);
}

/****************************************************************
Calculate metabolic cost of carrying a load.  

The formula is taken from Pandolf, Givoni, and
Goldman. 1977. Predicting energy expenditure with loads while standing
or walking very slowly.  Journal of Applied Physiology. 43. 577-581.

Units: weights are in grams, velocity is in meters per second.

On input:
  l = weight of load
  w = weight of person carrying the load
  v = velocity
  g = grade, measured as percent slope
 eta = terrain coefficient (how measured?)
 
On output:
  function returns the metabolic cost in kilowatts (Kjoules per second).
****************************************************************/
real kwatts(real l, real w, real v, real g, real eta)
{
  real m, r, gwt;

  r = l/w;
  gwt = w+l;
  m = 1.5*w + 2.0*gwt*r*r + eta * gwt*(1.5*v*v + 0.35*v*g);
  m *= 0.001; /* convert from watts to kilowatts */

  return(m);
}

/****************************************************************
Adjust the probabilities of the configurations in a .cfg file
using skeletal part weights specified in a .cfg file.
****************************************************************/
int main(int argc, char **argv)
{
  int i, j, carriers;
  real value, load, cost, net;
  real w, v, g, eta, seconds;
  real sum, prob;
  struct agent *agent=NULL;  /* pointer to agent */
  struct wgt *weight=NULL;   /* .wgt data */
  real joules_per_gram = 6600.0;
  /****************************************************************
  In 100 grams of roasted deer meat, there are 158 kcal of food
  energy or 660 kjoule.  Source: Anderson, B.A. 1989. Composition of
  Foods: Lamb, Veal, and Game Products. Agricultural Handbook Number
  8-17. US Department of Agriculture.
  ****************************************************************/

  /* header */
  header("WGTCFG", "Weight Configurations", stdout);

  /* echo command line */
  printf("\n#Cmd line: %s", argv[0]);
  for(i=1; i<argc; i++)
    printf(" %s", argv[i]);

  /* process command line arguments */
  for(i=1; i<argc; i++)
  {
    switch(classify(argv[i]))
    {
    case FLAG_ARG:
      switch(argv[i][1])
      {
      case 'c':
	if(++i >= argc)
	  usage("missing argument to -c");
	lambda = strtod(argv[i], NULL);
	break;
      case 'h':
	if(++i >= argc)
	  usage("missing argument to -h");
	hours = strtod(argv[i], NULL);
	break;
      case 'i':
	if(++i >= argc)
	  usage("missing argument to -i");
	inflate_size = strtod(argv[i], NULL);
	break;
      case 'l':
	if(++i >= argc)
	  usage("missing argument to -l");
	max_load = strtod(argv[i], NULL);
	break;
      case 'm' :
	multiply = !multiply;
	break;
      case 'u':
	utility_only = !utility_only;
	break;
      default:
	usage(argv[i]);
      }
      break;
    case CFG_ARG:
      if(agent != NULL)
	usage("Only one .cfg file is allowed.");
      agent = read_agent_cfg(argv[i]);
      printf("\n#Configured agent from file %s", argv[i]);
      break;
    case WGT_ARG:
      if(weight != NULL)
	usage("Only one .wgt file is allowed.");
      weight = read_wgt(argv[i]);
      printf("\n#Read weights from file %s", argv[i]);
      break;
    default:
      usage(argv[i]);
    }
  }

  if(agent == NULL)
    usage("Missing .cfg file");
  if(weight == NULL)
    usage("Missing .wgt file");

  /* make sure data are consistent */
  if(agent->npart != weight->npart)
    error("nparts is inconsistent in input files");
  if(lbls_differ(weight->label, agent->label, agent->npart))
    error("Labels disagree");

  /* initialize constants */
  w   = 55;     /* weight of forager in kilograms */
  v   = 1.0;   /* 1 m/s = 2.16 m/h */
  g   = 0.0;    /* percent slope */
  eta = 1.1;    /* walking along dirt trail */
  seconds = hours * 60.0 * 60.0; /* time required to transport load    */

  /* report constants */
  if(!utility_only)
  {
    printf("\n# %25s = %g", "weight of forager in kg", w);
    printf("\n# %25s = %g", "velocity (m/s)", v);
    printf("\n# %25s = %g", "percent slope", g);
    printf("\n# %25s = %g", "terrain coefficient (eta)", eta);
    printf("\n# %25s = %g", "distance (hours)", hours);
    printf("\n# %25s = 1 + %g", "mean carriers", lambda);
    printf("\n# %25s = %g", "max load", max_load);
    printf("\n# %25s = %g", "joules/gram", joules_per_gram);
  }
  printf("\n# %25s = %g", "size inflation factor", inflate_size);
  printf("\n# utility_only = %s", YES(utility_only));

  /* calculate net value */
  for(i=0; i<agent->nconfig; i++)
  {
    value = 0.0;
    load = 0.0;
    for(j=0; j<agent->npart; j++)
    {
      value += agent->cfg[j][i] * weight->value[j];
      load += agent->cfg[j][i] * weight->weight[j];
    }
    load *= 0.001;  /* load in kilograms */
    value *= 0.001; /* value in kilograms */
    value *= inflate_size;  /* adjust size */
    load *= inflate_size;   /* adjust size */
    value *= joules_per_gram; /* kgrams -> kjoules */
    if(utility_only)
      net = value;
    else
    {
      net = sum = 0.0;
      for(carriers=1; ; carriers++)
      {
	/* the number of carriers is a poisson r.v. */
	prob = poisson(carriers-1, lambda);
	if(prob < sum*tol)       /* stop when prob/sum < tol */
	  break;
	sum += prob;
	/* calculate cost per carrier */
	cost = kwatts(load/carriers, w, v, g, eta) * seconds; /* kjoules */
	cost *= carriers;        /* total cost */
	if(load <= max_load*carriers /* can the load be carried?*/
	   && value > cost)          /* is it worth carrying? */
	{
	  net += prob*(value - cost);
	}
      }
      net /= sum;   /* make probabilities sum to 1 */
    }
    if(multiply)
      agent->pr[i] *= net;
    else
      agent->pr[i] = net;
  }
  pr_agent(agent);
  putchar('\n');
  return(0);
}
