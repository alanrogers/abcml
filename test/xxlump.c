#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mytypes.h"
#include "likelihood.h"
#include "misc.h"
#include "io.h"
#include "getcic.h"
#include "amoeba.h"
#include "unirand.h"
#include "header.h"
#include "symeig.h"
#define FTOL 1.0e-6
#define YES(x) ((x) == 0 ? "No" : "Yes")
#define ON(x) ((x) == 0 ? "Off" : "On")

void usage(char *msg);
int rcompar(real *x, real *y);

int npart;                  /* number of parts */
int rdim;                   /* reduced number of parts */
int nagent=0;               /* number of agents */
int use_pc=1;               /* use principal components */
int randalpha=0;            /* random initial vector ? */
int attrition = 0;
int *y=NULL;

int nparam;                 /* dimension of parameter vector */
struct agent **agent=NULL;  /* vector of pointers to agents */
struct bonedef *bones=NULL; /* .bdf data */

void usage(char *msg)
{
  fflush(stdout);
  if(msg != NULL)
    fprintf(stderr,"\n\nCmd line error: %s", msg);
  fputs("\n\nusage: xlump [options] xxx.bdf xxx.cfg yyy.cfg [...]",
	stderr);
  fputs("\nwhere .cfg files configure agents; at least 2 are needed",
	stderr);
  fputs("\n      .bdf file defines characteristics of bones", stderr);
  fputs("\nFlag options may include:", stderr);
  fprintf(stderr,"\n -p    Principle components? Def: %s",
	  YES(use_pc));
  putc('\n', stderr);
  exit(1);
}

int rcompar(real *xx, real *yy)
{
  if(*xx < *yy)
    return(1);
  if(*xx > *yy)
    return(-1);
  return(0);
}

void main(int argc, char **argv)
{
  int i, j;
  real **Fmean;              /* mean cov matrix */
  real sum, *p0;
  real *alpha;               /* relevant portion of p0 */
  real *beta;                /* relevant entry in p0 */
  real *kappa;               /* relevant entry in p0 */
  real **c2;                 /* reduced matrix */
  real **evect1, **evect2;   /* eigenvector matrices */
  real *eigenval;            /* eigenvalues */
  real *work;
  real max, min;
  extern real **reducmat;
  extern int dim_reducmat;

  /* header */
  header("XLUMP", "Try out lump procedure", stdout);

  /* echo command line */
  printf("\n#Cmd line: %s", argv[0]);
  for(i=1; i<argc; i++)
    printf(" %s", argv[i]);

  /* count .cfg arguments to determine number of agents */
  nagent = 0;
  for(i=1; i<argc; i++)
  {
    if(classify(argv[i]) == CFG_ARG)
      nagent += 1;
  }

  agent = (struct agent **) mustalloc( nagent * sizeof(struct agent *));

  /* process command line arguments */
  nagent = 0;
  for(i=1; i<argc; i++)
  {
    switch(classify(argv[i]))
    {
    case FLAG_ARG:
      switch(argv[i][1])
      {
      case 'p':
	use_pc = !use_pc;
	break;
      default:
	usage(argv[i]);
      }
      break;
    case CFG_ARG:
      agent[nagent] = read_agent_cfg(argv[i]);
      if(configure_agent(agent[nagent]) != 0)
	error("bad return from configure_agent");
      nagent += 1;
      break;
    case BDF_ARG:
      if(bones != NULL)
	usage("Only one .bdf file is allowed.");
      bones = read_bonedef(argv[i]);
      break;
    case CNT_ARG:
      usage("No .cnt file is allowed");
      break;
    default:
      usage(argv[i]);
    }
  }

  if(bones == NULL)
    usage("Must specify a .bdf file");

  check_labels(nagent, agent, bones, NULL);

  nparam = nagent;

  /* echo parameters here (unfinished) */
  printf("\nNumber of agents    : %d", nagent);
  printf("\nInitial alpha       : %s",
	 randalpha ? "Random" : "gives equal weight to each agent");

  /* set externals */
  npart = rdim = bones->npart;

  /* allocations */
  p0 = (real *) mustalloc( nparam * sizeof(real));
  eigenval = (real *) mustalloc( npart * sizeof(real) );
  work = (real *) mustalloc( npart * sizeof(real) );
  Fmean = (real **) alloc2d(bones->npart, bones->npart, sizeof(real));
  c2 = (real **) alloc2d(bones->npart, bones->npart, sizeof(real));
  evect1 = (real **) alloc2d(bones->npart, bones->npart, sizeof(real));
  evect2 = (real **) alloc2d(bones->npart, bones->npart, sizeof(real));
  if(Fmean==NULL || p0==NULL || c2==NULL
     || evect1==NULL || evect2==NULL)
  {
    fprintf(stderr,"\nmain: alloc2d\n");
    exit(1);
  }

  /* point alpha, beta, and kappa at p0 */
  set_mnemonics(&alpha, &beta, &kappa, &p0);

  if(randalpha)
  {
    /* random initial vector */
    sum = uni();
    for(i=0; i<nagent-1; i++)
      sum += alpha[i] = uni();
    for(i=0; i<nagent-1; i++)
      alpha[i] /= sum;
  }else
  {
    /* initial vector assigns equal probabilities */
    for(i=0; i<nagent-1; i++)
      alpha[i] = 1.0/nagent;
  }

  /* get mean F matrix */
  getFmean(Fmean, agent, npart, alpha, nagent);
  prfmat(Fmean, npart, npart, "Fmean");

  /* copy Fmean into evect1 */
  for(i=0; i<npart; i++)
    for(j=0; j<npart; j++)
      evect1[i][j] = Fmean[i][j];

  /* get reducmat matrix */
  if(use_pc)
  {
    printf("\ncalling dopc");
    rdim = dopc(npart, evect1);
  }else
  {
    printf("\ncalling dolump");
    rdim = dolump(npart, evect1);
  }
  prfmat(reducmat, rdim, npart, "reducmat");

  printf("\nrdim=%d dim_reducmat=%d reducmat=%u", rdim, dim_reducmat,
	 (unsigned) reducmat);

  xaxT(reducmat, Fmean, c2, rdim, npart);
  prfmat(c2, rdim, rdim, "reduced Fmean");

  /* copy Fmean into evect1 again */
  for(i=0; i<npart; i++)
    for(j=0; j<npart; j++)
      evect1[i][j] = Fmean[i][j];

  /* copy c2 into evect2 */
  for(i=0; i<rdim; i++)
    for(j=0; j<rdim; j++)
      evect2[i][j] = c2[i][j];

  /* spectral decomposition of Fmean */
  (void) symeig(evect1, eigenval, work, npart);
  qsort(eigenval, npart, sizeof(real),
	(int (*)(const void *, const void *)) rcompar);
  printf("\n\nsorted Eigenvalues of Fmean:");
  max = min = fabs(eigenval[0]);
  for(i=0; i<npart; i++)
  {
    if(max < fabs(eigenval[i]))
      max = fabs(eigenval[i]);
    if(min > fabs(eigenval[i]))
      min = fabs(eigenval[i]);
    printf(" %g", eigenval[i]);
  }
  printf("\nmax=%g min=%g con=%g", max, min, max/min);
  
  /* spectral decomposition of c2 */
  (void) symeig(evect2, eigenval, work, rdim);
  qsort(eigenval, rdim, sizeof(real),
	(int (*)(const void *, const void *)) rcompar);
  printf("\n\nsorted Eigenvalues of c2:");
  max = min = eigenval[0];
  for(i=0; i<rdim; i++)
  {
    if(max < fabs(eigenval[i]))
      max = fabs(eigenval[i]);
    if(min > fabs(eigenval[i]))
      min = fabs(eigenval[i]);
    printf(" %g", eigenval[i]);
  }
  printf("\nmax=%g min=%g con=%g", max, min, max/min);

  putchar('\n');
  exit(0);
}
