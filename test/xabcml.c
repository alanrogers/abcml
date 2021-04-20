#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mytypes.h"
#include "likelihood.h"
#include "misc.h"
#include "io.h"
#include "getcic.h"
#include "amoeba.h"
#include "unirand.h"
#include "ninfomat.h"
#include "header.h"
#include "prlike.h"
#define YES(x) ((x) == 0 ? "No" : "Yes")
#define ON(x) ((x) == 0 ? "Off" : "On")
#define NQUANT 11

/*** external variables **/
real    qval[NQUANT] =
{0.001, 0.01, 0.025, 0.25, 0.5, 0.75, 0.95,
 0.975, 0.99, 0.999, 1.0};

void    usage(const char *msg);
struct results *alloc_results(void);
real    std_dev(real var);
void    pr_results(struct results *r, int whichdataset);
real    func(real * x);
void    residuals(int *yy, real * Ey, real * zscore, real * p, real ** cov);

int     npart;			/* number of parts */
int     rdim;			/* reduced number of parts */
int     nagent = 0;		/* number of agents */
real    size = 0.1;		/* size of initial simplex */
int     nparam;			/* dimension of parameter vector */
int     attrition = 1;		/* estimate attrition? */
struct agent **agent = NULL;	/* vector of pointers to agents */
struct bonedef *bones = NULL;	/* .bdf data */
int    *y;
real    scale_factor = -1.0;	/* for re-scaling attrition */

extern int constrain_lnL;

void    usage(const char *msg)
{
    fflush(stdout);
    if (msg != NULL)
	fprintf(stderr, "\n\nCmd line error: %s", msg);
    fputs("\n\nusage: abcml [options] xxx.bdf xxx.cfg [yyy.cfg ...] xxx.cnt",
	  stderr);
    fputs("\nwhere .cfg files configure agents; at least 1 is needed",
	  stderr);
    fputs("\n      .bdf file defines characteristics of bones", stderr);
    fputs("\n      .cnt file contains bone counts", stderr);
    fputs("\nFlag options may include:", stderr);
    fprintf(stderr, "\n -a    Estimate attrition (beta)? Def: %s",
	    YES(attrition));
    fprintf(stderr, "\n -i xxx Set initial parameter vector.");
    fprintf(stderr, "\n  (Here xxx=alpha0,alpha1,...kappa,beta)");
    exit(1);
}

/* create a new structure to hold results */
struct results *alloc_results(void)
{
    struct results *r;

    r = (struct results *) malloc(sizeof(struct results));
    if (r == NULL) {
	fprintf(stderr, "\nalloc_results: memory\n");
	exit(1);
    }
    r->p = (real *) mustalloc((size_t) (nparam * sizeof(real)));
    r->lnL = MISSING;
    r->ChiSq = MISSING;
    r->SCov = (real **) alloc2d(nparam, nparam, sizeof(real));
    if (r->SCov == NULL) {
	fprintf(stderr, "\nalloc_results: memory\n");
	exit(1);
    }
    r->SCov_err = 1;
    return (r);
}

/* calculate standard deviation from variance */
real    std_dev(real var)
{
    if (FINITE(var))
	return (sqrt(var));
    return (MISSING);
}

void    pr_results(struct results *r, int whichdataset)
{
    char    buff[20];
    real   *alpha, *beta, *kappa, lastalpha;
    int     i, j;

    set_mnemonics(&alpha, &beta, &kappa, &(r->p));

    /* print column labels */
    printf("\n%c", (whichdataset == 0 ? ' ' : '#'));
    printf("%7s %4s", "rowlbl", "mni");
    printf(" %12s", "kappa");
    if (attrition)
	printf(" %12s", "beta");
    for (i = 0; i < nagent; i++) {
	sprintf(buff, "alpha[%d]", i);
	printf(" %12s", buff);
    }
    if (print_lnL)
	printf(" %12s", "lnL");
    printf(" %12s", "ChiSq");

    /* print estimates */
    printf("\n%8s", "Estimate");
    printf(" %4d", r->mni);
    printf(" %12.6f", *kappa);
    if (attrition)
	printf(" %12.6f", *beta);
    lastalpha = 1.0;
    for (i = 0; i < nagent - 1; i++) {
	printf(" %12.6f", alpha[i]);
	lastalpha -= alpha[i];
    }
    printf(" %12.6f", lastalpha);
    if (print_lnL)
	printf(" %12.6f", r->lnL);
    printf(" %12.6f", r->ChiSq);

    /* print standard errors */
    if (r->SCov_err == 0) {
	printf("\n%8s", "StdErr");
	printf(" %4s", "***");	/* mni */
	printf(" %12.6f", std_dev(r->SCov[0][0]));
	if (attrition)
	    printf(" %12.6f", std_dev(r->SCov[1][1]));
	for (i = 0; i < nagent - 1; i++) {
	    j = i + 1 + attrition;
	    printf(" %12.6f", std_dev(r->SCov[j][j]));
	}
	printf(" %12s", "***");	/* lastalpha */
	if (print_lnL)
	    printf(" %12s", "***");	/* lnL */
	printf(" %12s", "***");	/* ChiSq */
    }
}
/* we want to maximize lnL, but amoeba is a minimizer.  So the */
/* objective function is -lnL                                  */
real    func(real * x)
{
    return (-lnL(x));
}
/****************************************************************
 Calculate residuals and squared Z values.
 On entry, y, Ey, and zscore should all be arrays of length npart,
 p should contain the parameter vector, and cov the covariance matrix.
 y is the data vector--the vector of skeletal part counts.  On return,
 Ey contains the expected value of y under the parameter values
 given in p and zscore is the vector of Z statistics.  The
 i'th entry of zscore equals (y[i] - Ey[i]) / sqrt(var(y[i])).
****************************************************************/
void    residuals(int *yy, real * Ey, real * zscore, real * p, real ** cov)
{
    real    resid;
    real   *alpha, *beta, *kappa;
    real    lastalpha;
    int     i, j;

    set_mnemonics(&alpha, &beta, &kappa, &p);

    /* check inequality constraints */
    lastalpha = 1.0;
    for (i = 0; i < nagent - 1; i++) {	/* alpha[i] >= 0 */
	if (alpha[i] < 0.0)
	    return;
	lastalpha -= alpha[i];
    }
    if (*kappa < 0.0 || (attrition && *beta < 0.0)
	|| lastalpha < 0.0)
	return;

    /* get vector x = yy - E{y} */
    for (i = 0; i < npart; i++) {
	Ey[i] = 0.0;
	for (j = 0; j < nagent - 1; j++)
	    Ey[i] += alpha[j] * agent[j]->m[i];
	Ey[i] += lastalpha * agent[nagent - 1]->m[i];
	Ey[i] *= *kappa;
	if (attrition)
	    Ey[i] *= exp(-*beta * bones->sensitivity[i]);
	resid = yy[i] - Ey[i];
	zscore[i] = resid / sqrt(cov[i][i]);
    }
    return;
}

int     main(int argc, char **argv)
{
    int     i, j, data, pass, nfunc, converged;
    real   *p0 = NULL;          /* vector of initial parameter values */
    real  **Fmean, sum;
    real  **p;			/* vector of parameter vectors */
    real   *alpha;		/* relevant entry in p0 */
    real   *beta;		/* relevant entry in p0 */
    real   *kappa;		/* relevant entry in p0 */
    real   *yyy;		/* function evaluations */
    real   *x;
    char    buff[20], *s;
    struct counts *counts = NULL;
    struct results *r;
    extern real lnL_ChiSq;
    extern real **lnL_cov;

    /* header */
    header("ABCML", "Debug", stdout);

    /* echo command line */
    printf("\n#Cmd line: %s", argv[0]);
    for (i = 1; i < argc; i++)
	printf(" %s", argv[i]);

    /* count .cfg arguments to determine number of agents */
    nagent = 0;
    for (i = 1; i < argc; i++) {
	if (classify(argv[i]) == CFG_ARG)
	    nagent += 1;
    }

    if (nagent < 1)
	usage("Need at least one .cfg file");

    if (nagent > 0)
	agent = (struct agent **)
	    mustalloc((size_t) (nagent * sizeof(struct agent *)));

    /* process command line arguments */
    j = 0;
    for (i = 1; i < argc; i++) {
	switch (classify(argv[i])) {
	case FLAG_ARG:
	    switch (argv[i][1]) {
	    case 'a':
		attrition = !attrition;
		break;
	    default:
		usage(argv[i]);
	    }
	    break;
	case CFG_ARG:
	    agent[j] = read_agent_cfg(argv[i]);
	    if (configure_agent(agent[j]) != 0)
		error("bad return from configure_agent");
	    j += 1;
	    break;
	case BDF_ARG:
	    if (bones != NULL)
		usage("Only one .bdf file is allowed.");
	    bones = read_bonedef(argv[i], &scale_factor);
	    break;
	case CNT_ARG:
	    if (counts != NULL)
		usage("Only one .cnt file is allowed.");
	    counts = read_counts(argv[i]);
	    break;
	default:
	    usage(argv[i]);
	}
    }

    if (bones == NULL)
	usage("Missing .bdf file");

    if (counts == NULL)
	usage("Missing .cnt file");
    if (verbose)
	printf("\n#Output is verbose.");

    check_labels(nagent, agent, bones, counts);
    nparam = (attrition ? nagent + 1 : nagent);

    /* echo parameters here */
    if (attrition)
	printf("\n#Estimating attrition parameter, beta");
    else
	printf("\n#Assuming that attrition is absent.");
    printf("\n#Number of agents           : %d", nagent);
    printf("\n#Number of parameters       : %d", nparam);
    printf("\n#Number of skeletal parts   : %d", counts->npart);

    if (nparam > counts->npart)
	error("Number of parameters > number of skeletal parts");

    /* allocations */
    r = alloc_results();
    x = (real *) mustalloc((size_t) (nparam * sizeof(real)));
    yyy = (real *) mustalloc((size_t) ((nparam + 1) * sizeof(real)));
    Fmean = (real **) alloc2d(counts->npart, counts->npart,
			      sizeof(real));
    p = (real **) alloc2d(nparam + 1, nparam, sizeof(real));
    if ( Fmean == NULL || p == NULL ) {
	fprintf(stderr, "\nmain: alloc2d\n");
	exit(1);
    }

    /* loop over datasets */
    for (data = 0; data < 1; data++) {
	free_reducmat();

	/* find minimum number of individuals */
	r->mni = 0;
	for (i = 0; i < bones->npart; i++) {
	    j = counts->y[data][i] / bones->live[i];	/* rounds down */
	    if (j * bones->live[i] < counts->y[data][i])	/* round up */
		j += 1;
	    if (j > r->mni)
		r->mni = j;
	}

	/* set externals */
	npart = counts->npart;
	y = counts->y[data];

	/* Set parameter vector */
	*kappa = 2.0 * r->mni;
	if (attrition)
	    *beta = 0.25;
	if (nagent > 1) {
	    for (i = 0; i < nagent - 1; i++)
		alpha[i] = 1.0 / nagent;
	}

	printf("\n#Params: kappa=%g", *kappa);
	if (attrition)
	    printf(" beta=%g", *beta);
	for (i = 0; i < nagent - 1; i++)
	    printf(" alpha[%d]=%g", i, alpha[i]);

	/* get mean F matrix */
	getFmean(Fmean, agent, npart, alpha, nagent);

	/* get reducmat matrix */
	printf("\n#Reducing dimension by Principal Components");
	rdim = dopc(npart, Fmean);
	printf("\n#Using %d / %d dimensions.", rdim, npart);

	r->lnL = lnL(r->p);

	constrain_lnL = 0;	/* turn off constraints */
	pr_results(r, data);

	if (print_SCov && r->SCov_err == 0)
	    prfmat(r->SCov, nparam, nparam,
		       "SCov (Sampling variances and covariances)");

	fflush(stdout);
    }				/* end dataset loop */

    putchar('\n');
    return (0);
}
