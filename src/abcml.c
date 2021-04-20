/****************************************************************
    abcml: Analysis of Bone Counts by Maximum Likelihood
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
int     use_pc = 1;		/* use principal components */
int     print_lnL = 0;		/* print log likelihood? */
int     print_SCov = 0;
int     print_residuals = 1;
int     do_chisq_quantiles = 0;	/* do quantiles? */
int     force_chisq_quantiles = 0;	/* force quantiles regardless of ndataset? */
int     quantile_threshold = 50;	/* do quantiles if ndaset>quantile_threshold */
int     likelihood_transects = 0;	/* print likelihood transects? */
int     randinit = 0;		/* random initial vector ? */
real    size = 0.1;		/* size of initial simplex */
int     nparam;			/* dimension of parameter vector */
int     attrition = 1;		/* estimate attrition? */
int     verbose = 0;
int     which = -1;		/* which column should I use ? */
struct agent **agent = NULL;	/* vector of pointers to agents */
struct bonedef *bones = NULL;	/* .bdf data */
int    *y;
real    scale_factor = -1.0;	/* for re-scaling attrition */
real    lkappa = 1.0;
real    hkappa = 200.0;
real    lbeta = 0.0;
real    hbeta = 2.0;

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
    fprintf(stderr, "\n -C    Print matrix of sampling covariances");
    fprintf(stderr, " Def: %s", YES(print_SCov));
    fprintf(stderr, "\n -D x  Set sensitivity to x/density. Def: automatic");
    fprintf(stderr, "\n -e x  Set size of initial simplex. Def: %g",
	    size);
    fprintf(stderr, "\n -i xxx Set initial parameter vector.");
    fprintf(stderr, "\n  (Here xxx=alpha0,alpha1,...kappa,beta)");
    fprintf(stderr, "\n -L    Print lnL? Def: %s", YES(print_lnL));
    fprintf(stderr, "\n -p    Principle components? Def: %s",
	    YES(use_pc));
    fprintf(stderr, "\n -q    Quantiles of ChiSq?");
    if (force_chisq_quantiles > 0)
	fprintf(stderr, " Def: %s", YES(1));
    else if (force_chisq_quantiles < 0)
	fprintf(stderr, " Def: %s", YES(0));
    else
	fprintf(stderr, " Def: only if ndatasets>%d", quantile_threshold);
    fprintf(stderr, "\n -R    Print residuals? Def: %s",
	    YES(print_residuals));
    fprintf(stderr, "\n -r    Random initial parameters?  Def: %s",
	    YES(randinit));
    fprintf(stderr, "\n -t lk,hk,lb,hb  Print likelihood transects?");
    fprintf(stderr, "  Def: %s", YES(likelihood_transects));
    fprintf(stderr, "\n    lk,hk give low and hi bounds on kappa. Def: %g %g",
	    lkappa, hkappa);
    fprintf(stderr, "\n    lb,hb give low and hi bounds on beta. Def: %g %g",
	    lbeta, hbeta);
    fprintf(stderr, "\n -w i  Use column i only of data in .cnt file.");
    fprintf(stderr, "\n -v    Toggle verbose mode.  Def: %s\n",
	    ON(verbose));
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
    r->p = (real *) mustalloc(((size_t) nparam) * sizeof(real));
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
    real   *Ey;			/* expected value of y */
    real   *zscore;		/* vector of Z values */
    real   *ChiSqVec;		/* ChiSqVec[i] = ChiSq for dataset i */
    real   *x, tbone, sbone;
    char    buff[20], *s, *paramString = NULL;
    struct counts *counts = NULL;
    struct results *r;
    extern real lnL_ChiSq;
    extern real **lnL_cov;

    /* header */
    header("ABCML", "Analysis of Bone Counts by Maximum Likelihood", stdout);

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
	    mustalloc(((size_t) nagent) * sizeof(struct agent *));

    /* process command line arguments */
    j = 0;
    for (i = 1; i < argc; i++) {
	switch (classify(argv[i])) {
	case FLAG_ARG:
	    switch (argv[i][1]) {
	    case 'a':
		attrition = !attrition;
		break;
	    case 'C':
		print_SCov = !print_SCov;
		break;
	    case 'D':
		if (++i >= argc)
		    usage("Missing arg for -D");
		if (bones != NULL)
		    error("-D must come before .bdf on command line");
		scale_factor = strtod(argv[i], NULL);
		break;
	    case 'e':
		if (++i >= argc)
		    usage("Missing arg for -e");
		size = strtod(argv[i], NULL);
		printf("\n#Reset size=%g", size);
		break;
	    case 'i':
		if (++i >= argc)
		    usage("Missing arg for -i");
		paramString = argv[i];
		break;
	    case 'L':
		print_lnL = !print_lnL;
		break;
	    case 'p':
		use_pc = !use_pc;
		break;
	    case 'q':
		switch (force_chisq_quantiles) {
		case -1:
		case -0:
		    force_chisq_quantiles = 1;
		    break;
		case 1:
		    force_chisq_quantiles = -1;
		    break;
		default:
		    error("bad value in -q switch");
		}
		break;
	    case 'R':
		print_residuals = !print_residuals;
		break;
	    case 'r':
		randinit = !randinit;
		break;
	    case 't':
		likelihood_transects = !likelihood_transects;
		if (++i >= argc)
		    usage("Missing arg for -t");
		s = strtok(argv[i], ",");
		if (s == NULL)
		    usage("Missing arg for -t");
		lkappa = strtod(s, NULL);
		s = strtok(NULL, ",");
		if (s == NULL)
		    usage("Missing arg for -t");
		hkappa = strtod(s, NULL);
		s = strtok(NULL, ",");
		if (s == NULL)
		    usage("Bad arg for -t");
		lbeta = strtod(s, NULL);
		s = strtok(NULL, ",");
		if (s == NULL)
		    usage("Bad arg for -t");
		hbeta = strtod(s, NULL);
		break;
	    case 'v':
		verbose = !verbose;
		break;
	    case 'w':
		if (++i >= argc)
		    usage("Missing arg for -w");
		which = strtol(argv[i], NULL, 10) - 1;
		if (which < 0)
		    error("-w must specify an integer > 0");
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

  /*****************************************************
    w/ attrition                   w/o attrition
    p[0] = kappa                  p[0] = kappa
    p[1] = beta                    p[1] = alpha[0]
    p[2] = alpha[0]                p[2] = alpha[1]
    ....                           ....
    p[nagent] = alpha[nagent-2]    p[nagent-1] = alpha[nagent-2]
    nparam = nagent+1              nparam = nagent
  ******************************************************/
    nparam = (attrition ? nagent + 1 : nagent);



    /*
     * If paramString != NULL, then parameters have been set on the
     * command line.  Parse them here.
     */
    if (paramString != NULL) {
	p0 = (real *) mustalloc(((size_t) nparam) * sizeof(real));
	/* point alpha, beta, and kappa at p0 */
	set_mnemonics(&alpha, &beta, &kappa, &p0);
	for (j = 0; j < nagent-1; j++) {
	    s = strtok(paramString, ",");
	    if (s == NULL)
		usage("Bad arg for -i");
	    alpha[j] = strtod(s, NULL);
	}
	s = strtok(NULL, ",");
	if (s == NULL)
	    usage("Bad arg for -i");
	s = strtok(NULL, ",");
	*kappa = strtod(s, NULL);  /* will be overwritten */
	if (s == NULL)
	    usage("Bad arg for -i");
	*kappa = strtod(s, NULL);  /* overwritten here */
	if (attrition) {
	    s = strtok(NULL, ",");
	    if (s == NULL)
		usage("Bad arg for -i");
	    *beta = strtod(s, NULL);
	}
    }

    /* decide whether to do quantiles */
    if (force_chisq_quantiles > 0)
	do_chisq_quantiles = 1;
    else if (force_chisq_quantiles < 0)
	do_chisq_quantiles = 0;
    else {
	if (which >= 0)
	    do_chisq_quantiles = 0;
	else if (counts->ndataset > quantile_threshold)
	    do_chisq_quantiles = 1;
	else
	    do_chisq_quantiles = 0;
    }

    /* echo parameters here */
    if (attrition)
	printf("\n#Estimating attrition parameter, beta");
    else
	printf("\n#Assuming that attrition is absent.");
    printf("\n#Output is                  : %s",
	   (verbose ? "verbose" : "not verbose"));
    printf("\n#Number of agents           : %d", nagent);
    printf("\n#Number of parameters       : %d", nparam);
    printf("\n#Number of skeletal parts   : %d", counts->npart);
    printf("\n#Sensitivity to attrition   : %0.20g / density", scale_factor);
    printf("\n#Parameter initializatiton  : %s", (randinit ? "random" :
						  "fixed"));
    if(paramString != NULL) {
	printf("\n#Initial parameter vector   :");
	for(i=0; i < nparam; i++)
	    printf(" %g", p0[i]);
    }

    if (nagent > 1)
	printf("\n#F mean gives %s weight to each agent",
	       randinit ? "random" : "equal");
    if (nparam > counts->npart)
	error("Number of parameters > number of skeletal parts");
    if (verbose) {
	for (i = 0; i < nagent; i++) {
	    printf("\n#Data for agent %d:", i);
	    pr_agent(agent[i]);
	}

	printf("\n#\n#Bone definition data:");
	pr_bonedef(bones);

	printf("\n#Archeological data:");
	pr_counts(counts);
    }
    /* initialize random number generator */
    initrand(0);

    /* allocations */
    r = alloc_results();
    x = (real *) mustalloc(((size_t) nparam) * sizeof(real));
    yyy = (real *) mustalloc(((size_t) (nparam + 1)) * sizeof(real));
    Ey = (real *) mustalloc((size_t) (counts->npart * sizeof(real)));
    zscore = (real *) mustalloc((size_t) (counts->npart * sizeof(real)));
    ChiSqVec = (real *) mustalloc((size_t) (counts->ndataset * sizeof(real)));
    Fmean = (real **) alloc2d(counts->npart, counts->npart,
			      sizeof(real));
    p = (real **) alloc2d(nparam + 1, nparam, sizeof(real));
    if ( Fmean == NULL || p == NULL ) {
	fprintf(stderr, "\nmain: alloc2d\n");
	exit(1);
    }

    set_mnemonics(&alpha, &beta, &kappa, &r->p);

    /* loop over datasets */
    for (data = 0; data < counts->ndataset; data++) {
	if (which >= 0 && which != data)
	    continue;
	printf("\n#\n### Dataset %d", data + 1);
	fflush(stdout);
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

	if(paramString != NULL) { /* use param values from command line */
	    memcpy(r->p, p0, nparam * sizeof(real) );
	}else if (randinit) {	  /* initial parameter vector is random */
	    *kappa = uni() * 2.0 * r->mni;	/* kappa */
	    if (attrition)
		*beta = -log(uni());	/* beta   */
	    if (nagent > 1) {
		sum = uni();
		for (i = 0; i < nagent - 1; i++)
		    sum += alpha[i] = uni();
		for (i = 0; i < nagent - 1; i++)
		    alpha[i] /= sum;
	    }
	} else {		/* initial parameter vector is fixed */
	    *kappa = 2.0 * r->mni;
	    if (attrition)
		*beta = 0.25;
	    if (nagent > 1) {
		for (i = 0; i < nagent - 1; i++)
		    alpha[i] = 1.0 / nagent;
	    }
	}
	printf("\n#Initial params: kappa=%g", *kappa);
	if (attrition)
	    printf(" beta=%g", *beta);
	for (i = 0; i < nagent - 1; i++)
	    printf(" alpha[%d]=%g", i, alpha[i]);

	/*
	 * Fmean is a weighted mean of per-agent F matrices.
	 * The weights are the entries of alpha.
	 */
	getFmean(Fmean, agent, npart, alpha, nagent);

	/* get reducmat matrix */
	if (use_pc) {
	    /*
	     * Principal components finds the eigenvectors of
	     * Fmean and discards those whose eigenvalues are small.
	     */
	    printf("\n#Reducing dimension by Principal Components");
	    rdim = dopc(npart, Fmean);
	} else {
	    printf("\n#Reducing dimension by heuristic lumping");
	    rdim = dolump(npart, Fmean);
	}

	printf("\n#Using %d / %d dimensions.",
	       rdim, npart);
	fflush(stdout);
#if 1
	if (nparam > rdim) {
	    printf("#\nParameter count > dimension...Skipping this data set");
	    continue;
	}
#endif

	converged = 1;
	for (pass = 1; pass <= 2; pass++) {
	    if (verbose) {
		printf("\n#\n#PASS %d...", pass);
		fflush(stdout);
	    }
	    for (i = 0; i <= nparam; i++) {
		for (j = 0; j < nparam; j++) {
		    p[i][j] = r->p[j] * (1.0 + (i==j ? size : -size));
		    x[j] = p[i][j];
		}
		yyy[i] = func(x);
	    }

#if 0
	    /* print initial simplex */
	    printf("\n#Initial simplex:");
	    printf("\n#\n#%3s", "i");
	    for (j = 0; j < nparam; j++) {
		sprintf(buff, "p[i][%d]", j);
		printf(" %12s", buff);
	    }
	    printf("%14s\n#", "function");
	    for (i = 0; i < (nparam + 1); i++) {
		printf("\n#%3d ", i);
		for (j = 0; j < nparam; j++)
		    printf("%12.6f ", p[i][j]);
		printf("%12.6f", yyy[i]);
	    }
#endif

	    constrain_lnL = 1;	/* do constrained optimization */

	    nfunc = amoeba(p, yyy, nparam, FTOL, func);
	    if (nfunc < 0) {
		if (verbose)
		    printf("\n#NO CONVERGENCE in amoeba pass %d", pass);
		nfunc = -nfunc;
		converged = 0;
	    } else {
		if (verbose)
		    printf("\n#amoeba pass %d converged", pass);
		converged = 1;
	    }
	    if (verbose)
		printf("\n#%d function evaluations", nfunc);
	    
	    /* make r->p equal max over simplex */
	    j = 0;
	    for (i = 1; i < (nparam + 1); i++)
		if (yyy[i] > yyy[j])
		    j = i;	/* j is index that maximizes lnL */
	    r->lnL = -yyy[j];
	    for (i = 0; i < nparam; i++)
		r->p[i] = p[j][i];
	}	/* end pass loop */

	if (verbose) {
	    printf("\n#Final simplex:");
	    printf("\n#\n#%3s", "i");
	    for (j = 0; j < nparam; j++) {
		sprintf(buff, "p[i][%d]", j);
		printf(" %12s", buff);
	    }
	    printf("%14s\n#", "function");
	    for (i = 0; i < (nparam + 1); i++) {
		printf("\n#%3d ", i);
		for (j = 0; j < nparam; j++)
		    printf("%12.6f ", p[i][j]);
		printf("%12.6f", yyy[i]);
	    }
	}

	if (converged) {
	    constrain_lnL = 0;	/* turn off constraints */
	    r->SCov_err = ninfomat(r->SCov, r->p, nparam, lnL, attrition);
	    if (r->SCov_err)
		printf("\n#Warning: could not calculate SCov");

	    (void) lnL(r->p);
	    r->ChiSq = ChiSqVec[data] = lnL_ChiSq;
	    (void) residuals(counts->y[data], Ey, zscore, r->p, lnL_cov);
	    pr_results(r, data);

	    if (attrition) {
		/* calculate expected survival fraction */
		tbone = 0;
		sbone = 0.0;
		for (i = 0; i < bones->npart; i++) {
		    tbone += bones->live[i];
		    sbone += bones->live[i]
			* exp(-*beta * bones->sensitivity[i]);
		}
		printf(
		"\n# %g%% of a whole animal would survive attrition",
		       100.0 * sbone / tbone);
	    }
	    if (print_residuals) {
		printf("\n\n# Residuals:");
		printf("\n#%-17s %5s %12s %12s %12s",
		       "label", "y", "mu", "y-mu",
		       "Z");
		for (i = 0; i < counts->npart; i++)
		    printf("\n#%-17s %5d %12.6f %12.6f %12.6f",
			   counts->label[i], counts->y[data][i], Ey[i],
			   counts->y[data][i] - Ey[i], zscore[i]);

	    }
	    if (print_SCov && r->SCov_err == 0)
		prfmat(r->SCov, nparam, nparam,
		       "SCov (Sampling variances and covariances)");
	    if (likelihood_transects)
		all_ltransects(r->p, nparam, 20, lkappa, hkappa, lbeta, hbeta);
	} else
	    printf("\nNo convergence");
	fflush(stdout);
    }				/* end dataset loop */

    if (do_chisq_quantiles) {
	qsort(ChiSqVec, (unsigned) (counts->ndataset), sizeof(real),
	      (int (*)(const void *, const void *)) compar);

	printf("\n# QUANTILES of ChiSq statistic:");
	for (i = 0; i < NQUANT; i++)
	    printf("\n# %6.4f: %11.6f",
		 qval[i], quantile(qval[i], ChiSqVec, counts->ndataset));
    }
    putchar('\n');
    return (0);
}
