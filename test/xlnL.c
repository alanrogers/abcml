#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mytypes.h"
#include "likelihood.h"
#include "io.h"
#include "misc.h"
#include "getcic.h"
#include "zeroin.h"

int *y;
int rdim=0, npart=0, nagent=2;
int attrition=0;
real *sensitivity=NULL;
real *p0;
struct agent **agent;
struct bonedef *bones;
struct counts *counts;

real func(real *x);
void usage(char *msg);

/* we want to maximize lnL, but amoeba is a minimizer.  So the */
/* objective function is -lnL                                  */
real func(real *x) {
    return(-lnL(x));
}

void usage(char *msg) {
    printf("\nCmd line error: %s", msg);
    printf("\nusage: xlnl residence.cfg kill.cfg bones.bdf data.cnt\n");
    exit(1);
}

int main(int argc, char **argv) {
    int i, j, nparam;
    real scale = 1.0;
  
    if(argc != 5)
        usage("wrong number of command line arguments");

    agent = (struct agent **) mustalloc( nagent * sizeof(struct agent *));

    j = 0;
    for(i=1; i<argc; i++) {
        switch(classify(argv[i])) {
        case CFG_ARG:
            agent[j] = read_agent_cfg(argv[i]);
            if(configure_agent(agent[j]) != 0)
                error("bad return from configure_agent");
            printf("\nConfigured agent %d from file %s",
                   j, argv[i]);
            pr_agent(agent[j]);
            j += 1;
            break;
        case BDF_ARG:
            if(bones != NULL)
                usage("Only one .bdf file is allowed.");
            bones = read_bonedef(argv[i], &scale);
            printf("\nConfigured bones from file %s", argv[i]);
            pr_bonedef(bones);
            break;
        case CNT_ARG:
            if(counts != NULL)
                usage("Only one .cnt file is allowed.");
            counts = read_counts(argv[i]);
            printf("\nRead faunal counts from file %s", argv[i]);
            pr_counts(counts);
            break;
        }
    }

    y = counts->y[0];
    npart = rdim = counts->npart;
    noreduc();
    check_labels(nagent, agent, bones, counts);
    nparam = (attrition ? nagent+1 : nagent);
    p0 = (real *) mustalloc( nparam * sizeof(real));

    for(;;) {
        printf("\nEnter alpha, lambda: ");
        scanf("%lf%lf", p0+1, p0+0);
        printf("\nlnL(a=%g, l=%g)=%g",
               p0[1], p0[0], lnL(p0));
    }

    putchar('\n');
    return 0;
}

