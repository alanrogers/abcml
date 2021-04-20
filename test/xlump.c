#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "mytypes.h"
#include "likelihood.h"
#include "misc.h"
#include "io.h"
#define NAGENT 1
#define NPART 4

int nagent = NAGENT;
int npart = NPART, npart_alloc = NPART;
int *y;
real **lump = NULL;
real ***F;
real **m;

int main(int argc, char **argv) {
    int i, j, k;

    /* allocate F */
    F = (real ***) mustalloc(nagent * sizeof(real));
    for(i = 0; i < nagent; i++) {
        F[i] = (real **) alloc2d(npart, npart, sizeof(real));
        if(F[i] == NULL) {
            fprintf(stderr, "\nmain: alloc2d\n");
            exit(1);
        }
    }

    // Next part is broken: FF isn't defined.
    // This is apparently a program I began writing in 2000 but
    // never finished. I'm leaving it as is for the moment, but
    // this file should probably be deleted.

    /* set F = FF */
    for(i = 0; i < nagent; i++)
        for(j = 0; j < npart; j++)
            for(k = 0; k < npart; k++)
                F[i][j][k] = FF[i][j][k];

    /* allocate m */
    m = (real **) alloc2d(nagent, npart, sizeof(real));
    y = (int *) mustalloc(npart * sizeof(int));

    /* set y = yy */
    for(i = 0; i < npart; i++)
        y[i] = yy[i];

    /* set m = mm */
    for(i = 0; i < nagent; i++)
        for(j = 0; j < npart; j++)
            m[i] = mm[i];

    prfmat(F[0], npart, npart, "F");
    prfvec(m[0], npart, "m");
    privec(y, npart, "y");
    printf("\n\nlumping indices 1 and 2");
    dolump(npart, F[0]);
    prfmat(F[0], npart, npart, "F");
    prfvec(m[0], npart, "m");
    privec(y, npart, "y");
    prfmat(lump, NPART, npart - 1, "lump");

    putchar('\n');
    return 0;
}
