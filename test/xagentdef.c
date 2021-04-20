#include <stdio.h>
#include <stdlib.h>
#include "mytypes.h"
#include "io.h"
#include "misc.h"

struct bonedef *alloc_bonedef(void);
struct agent *alloc_agent(void);
struct bonedef *read_bonedef(char *fname);
void pr_bonedef(struct bonedef *b);
char *copystring(char *a, int size);
void pr_agent(struct agent *a);
double sfun(double x);
struct agent *read_agent_cfg(char *fname);;

/* externals to keep bone.c happy */
int attrition=1;

void main(int argc, char **argv)
{
  struct agent *a;

  if(argc != 2)
  {
    fprintf(stderr,"usage: xagentdef inputfile.cfg\n");
    exit(1);
  }
  a = read_agent_cfg(argv[1]);

  if(configure_agent(a) != 0)
    error("bad return from configure_agent");

  pr_agent(a);
  exit(0);
}
