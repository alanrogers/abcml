#include <stdio.h>
#include <stdlib.h>
#include "mytypes.h"
#include "bone.h"

struct bonedef *alloc_bonedef(void);
struct bonedef *read_bonedef(char *fname);
void pr_bonedef(struct bonedef *b);

void main(int argc, char **argv)
{
  struct bonedef *b;

  if(argc != 2)
  {
    fprintf(stderr,"usage: xbonedef inputfile.bdf\n");
    exit(1);
  }

  b = read_bonedef(argv[1]);

  pr_bonedef(b);
  exit(0);
}
