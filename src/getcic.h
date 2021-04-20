/****************************************************************
    getcic.h: Prototypes for functions in getcic.c
****************************************************************/
int getcic(FILE *fp);
char *getwordic(char *buff, int bufsiz, FILE * ifp);
int getrealic(real *x, char *buff, int bufsiz, FILE *ifp);
int getintic(int *i, char *buff, int bufsiz, FILE *ifp);
