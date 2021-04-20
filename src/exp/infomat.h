/** infomat.h: prototypes for functions in infomat.c **/

real quadratic_form(real *v, real **m, real *u, int dim);
real l_kk(void);
real l_kb(void);
real l_ka(void);
real l_bb(void);
real l_ba(void);
real l_aa(void);
real lDetC_kk(void);
real lDetC_kb(void);
real lDetC_ka(void);
real lDetC_bb(void);
real lDetC_ba(void);
real lDetC_aa(void);
real Q_kk(void);
real Q_kb(void);
real Q_ka(void);
int infomat(real **imat, real **C, **S, int iparam, int P);
