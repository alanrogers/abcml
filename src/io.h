/** io.h: prototypes for functions in io.h **/

struct wgt *read_wgt(char *fname);
struct wgt *alloc_wgt(void);
struct bonedef *alloc_bonedef(void);
struct agent *alloc_agent(void);
struct counts *alloc_counts(void);
struct bonedef *read_bonedef(char *fname, real *scale_factor);
void pr_bonedef(struct bonedef *b);
void pr_agent(struct agent *a);
double sfun(double x);
char *copystring(char *a, int size);
struct agent *read_agent_cfg(char *fname);
struct counts *read_counts(char *fname);
void pr_counts(struct counts *c);
void prfmat(real **A, int nrow, int ncol, const char *label);
void primat(int **A, int nrow, int ncol, char *label);
void prfvec(real * x, int n, const char *label);
void privec(int *x, int n, char *label);
int configure_agent(struct agent *a);
FILE *mustopen(char *name, const char *mode);
int lbls_differ(char **lbl1, char **lbl2, int len);
void check_labels(int nagent, struct agent **agent,
		  struct bonedef *bones,
		  struct counts *counts);
