#ifndef lower_h
#define lower_h



static void sample_noreplace(int *x, int n, int k);

static void next_set(int *x, int n, int k);

static void mve_setup(int *n, int *p, int *ps);

static double mah(double *xr, int nnew, int p, double *x);

static int do_one(double *x, int *which, int n, int nnew, int p, double *det, double *d2);

void mve_fitlots(double *x, int *n, int *p, int *qn, int *mcd,
	    		 int *sample, int *nwhich, int *ntrials,
	    		 double *crit, int *sing, int *bestone);

#endif