#include <R.h>
#include <float.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>	/* for the QR	  routines */
#include <R_ext/Utils.h>	/* for the *sort() routines */
#include <stdio.h>
#include "lower_level.h"
#define BIG DBL_MAX


/* ---------------------------------------------------------------------------------------------
The code below is taken from the R package MASS
*/

static double *coef, *qraux, *work, *res, *yr, *xr, *means, *d2, *d2copy;
static int *pivot, *which, *which2;
static int *ind;



/*
   Sampling k from 0:n-1 without replacement.
 */
static void sample_noreplace(int *x, int n, int k){
    int i, j, nn=n;

    for (i = 0; i < n; i++) ind[i] = i;
    for (i = 0; i < k; i++) {
	j = (int)(nn * unif_rand());
	x[i] = ind[j];
	ind[j] = ind[--nn];
    }
}

/*
   Find all subsets of size k in order: this gets a new one each call
 */
static void next_set(int *x, int n, int k){
    int i, j, tmp;

    j = k - 1;
    tmp = x[j]++;
    while(j > 0 && x[j] >= n - (k - 1 -j)) tmp = ++x[--j];
    for(i = j+1; i < k; i++)  x[i] =  ++tmp;
}



static void mve_setup(int *n, int *p, int *ps) {
    xr = (double *) R_alloc((*ps)*(*p), sizeof(double));
    qraux = (double *) R_alloc(*p, sizeof(double));
    pivot = (int *) R_alloc(*p, sizeof(int));
    work = (double *) R_alloc(2*(*p), sizeof(double));
    d2 = (double *) R_alloc(*n, sizeof(double));
    d2copy = (double *) R_alloc(*n, sizeof(double));
    means = (double *) R_alloc((*p), sizeof(double));
    ind = (int *) R_alloc(*n, sizeof(int));
    which = (int *) R_alloc(*ps, sizeof(int));
    which2 = (int *) R_alloc(*ps, sizeof(int));
}


/* find the squared Mahalanobis distance to x via QR decomposition in xr. */
static double mah(double *xr, int nnew, int p, double *x){
    int i, j;
    double s, ss = 0.0;

    for(j = 0; j < p; j++) {
	s = x[j];
	if(j > 0) for(i = 0; i < j; i++) s -= work[i] * xr[i + nnew*j];
	work[j] = s / xr[j + nnew*j];
	ss += work[j] * work[j];
    }
    return(ss*(nnew-1));
}

/*
   Compute the squared Mahalanobis distances, in d2, to all points in x
   from the mean of the subset in which using the covariance of that
   subset.
*/
static int do_one(double *x, int *which, int n, int nnew, int p, double *det, double *d2){
    int i, j, k;
    int rank;
    double sum, tol = 1.0e-7;

    for(j = 0; j < nnew; j++)
	for(k = 0; k < p; k++) xr[j + nnew*k] = x[which[j] + n*k];
    for(k = 0; k < p; k++) {
	sum = 0.0;
	for(j = 0; j < nnew; j++) sum += xr[j + nnew*k];
	sum /= nnew;
	means[k] = sum;
	for(j = 0; j < nnew; j++) xr[j + nnew*k] -= sum;
    }

    F77_CALL(dqrdc2)(xr, &nnew, &nnew, &p, &tol, &rank, qraux, pivot, work);
    if(rank < p) return(1);

    sum = 0.0;
    for(k = 0; k < p; k++)
	sum += log(fabs(xr[k + nnew*k]));
    *det = sum;

    /* now solve R^T b = (x[i, ] - means) and find squared length of b */
    for(i = 0; i < n; i++) {
	for(j = 0; j < p; j++) qraux[j] = x[i + n*j] - means[j];
	d2[i] = mah(xr, nnew, p, qraux);
    }
    return(0);
}


void mve_fitlots(double *x, int *n, int *p, int *qn, int *mcd,
	    		 int *sample, int *nwhich, int *ntrials,
	    		 double *crit, int *sing, int *bestone) {
    int i, iter, j, nn = *n, quan = *qn, trial, this_sing;
    int nnew = *nwhich;
    double det, best = BIG, thiscrit, lim;

    if(*mcd != 1)
	mve_setup(n, p, nwhich);
    else
	mve_setup(n, p, n); /* could get ties */

    *sing = 0;
    if(!*sample) {
	for(i = 0; i < nnew; i++) which[i] = i;
    } else GetRNGstate();

    thiscrit = 0.0;		/* -Wall */

    for(trial = 0; trial < *ntrials; trial++) {

	R_CheckUserInterrupt();

	if(!(*sample)) {if(trial > 0) next_set(which, nn, nnew);}
	else sample_noreplace(which, nn, nnew);

	/* for(i = 0; i < nnew; i++) printf("%d ", which[i]); printf("\n");
	   fflush(stdout);*/


	/* Find the mean and covariance matrix of the sample. Check if singular.
	   Compute Mahalanobis distances of all points to the means using
	   this covariance matrix V, and find quantile largest. Volume is
	   then monotone in determinant(V * dist2). */

	this_sing = do_one(x, which, nn, nnew, *p, &det, d2);
	if(this_sing)  {(*sing)++; continue;}

	/*for(i = 0; i < nnew; i++) printf(" %d", which[i]); printf("\n");*/

	for(i = 0; i < nn; i++) d2copy[i] = d2[i];
	rPsort(d2copy, nn, quan-1);
	lim = d2copy[*qn-1];
	if(!*mcd) thiscrit = (*p) * log(lim) + 2*det;
	else {
	    for(iter = 0; iter < 4; iter++) {
		/*for(i = 0; i < nn; i++) printf(" %f", d2[i]); printf("\n");*/
		if(iter > 0) {
		    for(i = 0; i < nn; i++) d2copy[i] = d2[i];
		    rPsort(d2copy, nn, quan-1);
		    lim = d2copy[*qn-1];
		}
		j = 0;
		for(i = 0; i < nn; i++)
		    if(d2[i] <= lim) which2[j++] = i;
		/* note: we take all points that meet this limit:
		   there could be more than quan. */
		(void) do_one(x, which2, nn, quan, *p, &det, d2);
		if(iter > 0 && 2 * det >= 0.999*thiscrit) break;
		thiscrit = 2 * det;
		/* printf("iter %d %f", iter, thiscrit);
		   for(i = 0; i < quan; i++) printf(" %d", which2[i]);
		   printf("\n"); fflush(stdout);*/
	    }

	}
	/*   printf("this %f\n", thiscrit);*/


	if(thiscrit < best) { /* warning: first might be singular */
	    best = thiscrit;
	    for(i = 0; i < nn; i++) bestone[i] = (d2[i] <= lim);
	}
    }
    *crit = best;
    if(*sample) PutRNGstate();
    /* mve_free(); */
}
