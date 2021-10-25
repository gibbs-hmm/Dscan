#include "betaprior.h"

bp_type MkBPrior(long expect, double weight, double N)
/*******************************************************************
  N = total sites; 	- Pseudo and Real sites.
  tot_sites  = number of sites;
  alpha = number site pseudo counts;
  beta = number site pseudo counts;
 *******************************************************************/
{
	bp_type	B;

	if(N <=0) bprior_error("total sites must be > 0");
	NEW(B,1,beta_prior_type);
        B->expect = expect;
        B->weight = weight; 
	SetBPriorN(N, B);
	B->success = 0;
        B->p = (B->alpha + (double)B->expect)/B->T;
	B->calc = FALSE;
	return B;
}

double	BPriorNullMAP(register bp_type B)
{
	register double      n;
        return LnBeta(B->alpha,B->N + B->beta) - LnBeta(B->alpha,B->beta);
}

double	RatioBPrior(register bp_type B1, register bp_type B2)
/** Return the ratio of new to old. **/
{
	register double      v1,v2;

        v1 = (double) B1->success;
        v2 = (double) B2->success;
        v1 = LnBeta(v1+B1->alpha,B1->N - v1 + B1->beta)
		- LnBeta(B1->alpha,B1->beta);
        v2 = LnBeta(v2+B2->alpha,B2->N - v2 + B2->beta)
		- LnBeta(B2->alpha,B2->beta);
	return exp(v1 - v2);
}

double	BPriorMAP(register bp_type B)
{
	register double      n;

        n = (double) B->success;
        return LnBeta(n+B->alpha,B->N - n + B->beta) - LnBeta(B->alpha,B->beta);
}

double	LogLikeBPrior(bp_type B)
/** in bits of information **/
{
	double	n,p,L,N;

        n = (double) B->success;
	N = (double)B->N;
	p = n/N;
	L = n*log(p) + (N-n)*log(1.0 - p);
        return 1.4427*L;
}

void    SetBPriorN(double N, bp_type B)
{ 
        double  ratio;

        ratio = (B->weight/(1.0 - B->weight));
	B->N = N;
        B->A =N*ratio;		/* A = N*(w/(1-w)) */
        B->alpha = (double) B->expect*ratio;
        B->beta= B->A - B->alpha;
        B->T=(B->A + B->N);
	B->calc = TRUE; 
}

void    SetBPriorS(long success, bp_type B)
{ B->success = success; B->calc = TRUE; }

void    ClearBPrior(bp_type B)
{
	B->success = 0;
        B->p = (B->alpha + (double)B->expect)/B->T;
        B->calc = FALSE;
}

void	PutBPrior(FILE *fptr, bp_type B)
{
        double  variance,stdev;

        variance = (B->alpha + (double)B->expect);
        variance *= (B->beta + B->N - (double)B->expect);
        variance /= (B->A+B->N)*(B->A+B->N)*(B->A+B->N+1.0);
        stdev = sqrt(variance);
        stdev *= 2*B->N;
        fprintf(stdout, "%d (+/-%4.1f) out of %.0f ",
                        B->expect, stdev,B->N);
        fprintf(stdout,"  a = %g; b = %g; p = %g\n",
                B->alpha, B->beta, B->p);
}

void    bprior_error(char *s) {fprintf(stderr,"betaprior: %s\n",s); exit(1);}

