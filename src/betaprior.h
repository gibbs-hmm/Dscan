/****************** betaprior.h - .***************

  beta prior distribution for Bayesean statistics:

	mean = m = a/(a+b)
	var = S^2 = a*b/((a+b+1)(a+b)^2)

	a = -m - m^2*(m-1)/s^2
	b = -1 +m + m*(m-1)^2/s^2

*************************************************************************/
#if !defined(BETAPRIOR)
#define BETAPRIOR
#include "afnio.h"
#include "stdinc.h"
#include "probability.h"

/********************************* PRIVATE ********************************/
typedef struct {
	double		A;		/* total pseudo trials */
	double		N;		/* total real trials */
	double		T;		/* Total trials */
	double  	alpha;		/* #successful pseudo trials */
	double		beta;		/* #failed pseudo trials */
	unsigned long	success;	/* #successful real trials */
	unsigned long	expect; 	/* expected number of success */
	double		p;		/* posterior probability of sucess */
	double		weight;		/* A/N = w/(1-w) - fractional weight */ 
	Boolean		calc;
} beta_prior_type;

typedef beta_prior_type *bp_type;

/******************************** private ********************************/
void    bprior_error(char *s);

/********************************* PUBLIC ********************************/
bp_type MkBPrior(long expect, double weight, double N);
void	ClearBPrior(bp_type B);
void	SetBPriorS(long success, bp_type B);
void	SetBPriorN(double N, bp_type B);
void    PutBPrior(FILE *fptr, bp_type B);
double  BPriorMAP(register bp_type B);
double  LogLikeBPrior(bp_type B);
double  BPriorNullMAP(register bp_type B);
double  RatioBPrior(register bp_type B1, register bp_type B2);
/********************************* MACROS ********************************/
#define NilBPrior(B)	free((B))
#define SuccessBPrior(B) ((B)->success)
#define BPriorN(B)	((B)->N)
#define ExpectBPrior(B)	((B)->expect)
#define alphaBPrior(B)	((B)->alpha)
#define betaBPrior(B)	((B)->beta)
#define PostProbBPrior(B)  (((B)->calc)? ((B)->calc)=FALSE,\
			((B)->p=(B->alpha+(double)B->success)/(B->T-1)):\
			(B)->p)
#define AddBPrior(B)	(((B)->calc)=TRUE, (B)->success++)
#define RmBPrior(B)	(((B)->calc)=TRUE, (B)->success--)

#endif

