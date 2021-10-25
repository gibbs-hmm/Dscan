#if !defined (SMATRIX)
#define	SMATRIX
#include "stdinc.h"
#include "alphabet.h"
#include "probability.h"
/*************************** ADT PROFILE ***************************
	
			scoring matrix data type

**********************************************************************/

/*************************** SMatrix Type **************************/
typedef struct {
	long	*cmax,*max;		/* maximum */
	long	*cmin,*min;		/* minimum */
	long	K;			/* length */
	char	*maxseg;		/* maximum scoring segment */ 
	long	**score;		/* scoring matrix */
	double	mean,var,sd;
	double	nsd;			/* minimum # std. dev. above mean */
	double	*f,*f0,*fx;
	double	*freq;			/* freq[r] */
	a_type	A;			/* alphabet */
	long	nlet;
	long	neginf;			/* lowest permissible score */
	Boolean	changed,calc_stats,calc_prob;	/* update data */
} smatrix_type;
typedef smatrix_type *smx_typ;

/******************************* private *****************************/
Boolean smatrix_prob(register long k, register double *f, smx_typ M);
long     min_max_smatrix(smx_typ M);
long     stats_smatrix(smx_typ M);
double  smatrix_prob_fast(register long k, register double *f, smx_typ M);
void    smatrix_error(char *s);

/******************************* Public ******************************/
/***************************** operations ****************************/
void    PutSMatrix(FILE *fptr, smx_typ M);
smx_typ MkSMatrix(double nsd, long K, double *freq, a_type A);
void    NilSMatrix(smx_typ M);
long     MaxScoreSMatrix(smx_typ M);
char    *MaxSegSMatrix(smx_typ M);
long     SetSMatrix(long r, long row , long score, smx_typ M);
long     ScoreSMatrix(register char *seq, register long start,
			register smx_typ M);
long     SplitScoreSMatrix(char **seq, long n, long *start, long *leng, smx_typ M);
double  SMatrixProb(long score, smx_typ M);
double  SMatrixProbFast(long score, smx_typ M);
/**************************** macro operations **********************/
#define kSMatrix(M)		((M)->K)
#define meanSMatrix(M)	(((M)->changed)? stats_smatrix(M),(M)->mean:\
				(M)->mean)
#define sdSMatrix(M)	(((M)->changed)? stats_smatrix(M),(M)->sd:(M)->sd)
#define ValSMatrix(pos,r,M)	((M)->score[(pos)][(r)])
#define NegInfSMatrix(M)	((M)->neginf)

#endif

