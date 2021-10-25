/**************************** gibbs.h - *******************************
   Gibbs Sampling algorithms for local multiple alignment.
*************************************************************************/
#if !defined(GIBBS)
#define GIBBS
#include "time.h"
#include "afnio.h"
#include "stdinc.h"
#include "fmodel.h"
#include "betaprior.h"
#include "sites.h"
#include "dheap.h"
#include "mheap.h"
#include "seqset.h"
#include "histogram.h"
#include "probability.h"
#include "residues.h"
#include "order.h"
#include "random.h"

/********************************* PRIVATE ********************************/
typedef struct {
        ss_type		data;			/** input data **/
	char		*name;			/** input file name **/
        a_type		A;			/** alphabet used **/
	Boolean		test[10];		/** test new items -T option **/
	unsigned long	seed;			/** random seed **/
	/****** sites ****/
        st_type		sites;			/** working sites **/
	long		ntyps,*order;
	long		*site_len,*pos;
	sti_typ		best,map,start;		/** archived sites **/
	/****** fmodel ****/
        fm_type		*model;
	long		*ncol,*maxlen;
	Boolean 	*null,**bestnull,**mapnull;
        double		pseudo;
	double		qseudo;
	/****** sampling parameters ********/
        long		nruns,nconverge,ncycles;
        long		**tmpfreq;
        double		*tmpratio;
	Boolean		move,fragment,verbose;
	long		wilcoxon;		/* = 0 if not wilcoxon */
        FILE		*fptr,*ifptr,*sfptr,*mfptr;
	/******** fgibbs items *********/
	long		limit;
	Boolean		use_order;
	/******** motif sampler items *********/
	long		stop;
	double		***readprob,readcutoff,*sumprob,nread;
	double		*p, *map_p, *best_p;
	bp_type		*prior;
	long		*expect;
	double		weight;
} gibbs_sampler_type;

typedef gibbs_sampler_type *gs_type;

/******************************** private ********************************/
Boolean move_column_gibbs(gs_type G, long lemon, long t);
long	*GetSiteFreq(gs_type G, long t, long d);
fm_type *InitGibbs(gs_type G);
fm_type *InitMAPGibbs(gs_type G);
Boolean MoveColumn(gs_type G, long t);
Boolean MoveMultiCols(gs_type G, long t, long num);
Boolean Metropolis(gs_type G, long t, fm_type M);
long     SaveBestGibbs(gs_type G);
Boolean	SaveFinalGibbs(gs_type G);
Boolean ShiftGibbs(gs_type G, long t, fm_type M);
Boolean OptionsGibbs(long argc, char *argv[], gs_type G);
long     GetFreqProb(long t, long n, fm_type M, gs_type G, o_type R);
double	motif_sampler(gs_type G);
double  MapBGibbs(gs_type G);
double  NetMapBGibbs(gs_type G);
double  NetMapBGibbs1(long t, fm_type M, gs_type G);
double  NullMapBGibbs(gs_type G);
Boolean TransferColumn(gs_type G, long ntyp, fm_type *model);

/********************************* PUBLIC ********************************/
gs_type MkGibbs(long nopt, char *options[], st_type S);
st_type StartSitesGibbs(long argc, char *argv[]);
long     RunGibbs(FILE *fptr, gs_type G);
void    PutWilcoxBGibbs(FILE *fptr,long first, long last, long t, gs_type G);
void    NilGibbs(gs_type G);

/**************************** site sampler *******************************/
void	SiteSampler(FILE *fptr, gs_type G);

/***************************** motif sampler *****************************/
long	MotifSampler(FILE *fptr, gs_type G);

/*********************** near optimum sampler ****************************/
double  ***NearOptimumSampler(long niter, gs_type G);
long     NearOptFit(FILE *fptr, gs_type G);

/********************************* MACROS ********************************/
#define GIBBS_USAGE0 "\nUsage(sites sampler): gibbs file lengths [options] \n\
\nUsage(motif sampler): gibbs file lengths expect [options] \n\n\
  lengths = <long>[,<long>]: lengths of elements for each type\n\
  expect = <long>[,<long>]: expected number of elements for each type \n\
  <long>[,<long>] = numbers for each element (e.g. \"10,12,15\" for 3 types)\n\
  options:\n\
\t-C<real>    - prob. cutoff (0 < C <= 1) for near optimum sampling\n\
\t-c<int>     - number of cycles between shifts (sites sampler)\n\
\t              or maximum number of cycles per run (motif sampler) \n\
\t-d          - DON'T use fragmentation (i.e., column sampler)\n\
\t-f          - create a scan output file (file.sn)\n\
\t-I          - interactively specify #sites/sequence (sites sampler)\n\
\t              (default: one site of each type per sequence)\n\
\t-L<int>     - set rapid convergence limit (higher = longer to converge)\n\
\t-m<int>     - set maximum number of cycles in each run (sites sampler)\n\
\t-n          - use nucleic acid alphabet\n\
\t-o          - use element order in probabilities (sites sampler)\n\
\t              (each sequence must contain the same number of elements)\n\
\t-p<float>   - number of pseudo counts for product multinomial model\n\
\t-q<float>   - pseudo counts for ordering model (sites sampler)\n\
\t-R<int>     - set number of near-optimum readings taken\n\
\t-r          - randomly shuffle input sequences\n\
\t-s<int>     - give seed for random number generator\n\
\t-t<int>     - maximum number of sampling runs\n\
\t-x          - DON'T remove protein low complexity regions\n\
\t-w          - output wilcoxon rank test information (motif sampler)\n\
\t-W<float>   - set fractional weight (0 to 1.0) on priors (motif sampler)\n\
\n"

#endif

