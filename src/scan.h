/* scan.h - codes and constants for scan program. 
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/scan.h,v 2.0 2003/11/20 15:13:54 palumbo Exp $
 *
 */

#if !defined (SCAN)
#define SCAN
/*** #include <time.h> /****/
#include <limits.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "afnio.h"
#include "histogram.h"
#include "probability.h"
#include "sequence.h"
#include "smatrix.h"
#include "smodel.h"
#include "mheap.h"
#include "sites.h"
#include "residues.h"

typedef struct {
	a_type  A;
	char	*snfile;
	Boolean	seqfile,repeats,shuffle,new_way;
	Boolean	neg_mask;
	long	min_segs, max_segs;
	long	N;			/** number of motif models **/
	long	hpsz;
	long	**pos;
	double	total;
	double	*freq;
	double	maxEval,logmaxEval;
	double	singlePval;	/** single Pval min **/
	sm_type  *M;
	mh_type H;
	h_type	HG;
	e_type	*E;
	char	method;		/** method used for model: g,c,r **/
        short   nCountFlag;     /* BT 10/17/97 */
        int     nCount;         /* override -e and print top nCount items */
        int     nNoOverlap;     /* BT 10/4/2000 *//* don't allow overlapping motif with multiple models */
        int     *numSites;      /* effective number of search sites for each model */
  int     adjustN;        /* adjust for actaula size of data base */
} scan_type;
typedef scan_type *scn_typ;

#define MAX_BLOCK_SIZE  10000      /* BT 7/30/2001 */
#define MAX_BLOCK_LENGTH 200
#define MAX_NUM_MODELS	50
#define	ID_LENGTH	200

#define DEFAULT_FREQ_MODEL_NUM_SITES   (100)   /* mjp 10-31-03 */

/********************************* PRIVATE ********************************/
long     foo_scan(long start, long m, long N, long **Items, long **scores, 
		long **site,long *path, long *w, long *tmp, long tmp_score, 
		long *best);

/********************************* PUBLIC ********************************/
scn_typ  MkScan(double pseudo, char *snfile, long *counts, a_type A, 
		long hpsz, char method, short nFlag, int nCount, 
		int nNoOverlap,
		int adjustN,    /* BT 11/06/2001 */
		long *freqModelNumSites,  /* mjp, 11-03-03 */
		int   numNumSitesArgs);   /* mjp, 11-03-03 */


void     NilScan(scn_typ F);
long     MaskScan(FILE *fptr,long number, long *nsize, scn_typ F);    /* BT 9/26/97 */
long     ScanScan(FILE *fptr,long number, long *nsize, scn_typ F);   /* BT 9/26/97 */
long     ScanRScan(FILE *fptr,long number, long *nsize, scn_typ F);  /* BT 9/26/97 */
long     OrderScanScan(FILE *fptr,long number, long *nsize, scn_typ F);  /* BT 9/26/97 */
long     OrderScanMtf(char *name, scn_typ F);
long     PutScan(FILE *fptr, scn_typ F);
long     PutRScan(FILE *fptr, scn_typ F);
long     ReadScan(scn_typ F, double pseudo,char *snfile, long *counts,
		  long *freqModelNumSites,   /* mjp, 11-03-03 */
		  int   numNumSitesArgs);    /* mjp, 11-03-03 */

long     SetMaxEvalScan(double maxEval, scn_typ F);
int      Overlap( int m, int i, long *p, smx_typ *M);
int      CountNumSites( FILE *fptr, long number,  long *nsize, scn_typ F, long *w, int *numSites );


/********************************* MACROS ********************************/
#define OpenSeqFileScan(b,F)	((F)->seqfile = (b))
#define ShuffleSegsScan(F)	((F)->shuffle = TRUE)
#define SetMinSegsRScan(x,F)    ((F)->min_segs = (x))
#define SetMethodScan(x,F)	((F)->method = (char)(x))
#define SetMaxSegsRScan(x,F)    ((F)->max_segs = (x))
#define SetSinglePvalScan(x,F)  ((F)->singlePval = (double)(x))
#define UseNegMaskScan(F)	((F)->neg_mask = TRUE)
#define NumModelsScan(F)    	((F)->N)
#define SMatrixScan(m,F)	(((m) <= (F)->N && (m) > 0)?\
					GetSMatrixSModel((F)->M[(m)]): NULL)

#endif

