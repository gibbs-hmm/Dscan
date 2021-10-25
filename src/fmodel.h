/************** fmodel.h - fragmented model abstract data type.***********/
#if !defined(FMODEL)
#define FMODEL
#include "stdinc.h"
#include "alphabet.h"
#include "smatrix.h"
#include "random.h"
#include "dheap.h"
#include <float.h>

/**************************** model ADT *****************************
        Product multinomial model for an ungapped element allowing 
	fragmentation.

r:  counts[r] freq[r]           
X( 0): 0      0.000             (X is a dummy residue)
A( 1): 806    0.064             freq[r] = counts[r]/tot_cnts
R( 2): 756    0.060             seqs:   H..GD..I.K
N( 3): 629    0.050                     H..AD..L.K
D( 4): 671    0.054                     H..RD..L.K
C( 5): 241    0.019                     H..RD..V.K
Q( 6): 577    0.046                     H..FD..I.T
E( 7): 701    0.056                     H..SD..I.S
G( 8): 691    0.055                     Y..RD..I.K
H( 9): 336    0.027                     C..RD..I.C
I(10): 680    0.054                     H..RD..L.A
L(11): 1376   0.110                     P..PE..L.K
K(12): 615    0.049			----+----+
M(13): 235    0.019			1   5   10
F(14): 471    0.038
P(15): 627    0.050                     totsites = 10
S(16): 1071   0.086                     length = 10
T(17): 704    0.056                     npseudo = 2
W(18): 158    0.013                     Ps[r] = npseudo * freq[r]
Y(19): 434    0.035
V(20): 730    0.058
tot_cnts: 12509

observed[pos][r]:

        pos     A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
     |   0    (nonsite observed frequencies)
     |   1    (null)
     :
     |  x-1   (null)
start+-> x      0  0  0  0  1  0  0  0  7  0  0  0  0  0  1  0  0  0  1  0
     |  x+1   (null)
     |	x+2   (null)
     |  x+3     1  5  0  0  0  0  0  1  0  0  0  0  0  1  1  1  0  0  0  0
     |  x+4     0  0  0  9  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
     |  x+5   (null)
     |  x+6   (null)
     |  x+7     0  0  0  0  0  0  0  0  0  5  4  0  0  0  0  0  0  0  0  1
     |  x+8   (null)
 end +->x+9     1  0  0  0  1  0  0  0  0  0  0  6  0  0  0  1  1  0  0  0
     |  x+10  (null)
     :
     | maxlen (null)

********************************************************************/

typedef struct {
	a_type	A;		/* alphabet */
	long	length;		/* length of motif model */
	long	ncols;		/* number columns in motif model */
	long	maxlen;		/* maximum possible length of model */
	long	start;		/* pointer to first column */
	long	end;		/* pointer to last column */
	long	*observedS;	/* observedS(b) counts for all sites */
	long	**observed;	/* observed(j,b)=# b's at pos j of site*/
	double	**target;	/* target frequencies for model */
	double	**likelihood;	/* likelihood ratio for being in model */
	double	npseudo;	/* npseudo = # pseudo priors */
	double	totalS,total0;	/* total counts = real + pseudo counts */
	double	*Ps;		/* Ps(nres) = # pseudo residues by type */
        double	*freq;          /* freq(b)= total frequency of b's */ 
        double	*sumlgamma;     /* sumlgamma(i)= column i sum of lgamma b's */ 
	double	*tmp_val;	/* temp value(maxlen) = temp double array */
	long	totsites;	/* number of sites in model */
	long	tot_cnts;	/* total number of residues */
	long	*counts;	/* total number of b residues in sequences */
	Boolean	update,recalc;  /* update the normalized freq */
	Boolean	*mark;		/* for marking columns */
} fmodel_type;
typedef fmodel_type *fm_type;

#define FMODEL_UNDEF	-9999
/********************************* PRIVATE ********************************/
void    ClearMarksFModel(fm_type M);
long	center_model(long site, fm_type M);
void	fmodel_error(char *s);
void    update_fmodel(fm_type M);

/********************************* PUBLIC ********************************/
long	AddColumnFModel(long *observed, long pos, fm_type M);
long     AddColumnFModel2(long *observed, long pos, fm_type M);
void    Add2FModel(char *seq, long site, fm_type M);
long     ChoiceLemonFModel(fm_type M);
long     ChoiceOrangeFModel(fm_type M);
fm_type CopyFModel(fm_type M);
Boolean GetSegFModel(char *seg, fm_type M);
double  InfoColFModel(long *observed, fm_type M);
double  InfoFModel(fm_type M);
void    InitFModel(fm_type M);
long     LemonFModel(fm_type M);
double	LikelihoodFModel(register char *seq, register long pos, 
	register fm_type M);
double  LnMapFModel(fm_type M);
double  LnMapFModelNull(fm_type M);
double  LogLikeFModel(fm_type  M);
fm_type MkFModel(Boolean *null,long length,double npseudo,long *counts,a_type A);
long     MvColumnFModel(long *observed, long lemon, long pos, fm_type M);
fm_type  NilFModel(fm_type M);
Boolean NonSiteSegFModel(char *seg, fm_type M);
Boolean NullSiteFModel(long s,fm_type M);
long     NullSitesFModel(Boolean *null, fm_type M);
long     OrangeFModel(fm_type M);
double  ProbFModel(register char *seq, register long pos, register double p,
        register fm_type M);
void    PutFModel(FILE *fptr, fm_type M);
void    PutLogFModel(FILE *fptr, fm_type M);
long	*RandColFModel(fm_type M);
double  RatioFModel(long *observed, long d, fm_type M);
Boolean	RmColumnFModel(long pos, fm_type M);
long	RmColumnFModel2(long pos, fm_type M);
void    RmFModel(char *seq, long site, fm_type M);
void    SetPseudoFModel(double npseudo, fm_type M);
void    ShiftFModel(long *observed, Boolean left, fm_type M);

/********************************* MACROS ********************************/
#define IsColumnFModel(j,M)	((M)->observed[(j+(M)->start - 1)]!=NULL)
#define LenFModel(M)		((M)->length)
#define ContigFModel(M)		((M)->ncols==(M)->length)
#define nPseudoFModel(M)	((M)->npseudo)
#define nColsFModel(M)		((M)->ncols)
#define TotSitesFModel(M)	((M)->totsites)

#endif

