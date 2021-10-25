/************** smodel.h - score model abstract data type.***********
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/smodel.h,v 2.0 2003/11/20 15:14:10 palumbo Exp $
 *
 */



#if !defined(SMODEL)
#define SMODEL
#include "stdinc.h"
#include "alphabet.h"
#include "smatrix.h"
#include "blosum62.h"
#include "random.h"
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
S(16): 1071   0.086                     length = 5
T(17): 704    0.056                     npseudo = 2
W(18): 158    0.013                     N0[r] = npseudo * freq[r]
Y(19): 434    0.035
V(20): 730    0.058
tot_cnts: 12509

        site_freq[pos][r] = N0[r] + counts[r]
        pos   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
     |  1    0  0  0  0  1  0  0  0  7  0  0  0  0  0  1  0  0  0  1  0
     |	2 (null)
     |	3 (null)
     |  4   1  5  0  0  0  0  0  1  0  0  0  0  0  1  1  1  0  0  0  0
     |  5   0  0  0  9  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
     |  6 (null)
     |  7 (null)
     |  8   0  0  0  0  0  0  0  0  0  5  4  0  0  0  0  0  0  0  0  1
     |  9 (null)
     |  10   1  0  0  0  1  0  0  0  0  0  0  6  0  0  0  1  1  0  0  0

        site_freq0[r] (nonsite model) == # r residues not in model + N0[r].

        site_freqN[0][r] = site_freq[0][r]/sum(site_freq[0][r=A..V]-N0[r=A..V])

 	factor = totsites 

        site_freqN[pos][r] = site_freq[pos][r]/(factor*site_freqN[0][r])

********************************************************************/

typedef struct {
	a_type	A;		/* alphabet */
	smx_typ	smx;		/* smatrix for p-values */
	double	maxscore;	/* maximum possible score */
	long	length;		/* length of motif model */
	long	**site_freq;	/* site_freq(j,b)=# b's at pos j of site*/
	double	**likelihood;	/* likelihood ratio for being in model */
	double	pseudo;		/* pseudo = base pseudo count priors */
	double	npseudo;	/* npseudo = # pseudo priors */
	double	*N0;		/* N0(nres) = # pseudo residues by type */
        double	*freq;          /* freq(b)= total frequency of b's */ 
	double	*temp;		/* temp(nres) = temp double array */
	double	*tmp_val;	/* temp value(maxlen) = temp double array */
	long	totsites;	/* number of sites in model */
	long	tot_cnts;	/* total number of residues */
	long	*counts;	/* total number of b residues in sequences */
	char	method;		/* method to calc Staden pval (default: g)*/
	Boolean	*null;		/* null site or not?? */
	Boolean	update;		/* update the normalized freq */
} smodel_type;
typedef smodel_type *sm_type;

/********************************* PRIVATE ********************************/
#define SMODEL_MAX_SCORE 1000.0
void    update_smodel_freqN(sm_type M);
void    get_smx_smodel(register sm_type M);
void	smodel_error(char *s);

/********************************* PUBLIC ********************************/
sm_type SModel(long length,double npseudo,long *counts,a_type A);
sm_type MkSModel(Boolean *null, long length,double npseudo,long *counts,a_type A);
void    InitSModel(sm_type M);
void    Add2SModel(char *seq, long site, sm_type M);
void    RmSModel(char *seq, long site, sm_type M);
sm_type  NilSModel(sm_type M);
char    *MaxSegSModel(sm_type M);
double	RelProbSModel(register char *seq, register long pos, register sm_type M);
long	ScoreSModel(register char *seq, register long pos, register sm_type M);
double  PvalSModel(register char *seq, register long pos, register sm_type M);
smx_typ GetSMatrixSModel(sm_type M);
double  ProbSModel(register char *seq, register long pos, register double p,
        register sm_type M);
double	PutSModel(FILE *fptr, sm_type M);
Boolean NullSiteSModel(long s,sm_type M);
long	CellScoreSModel(long r, long pos, sm_type M);
double  VarianceInfoSModel(double *average, long N, sm_type M);
double  InfoSiteSModel(long site, sm_type M);
double  InfoSModel(sm_type M);
double  InfoSModel2(sm_type M);
double  **RealScoresSModel(sm_type M);
long    DrawSeqSModel(char *seq, sm_type M);
void    SetSModelFromProbCounts(sm_type M, double **modelProbMatrix, 
			     long modelNumSites);


/********************************* MACROS ********************************/
#define LenSModel(M)		((M)->length)
#define SetMethodSModel(x,M)	((M)->method = (char) (x))
#define NumResSModel(j,b,M)	((M)->site_freq[(j)][(b)])
#define IsColumnSModel(j,M)	(!(M)->null[(j)])
#define TotSitesSModel(M)	((M)->totsites)

#endif
