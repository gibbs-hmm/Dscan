/****************** sites.h - sites abstract data type.***************/
#if !defined(SITES)
#define SITES
#include "stdinc.h"
#include "olist.h"
#include "seqset.h"
#include "dheap.h"
#include "alphabet.h"
#include "random.h"
/*************************** ADT sites *********************************
Defines regions of sequences to different types of elements.  

ntyp = 6; 	t = ( A  B  C  D  E  F)
len_elem = 	      5  5  3  3  5  6
maxinc = min(len_elem) - 1 = 2

type[1]   =    A  F  B  E  D
pos[1]    =    3 18 30 40 50 
 
type[2]   =    E  D  C  B  F  A  D  D
pos[2]    =    6 13 21 24 33 42 52 55

type[3]   =    C  E  A  F  D  B
pos[3]    =    4 10 18 23 32 38

nseq = 3	( 1   2   3)
len_seq =	 57  59  51
nsites[A] =       1   1   1
nsites[B] =       1   1   1
nsites[C] =       0   1   1
nsites[D] =       1   3   1
nsites[E] =       1   1   1
nsites[F] =       1   1   1

type[1]:     ..Aaaaa..........Ffffff......Bbbbb.....Eeeee.....Ddd.....*
type[2]:     .....Eeeee..Ddd.....CccBbbbb....Ffffff...Aaaaa.....DddDdd..*
type[3]:     ...Ccc...Eeeee...AaaaaFfffff...Ddd...Bbbbb.........*
 	     ....:....|....:....|....:....|....:....|....:....|....:....|
                 5   10   15   20   25    30  35   40   45   50   55   60

Note: '.' = VACANT; lowercase  = BLOCKED(t); '*' = ENDTYPE_SITE.

pos_prob(t,n,i) = probability of a t site in sequence n at position i.
	NOTE: pos_prob[t][n][0] = sum of probs. 
 ********************************************************************/

typedef struct {
	ss_type	data;		/* sequence set */
	long	ntyp;		/* ntyp = number of types of elements. */
	long	*len_elem;	/* len_elem(t) = the length of element t */
	long	**nsites;	/* nsites(t,n) = # of type t sites in seq n */
	unsigned short	***site_pos;	/* site_pos[t][n][s]: site positions */
	long	nseq;		/* number of sequences */
	long	maxinc;		/* maxinc = length of shortest element - 1 */
	long	*len_seq;	/* length of sequences */
	long	*totsites;	/* totsites(t) = total # of type t sites */
	long	*tmp;		/* temp. buffer for site positions */
	ol_type *pos;		/* ordered list of site positions in seq n */
	char	**type;		/* type(n,i) = type of ith element in seq n */
	double	***pos_prob;	/* probability of a t site at pos i in seq n*/
} sites_type;
typedef sites_type *st_type;

typedef struct {
	ss_type	data;			/* sequence set */
	long	ntyp;			/* number of types of elements. */
	long	*len_elem;		/* the length of element t */
	long	**nsites;		/* number of type t sites in seq n */
	unsigned short	***site_pos;	/* positions of sites */
} sites_info_type;
typedef sites_info_type *sti_typ;

/********************************* PRIVATE ********************************/
#define VACANT			0
#define BLOCKED(t)		-(t)
#define MAX_NO_TYPESITES	100
#define	ENDTYPESITE		-120

long     print_sites(unsigned short *L);
long     bubble_sort_sites(unsigned short *L);

void    sites_error(char *s);
/********************************* PUBLIC ********************************/
/*-------------------------------- archive -----------------------------*/
sti_typ ArchiveSites(st_type S);
Boolean SameArchiveSites(sti_typ X1, sti_typ X2);
st_type ExtractSites(sti_typ X);
void    NilArchiveSites(sti_typ X);
/*-------------------------------- Sites -------------------------------*/
long     AddRandomSite(long t,long n, long ntries,st_type S);
void    AddSite(long t, long n, long site, st_type S);
long     ChooseSite(long t, long n, st_type S);
st_type	CopySites(st_type S);
void	GrowSites(long t, st_type S);
void    InitSites(st_type S);
double	MissInfoSites(long typ, st_type S);
st_type MkSites(long ntyps, long *len_elem, ss_type data);
void	NilSites(st_type S);
Boolean NRandomSites(long t, long N, long max, st_type S);
Boolean OccupiedSite(register long t, register long n, register long site, 
        register st_type S);
void	PosSites(long n, long *pos, st_type S);
void    PosTSites(long t, long n, long *pos, st_type S);
void    PutScanSites(FILE *fptr, long t, st_type S, Boolean *off);
void    PutSitesMtfDBS(FILE *fptr,long t,st_type S,double **prob,
	Boolean *off);
void    PutSites(FILE *fptr,long t,st_type S,double **site_prob, Boolean *off);
void	PutTypeSites(FILE *fptr, st_type S);
void	OrderSites(long n, long *order, st_type S);
st_type ReadRandomSites(long ntyps, long *len_elem, ss_type P);
void    ShiftSites(st_type S, long t, Boolean left);
void    ShiftSitesM(st_type S, long t, long d);
void	ShrinkSites(long t, st_type S);
Boolean ShuffleSites(st_type S);
void	VacateSite(long t, long n, long site, st_type S);
/********************************* MACROS ********************************/
#define SitesSeqSet(S)		((S)->data)
#define SitePos(t,n,k,S)	((S)->site_pos[(t)][(n)][(k)])
#define TypeSite(n,s,S)		((S)->type[(n)][(s)])
#define StartSite(n,s,S)	((S)->type[(n)][(s)] > 0)
#define nSites(t,n,S)		((S)->nsites[(t)][(n)])
#define TotalSites(t,S)		((S)->totsites[(t)])
#define SiteLen(t,S)		((S)->len_elem[(t)])
#define nTypeSites(S)		((S)->ntyp)
#define OpenPos(n,s,S)		(!(S)->type[(n)][(s)])
#define PosProbSite(t,n,S)	((S)->pos_prob[(t)][(n)])
#define ProbSite(t,S)		((S)->pos_prob[(t)])
#define BlockedSite(t,S)      	-(t)

#endif

