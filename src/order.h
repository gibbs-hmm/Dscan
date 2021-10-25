#if !defined (ORDER)
#define	ORDER
#include "stdinc.h"
/************************ ADT Order ***********************************
Defines a model for the orderings of ntyp different types of linear 
elements having n (or npos) total positions in each of nseq sequences 
(of elements).

  pos:           0     1     2     3     4     n-1    n
  order:        ---[1]---[2]---[3]---[4]---....---[n]---

  seq[1]:	    A     B     A     C            B
  seq[2]:	    C     B     A     B            A
    :
    :
  seq[nseq]:	    B     C     A     B            A

 model[pos]:     ---[1]---[2]---[3]---....---[n]---
   A		  #A[1] #A[2] #A[3]  ....  #A[n] 
   B		  #B[1] #B[2] #B[3]  ....  #B[n] 
   C		  #C[1] #C[2] #C[3]  ....  #C[n] 

  pseudo:	   A*N0  ...etc...  (all positions have N0 of each element)
		   B*N0
		   C*N0 
N0 defines the prior probability N0 (in pseudocounts) that any 
element type will occur at any position in the ordering.
**********************************************************************/
/*************************** Order type **************************/
typedef struct {
	long	ntyp;		/* number of types of elements */
	long	npos;		/* total number of postions for elements */
	long	nseq;		/* number of sequences in the model */
	double	N0;		/* number of pseudo elements at each pos */
	double	**model; 	/* model[npos][ntyp] = # each type @ pos */
} order_type;
typedef order_type *o_type;

/******************************* private *****************************/
void    order_error(char *s);

/******************************* Public ******************************/
/***************************** operations ****************************/
o_type  Order(long ntyp, long npos, double npseudo);
void	PutOrder(FILE *fptr, o_type R);
long    *ConcensusOrder(o_type R);
void    InitOrder(o_type R);
void    Add2Order(long *order, o_type R);
void    RmOrder(long *order, o_type R);
double	RelProbOrder(long *order, long typ, long pos, o_type R);
o_type  NilOrder(o_type R);

/**************************** macro operations **********************/
#define MAX_NUMBER_MOTIFS	26
#define nMotifsOrder(R)		((R)->ntyp)
#define nPosOrder(R)		((R)->npos)

#endif

