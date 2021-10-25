#if !defined (OLIST)
#define OLIST
#include "stdinc.h"
/**************************  ADT olist ********************************
 Header file for data structure representing ordered list of small 
positive integers.

	list: x  x  x   ... x   where x < x < x < ... < x
               1  2  3       n         1   2   3         n

	and 0 < n <= N.	(note: next[x ] = x   .)
				     i     i+1
**********************************************************************/

typedef struct {
	long	N;		/* list defined on ints in {1,...,N} */
	long	n;		/* number of elements in list */
	long	*next;		/* next[i] is successor of i in list */
} orderlist_type;

typedef orderlist_type *ol_type;

void olist_error(char *s);
/****************************** PUBLIC ***************************/
ol_type	Olist(long N);
long	RmOlist(long i, ol_type L);
long	InsertOlist(long i, ol_type L); /* access item i in THETA(i) time */
void	PutOlist(FILE *fptr, ol_type L);/* print item i on list */
void    GetOlist(long *array, ol_type L);
void	NilOlist(ol_type L);

/********************** MACROS **********************************/

/* Return number of elements in list */
#define LengthOlist(L)	((L)->n)
#define ClearOlist(L)	((L)->next[0]=0)

#endif
