/* olist.c 
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/olist.c,v 2.0 2003/11/20 15:13:31 palumbo Exp $
 *
 */

#include "olist.h"

ol_type	Olist(long N) 
/* Create an ordered list of zero items. */
{
	ol_type L;
	NEW(L,1,orderlist_type);
	L->N = N; 
	L->n = 0;
	NEW(L->next,N+1,long);
	L->next[0] = 0;		/* 0 = null */
	return L;
}

void	NilOlist(ol_type L) { free(L->next); free(L); }

long	RmOlist(long i, ol_type L)
/* Remove item i from list L.  */
{
	long	x=0;
	do{
		if( L->next[x] == i){
			L->next[x] = L->next[L->next[x]];
			L->n--;
			return i;
		} else x =  L->next[x];
	} while( x != 0);
	return 0;
}

long	InsertOlist(long i, ol_type L)
/* Add item i to list L in theta(i) time. */
{
	long x,t;
	if(i > L->N) olist_error("item larger than N");
	for(x=0; L->next[x] != 0; x =  L->next[x]) {
	   if(L->next[x] >= i){
		if(L->next[x] == i) return 0;
		break;
	   }
	}
	t = L->next[x]; L->next[x] = i; L->next[i] = t;
	L->n++;
	return i;
}

void	GetOlist(long *array, ol_type L)
/*  put the items into the array (starting from i=1) as ordered on list L. */
{
	long x,i;
	for(x=0,i=1; L->next[x] != 0; x=L->next[x],i++)
			array[i]=L->next[x];
	array[i] = 0;
}

void PutOlist(FILE *fptr, ol_type L)
/* Print the list L. */
{
	long i;
	for(i=0; L->next[i] != 0; i=L->next[i])
		fprintf(fptr,"%ld ", L->next[i]);
}

void olist_error(char *s){fprintf(stderr,"olist: %s\n",s); exit(1); }


