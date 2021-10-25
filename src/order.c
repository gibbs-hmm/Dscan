#include "order.h"

o_type  Order(long ntyp, long npos, double npseudo)
/* create and return a null ordering of ntyp elements in npos positions;
   nseq = 0, only psuedo sequences counts are present (=prior opinion). */
{
	o_type R;
	long	i;
	long	t;

	if(ntyp > MAX_NUMBER_MOTIFS) order_error("too many motifs");
	NEW(R,1,order_type);
	R->ntyp = ntyp;
	R->npos = npos;
	R->nseq = 0;
	R->N0 = npseudo;
	NEWP(R->model,npos+1,double);
	for(i=1; i<= npos; i++){
		NEW(R->model[i],ntyp +1,double);
	}
	InitOrder(R);
	/** fprintf(stderr,"order pseudo counts = %g\n",R->N0); /***/
	return R;
}

void	InitOrder(o_type R)
/* set the order to contain no sequences (except pseudo sequences) */
{
	long i,j;
	for(i=1; i<=R->npos; i++){
		for(j=1; j<=R->ntyp; j++){
			R->model[i][j] = R->N0;
		}
	}
	R->nseq = 0;
}

void	Add2Order(long *order, o_type R)
/* add a sequence with ordering  *order to the model; assumes that 
   order is an array of npos elements (integers). */
{
	long i;
	R->nseq++;
	for(i=1; i<=R->npos; i++){
		R->model[i][order[i]] += 1.0;
	}
}

void	RmOrder(long *order, o_type R)
/* Remove a sequence with ordering *order from the model;
   assumes that order is an array of npos elements (integers). */
{
	long i;
	R->nseq--;
	for(i=1; i<=R->npos; i++){
		R->model[i][order[i]] -= 1.0;
	}
}


double	RelProbOrder(long *order, long t, long pos, o_type R)
/**********************************************************
   Return the relative probability that the order (given by inserting 
   an element of type t at position i in *order) belongs to the 
   model.  (*order is assumed to be an array of n-1 elements.) 

  pos=2:		    t=A
  type:		    B     A     C     B             A
  *order:	---[1]---[2]---[3]---[4]---....---[n-1]---
  position:   	 0     1     2     3     4     n-2      n-1

  creates:	---[1]---[2]---[3]---[4]---[5]....---[n]---
  type:		    B     A     A     C     B         A

***********************************************************/
{
	long i;
	double	P;	/* probability */

	if(pos >= R->npos) order_error("not that many positions in order");
	for(P=1.0,i=1; i <= R->npos; i++){
		if(i-1 == pos){
			P *= R->model[i][t];
		} else if(i-1 < pos){
			P *= R->model[i][order[i]];
		} else {  /* i-1 > pos */
			P *= R->model[i][order[i-1]];
		}
	}
	return P;
}

long	*ConcensusOrder(o_type R)
/********************************************************************
  Return an array of the consensus order for the types.

    type |   1   2
   ------+--------
       A |   7  16
       B |  16   7
	    
  order: B A.
 ********************************************************************/
{
	long i,j,k,*order,n,max_typ,t,*array,N,jmax;
	double	max;

	N = R->npos;
	NEW(order,N+3,long);
	printf("npos = %d; ntyps = %d\n",N,R->ntyp);
	if(N != R->ntyp) {
		for(i=1; i<=N; i++) order[i] = i;
		return order;
	}
	NEW(array,N+3,long);
	for(i=1; i<=N; i++) array[i] = i;
	for(n=N,i=1; i<=N; i++){
		for(max=0.0,j=1; j<=n; j++){
		   t = array[j];
		   if(max < R->model[i][t]){
			max = R->model[i][t];
			max_typ = t;
			jmax = j;
		   }
		}
	     	order[i] = max_typ;
		array[jmax] = array[n]; n--;
	}
	for(i=1; i<=N; i++) printf("%c ", 'A' + order[i] - 1);
	printf("\n");
	free(array);
	return order;
}

void	PutOrder(FILE *fptr, o_type R)
/*  Print out a table showing the number of elements at each type 
    in the model */
{
	long i,j,k;

	fprintf(fptr,"\n%5s |", "type");
	for(i=1; i<=R->npos; i++){ fprintf(fptr,"%4d", i); }
	fprintf(fptr,"\n------+");
	for(i=1; i<=R->npos; i++){ fprintf(fptr,"----"); }
	fprintf(fptr,"\n");
	for(j=1; j<=R->ntyp; j++){
	     fprintf(fptr,"%5c |", j+'A'-1);
	     for(i=1; i<=R->npos; i++){
		fprintf(fptr,"%4d", (long)(R->model[i][j]-R->N0+0.00001));
	     }
	     fprintf(fptr,"\n");
	}
}

o_type  NilOrder(o_type R)
/*  Destroy Order data object R. */
{
	long i;
	for(i=1; i<= R->npos; i++){ free(R->model[i]); }
	free(R->model);
	free(R);
	return (o_type) NULL;
}

void	order_error(char *s){ fprintf(stderr,"Order: %s\n",s); exit(1); }


