/****
A Gibbs Sampler algorithm for finding multiple sites in multiple sequences 
****/
#include "gibbs.h"

void	PutWilcoxBGibbs(FILE *fptr,long first, long last, long t, gs_type G)
/* Wilcoxon using MAP sampled probabilities and sites */
{
	ss_type	P=G->data;
	long	n,s,N;
	double	***prob,p;
	e_type	E;
   
	if(G->readprob==NULL) print_error("Wilcox() error");
	prob = G->readprob; N = NSeqsSeqSet(P);
	fprintf(fptr,"\n\n");
	first = MAX(long,1,first); first = MIN(long,first,N);
	last = MAX(long,1,last); last= MIN(long,last, N);
	for(n=first; n<= last; n++){
	   E=SeqSetE(n,P);
	   for(s=1; s<=(long)LenSeq(E); s++){
	     p=prob[t][n][s];
	     if(p >= G->readcutoff){
		if(first == 1) fprintf(fptr,"0 %.4f\n",p);
		else fprintf(fptr,"%.4f 0\n",p);
	     }
	   }
	}
	fprintf(fptr,"\n\n");
}

