/****
A Gibbs Sampler algorithm for finding multiple sites in multiple sequences 
****/
#include "gibbs.h"

gs_type MkGibbs(long nopt, char *options[], st_type S)
{
	gs_type G;
	ss_type	data=SitesSeqSet(S);
	a_type	A=SeqSetA(data);
	long	t,n,max;
	extern char GIBBS_USAGE[];

	NEW(G,1,gibbs_sampler_type);
	G->sites = S;           G->A = A;               G->data = data; 
    	G->ntyps=nTypeSites(S);
	G->best = NULL;         G->map = NULL; 		G->expect = NULL;
	G->nruns=10;            G->nconverge = 500;
	G->verbose = FALSE;
	for(n=0; n< 10;  n++) G->test[n]=FALSE;
	G->wilcoxon = 0;
	G->name=NULL;
	G->ncycles = 1;
	G->pseudo = 0.10;      G->qseudo = 0.05; 
	G->nread = 500;
	G->start = NULL;
	G->fragment = TRUE; 
	G->move = FALSE;
	G->use_order = FALSE;
	G->weight=0.8;
	G->ifptr=NULL;         G->fptr=NULL;           G->sfptr=NULL;
	G->mfptr = NULL;
	G->readprob = NULL; 
	G->sumprob = NULL; 
	G->readcutoff = 0.5;
	G->model= NULL; 
	G->seed = 0;
	G->p = NULL;            G->map_p = NULL;        G->best_p = NULL;
	G->start = NULL; 
	G->prior = NULL;
	G->limit = 10;
	for(G->stop=0, t = 1; t <= G->ntyps; t++) { 
		G->stop += NSeqsSeqSet(data) * G->limit;
	}
	if(!OptionsGibbs(nopt, options, G)) print_error(GIBBS_USAGE);
	NEW(G->order,MaxSeqSeqSet(data),long);
	NEW(G->pos,MaxSeqSeqSet(data),long);
	NEW(G->null, MaxSeqSeqSet(data)+1,Boolean);
	NEW(G->ncol, G->ntyps +1, long);
	NEW(G->maxlen,G->ntyps +1,long);;
	NEW(G->site_len, G->ntyps+1, long);
	NEWP(G->bestnull, G->ntyps+1, Boolean);
	NEWP(G->mapnull, G->ntyps+1, Boolean);
	for(max=0,t = 1; t <= G->ntyps; t++){
	        MEW(G->bestnull[t], MaxSeqSeqSet(data)+1, Boolean);
	        MEW(G->mapnull[t], MaxSeqSeqSet(data)+1, Boolean);
		G->ncol[t] = G->site_len[t] = SiteLen(t,S);
		G->maxlen[t] = 5*G->ncol[t];
		G->maxlen[t] = MIN(long,MaxSeqSeqSet(data),G->maxlen[t]);
		max= MAX(long,G->maxlen[t],max);
	}
	NEWP(G->tmpfreq, max+5,long);
	NEW(G->tmpratio, max+5,double);
	return G;
}

void	NilGibbs(gs_type G)
{
	long	t,n,N;

	N = NSeqsSeqSet(G->data);
   	if(G->name != NULL)  free(G->name);
	if(G->sfptr != NULL) fclose(G->sfptr);
	if(G->p != NULL) free(G->p);
	if(G->map_p != NULL) free(G->map_p);
	if(G->best_p != NULL) free(G->best_p);
	if(G->sumprob != NULL) free(G->sumprob);
	if(G->readprob != NULL){
	   for(t=1; t<=G->ntyps; t++){
		for(n=1; n<=N; n++) free(G->readprob[t][n]);
		free(G->readprob[t]);
	   } free(G->readprob);
	}
	for(t = 1; t <= G->ntyps; t++){
	   if(G->bestnull[t]!=NULL) free(G->bestnull[t]);
	   if(G->mapnull[t]!=NULL) free(G->mapnull[t]);
	   if(G->model[t] != NULL) NilFModel(G->model[t]);
	   if(G->prior != NULL) NilBPrior(G->prior[t]);
	}
	if(G->model != NULL) free(G->model);
	if(G->order!= NULL) free(G->order); 
	free(G->pos);
	if(G->sites != NULL) NilSites(G->sites);
	free(G->null); free(G->tmpfreq); free(G->tmpratio);
	free(G->ncol); free(G->maxlen); 
	free(G->site_len); 
	free(G->bestnull); free(G->mapnull);
	if(G->best != NULL) NilArchiveSites(G->best);
	if(G->map != NULL) NilArchiveSites(G->map);
	if(G->start != NULL) NilArchiveSites(G->start);
	if(G->prior != NULL) free(G->prior);
	if(G->expect != NULL) free(G->expect);
	free(G);
}

long	GetFreqProb(long t, long n, fm_type M, gs_type G, o_type R)
{
	st_type S=G->sites;
	char	*seq  = SeqSeqSet(n,G->data);
	long	i,end,pos,*order=G->order;
	double	denom,*freq_prob = PosProbSite(t,n,G->sites);

	end = SqLenSeqSet(n,G->data) - LenFModel(M) + 1;
	if(R!=NULL){
		OrderSites(n,order,S); pos = 0;
		denom = RelProbOrder(order,t,pos,R);
	} else denom = 1.0;
	for(freq_prob[0]=0.0, i= 1; i<= end; i++){
	    if(!OccupiedSite(t,n,i,S)){
		freq_prob[i] = (double)(denom * LikelihoodFModel(seq, i, M));
		freq_prob[0] += freq_prob[i];
	    } else {
		freq_prob[i] = 0.0;
	    	if(R != NULL && StartSite(n,i,S)) {
			pos++; /* == next site? */
			denom = RelProbOrder(order,t,pos,R);
		}
	    }
	}
}

long	*GetSiteFreq(gs_type G, long t, long d)
/*********************************************************************
   returns the residue frequencies for a site d residues right (or left
   if d negative) of sites of type t in S.  If position is blocked or 
   off ends of a sequence NULL is returned. 
   d = 0..len-1 means column is within the site.
 *********************************************************************/
{
	long	*site_freq,n,k,b,p,len;
	st_type	S=G->sites;
	ss_type	P=G->data;
	e_type	E;

	len = SiteLen(t,S);
	MEW(site_freq, nAlpha(SeqSetA(P))+2,long);
	for(b=0;b<= nAlpha(SeqSetA(P)); b++) site_freq[b] = 0;
	for(n=1; n <=NSeqsSeqSet(P); n++){
	   E = SeqSetE(n,P);
	   for(k=1; k<=nSites(t,n,S); k++){
		p = SitePos(t,n,k,S) + d;
		if(p > (long) SqLenSeqSet(n,P) || p < 1){
			free(site_freq); return (long*) NULL;
		}
		if((d >= len || d < 0) && !OpenPos(n,p,S)){
			free(site_freq); return (long*) NULL;
		}
		b = XnuSeq(p,E); site_freq[b]++;
		/** b = ResSeq(p,E); site_freq[b]++; /***/
	   }
	}
	return site_freq;
}

