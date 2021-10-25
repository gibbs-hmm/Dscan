/****
A Gibbs Sampler algorithm for finding multiple sites in multiple sequences 
****/
#include "gibbs.h"

double  NetMapBGibbs1(long t, fm_type M, gs_type G)
{
	long	N,n,w;
	double	m1,m2;

        for(N=0,n=1; n <= NSeqsSeqSet(G->data); n++){
                w = SqLenSeqSet(n,G->data) - LenFModel(M) + 1;
                if(w > 0) N += w;
        }
	SetBPriorN(N, G->prior[t]);
	m1 = LnMapFModel(M) - LnMapFModelNull(M);
	m2 = BPriorMAP(G->prior[t]) - BPriorNullMAP(G->prior[t]);
	fprintf(stdout,"model map = %g; betaprior map = %g\n", m1,m2);
	return (m1 + m2);
}

double  MapBGibbs(gs_type G)
{
    long	N,n,w,t;
    double	MAP=0.0;

    for(t = 1; t <= nTypeSites(G->sites); t++) {
        for(N=0,n=1; n <= NSeqsSeqSet(G->data); n++){
                w = SqLenSeqSet(n,G->data) - LenFModel(G->model[t]) + 1;
                if(w > 0) N += w;
        }
	SetBPriorN(N, G->prior[t]);
        MAP += BPriorMAP(G->prior[t]);
        MAP += LnMapFModel(G->model[t]);
    }
    return MAP;
}

double  NullMapBGibbs(gs_type G)
{
    long	N,n,w,t;
    double	NullMap=0.0;

    for(t = 1; t <= nTypeSites(G->sites); t++) {
        for(N=0,n=1; n <= NSeqsSeqSet(G->data); n++){
                w = SqLenSeqSet(n,G->data) - LenFModel(G->model[t]) + 1;
                if(w > 0) N += w;
        }
	SetBPriorN(N, G->prior[t]);
        NullMap += BPriorNullMAP(G->prior[t]);
        NullMap += LnMapFModelNull(G->model[t]);
    }
    return NullMap;
}

double  NetMapBGibbs(gs_type G){return (MapBGibbs(G)-NullMapBGibbs(G));}

long	SaveBestGibbs(gs_type G)
{
	long	t,n;

	if(G->best != NULL) NilArchiveSites(G->best);
	G->best = ArchiveSites(G->sites);
	for(t = 1; t <= nTypeSites(G->sites); t++){
	   NullSitesFModel(G->bestnull[t],G->model[t]);
	   if(G->best_p != NULL){	/** i.e., is motif sampler **/
	      G->best_p[t] = G->p[t];
	   }
	}
}

Boolean	SaveFinalGibbs(gs_type G)
{
	long	t,n,k,i,*tmpleng,***tmpsites,**tmpnsite;
	Boolean	converged=TRUE,**tmpnull;

	if(G->map !=NULL) {
		converged = SameArchiveSites(G->map, G->best);
		NilArchiveSites(G->map);
	} else converged = FALSE;
	G->map = G->best; G->best = NULL;
	if(G->best_p != NULL){
	   for(t=1; t<=G->ntyps; t++){ G->map_p[t]=G->best_p[t]; }
	}
	tmpnull = G->mapnull; G->mapnull = G->bestnull;
	G->bestnull = tmpnull;
	return converged;
}


