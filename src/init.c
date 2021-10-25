/****
A Gibbs Sampler algorithm for finding multiple sites in multiple sequences 
****/
#include "gibbs.h"

fm_type	*InitMAPGibbs(gs_type G)
{
	long	i,b,n,k,s,t,ntyps,N,w,expect;
	double	npseudo,p,T;
	ss_type	P=G->data;
	a_type  A=SeqSetA(P);
	st_type S;
	bp_type *OldPriors;

  N = NSeqsSeqSet(P); ntyps = G->ntyps;
  /******************* Get MAP Sites **********************/
  if(G->sites != NULL)  NilSites(G->sites);
  G->sites=S=ExtractSites(G->map); 
  /******************* Get MAP priors **********************/
  if(G->prior != NULL){
	OldPriors = G->prior;
	NEW(G->prior,ntyps+2,bp_type);
	for(t = 1; t <= G->ntyps; t++){
           for(T=0.0,expect=0, n = 1; n <= N; n++) {
		expect += nSites(t,n,S);
                w = SqLenSeqSet(n,G->data) - SiteLen(t,S) + 1;
                if(w > 0) T += (double) w;
           }
	   T -= SiteLen(t,S)*expect; 
           G->prior[t] = MkBPrior(expect, G->weight, T);
           /***/ PutBPrior(stdout, G->prior[t]); /***TEST***/
	   NilBPrior(OldPriors[t]);
	   G->p[t]=G->map_p[t];
	}
	free(OldPriors);
  }
  /******************* Get MAP models **********************/
  if(G->model == NULL) NEW(G->model,ntyps+1,fm_type); 
  else for(t=1; t<= ntyps; t++){
	if(G->model[t]!= NULL) NilFModel(G->model[t]);
  }
  for(t=1; t<= ntyps; t++){
    if(G->prior!=NULL){
	npseudo=(double)ExpectBPrior(G->prior[t]);
    	npseudo = npseudo * G->pseudo;
    } else {
	for(npseudo=0.0,n=1; n<= N; n++)
	    	npseudo += (double) nSites(t,n,S);
	npseudo = npseudo * G->pseudo; 
    }
    G->model[t]=MkFModel(NULL,SiteLen(t,S),npseudo,CntsSeqSet(P),A);
    for(n =1; n <= LenFModel(G->model[t]); n++){
	   if(G->mapnull[t][n]) RmColumnFModel(n,G->model[t]);
    }
    for(n =1; n <= N; n++){
	PosTSites(t,n,G->pos,S);
	for(k = 1; k <= nSites(t,n,S); k++){
		Add2FModel(SeqSeqSet(n,P), G->pos[k],G->model[t]);
		if(G->prior != NULL) AddBPrior(G->prior[t]);
	}
    }
  }
  return G->model;
}

fm_type	*InitGibbs(gs_type G)
{
	long	w,i,b,n,k,s,t,ntyps,N;
	double	npseudo,p,T;
	ss_type	P=G->data;
	a_type  A=SeqSetA(P);
	st_type S;

     N = NSeqsSeqSet(P); ntyps = G->ntyps;
    /******************* Initialize priors **********************/
    if(G->prior != NULL){
	for(t = 1; t <= G->ntyps; t++){
	   ClearBPrior(G->prior[t]);
           for(T=0.0, n = 1; n <= N; n++) {
                w = SqLenSeqSet(n,G->data) - G->site_len[t] + 1;
                if(w > 0) T += (double) w;
           }
	   SetBPriorN(T,G->prior[t]);
	}
    }
    /******************* Initialize sites **********************/
    if(G->sites != NULL)  NilSites(G->sites);
    if(G->prior != NULL){
	S = MkSites(ntyps, G->site_len, P); G->sites = S;
        for(t=1; t<= ntyps; t++){
	    if(!NRandomSites(t, G->expect[t], 500, S)){
		PutSites(stderr, t, S, NULL,NULL);
		fprintf(stderr,"%d expected\n",G->expect[t]);
		print_error("too many sites; try a smaller number");
	    }
	}
    } else {
	S = G->sites = ExtractSites(G->start);
	if(!ShuffleSites(G->sites)){
		print_error("too many sites; try a smaller number");
	}
    }
    /******************* Initialize models **********************/
  if(G->model == NULL) NEW(G->model,ntyps+1,fm_type); 
  else for(t=1; t<= ntyps; t++){
	if(G->model[t]!= NULL) NilFModel(G->model[t]);
  }
  for(t=1; t<= ntyps; t++){
    if(G->prior!=NULL) {
	npseudo=(double)G->expect[t];
    	npseudo = npseudo * G->pseudo;
    } else {
	for(npseudo=0.0,n=1; n<= NSeqsSeqSet(P); n++)
	    npseudo += (double) nSites(t,n,S);
    	npseudo = npseudo * G->pseudo;
    }
    G->model[t]=MkFModel(NULL,G->site_len[t],npseudo,CntsSeqSet(P),A);
    for(n =1; n <= N; n++){
	PosTSites(t,n,G->pos,S);
	for(k = 1; k <= nSites(t,n,S); k++){
		Add2FModel(SeqSeqSet(n,P), G->pos[k],G->model[t]);
		if(G->prior != NULL) AddBPrior(G->prior[t]);
	}
    }
    if(G->prior != NULL) G->p[t]=PostProbBPrior(G->prior[t]);
  }
  return G->model;
}

