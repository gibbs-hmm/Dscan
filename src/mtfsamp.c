#include "gibbs.h"

long	MotifSampler(FILE *fptr, gs_type G)
{
	long	nruns=G->nruns, nread=G->nread;
	long	run,t,n,s;
	st_type	S;
	double	***prob,best,oldbest;
	char	str[200];
	fm_type	M;
	Boolean	converged,okay;
	FILE    *tfptr;

    fprintf(stderr,"nconv = %d\n", G->nconverge);
    oldbest = -DBL_MAX; converged = FALSE;
    s=0;
    for(run = 1; run <= nruns && !converged; run++){
	fprintf(stderr,"\r** %d **\n", run);
	if(oldbest < (best=motif_sampler(G))){ 
		s=run;
		oldbest = best;
		converged = SaveFinalGibbs(G); 
		for(t = 1; t <= nTypeSites(G->sites); t++){ 
			fprintf(stderr,"\tMAP[%d] = %g\n",t,best);
		}
	} 
        fprintf(stderr,"run %d MAP = %g\n",run,best);
	
    }
    fprintf(stderr,"best run = %d\n",s);
    InitMAPGibbs(G);
  if(!G->test[2] && !G->test[3]){
    okay = FALSE;
    for(t = 1; t <= nTypeSites(G->sites); t++) {
     sprintf(str,"MOTIF %c",t + 'A' -1);
     fprintf(fptr,"\n\n%*s%*s%s",
		23,"",(SiteLen(t,G->sites)-7)/2,"",str);
     if(TotalSites(t,G->sites) > 0) {
        if(G->fragment) {
                NullSitesFModel(G->null, G->model[t]);
                PutSites(fptr,t,G->sites,NULL,G->null);
        } else { PutSites(fptr,t,G->sites,NULL,NULL); }
	PutFModel(fptr,G->model[t]);
        okay = TRUE;
     } else fprintf(fptr," (failed: no sites)\n\n");
    }
    if(!okay) return;
    fprintf(fptr,"\n NetMAP = %f\n\n", NetMapBGibbs(G));
  }
  if(nread > 0){ 		/*** NEW ***/
	NearOptFit(fptr, G);
        if(G->wilcoxon > 0){
    	   for(t=1; t<=nTypeSites(G->sites); t++) {
		if(G->name != NULL) {
                    sprintf(str,"%s%d",G->name,t);
		} else sprintf(str,"temp_file%d",t);
                tfptr = open_file(str,".wc","w");
                PutWilcoxBGibbs(tfptr,1, G->wilcoxon, t, G);
		PutWilcoxBGibbs(tfptr,G->wilcoxon+1,
			NSeqsSeqSet(G->data),t,G);
                fclose(tfptr);
           }
        }
	fprintf(fptr,"\tseed: %d\n", G->seed);
  }
}

long	NearOptFit(FILE *fptr, gs_type G)
/** fit the sites into the sequences - most probably first **/
{
	long	nruns=G->nruns, nread=G->nread;
	long	run,t,n,s;
	st_type	S;
	double	***prob,best,oldbest;
	keytyp	key;
	char	str[50];
	long	*item_s,*item_n,*item_t,item;
	fm_type	*M;
	Boolean	converged;
	dh_type	H;
	double	*missing_info,p,q,minfo,binfo,L0,df;

    L0 = LogL0SeqSet(G->data);
    S=CopySites(G->sites);
    prob = NearOptimumSampler(nread, G);
    n = TotalSeqSet(G->data);
    NEW(M,nTypeSites(G->sites)+2,fm_type);
    NEW(item_t,n+3,long);
    NEW(item_n,n+3,long);
    NEW(item_s,n+3,long);
    NEW(missing_info,nTypeSites(G->sites)+2,double);
    H = dheap(n,3);
    for(item=1,t = 1; t <= nTypeSites(G->sites); t++){ 
	ClearBPrior(G->prior[t]);
	M[t] = CopyFModel(G->model[t]); 
	for(n=1; n <= NSeqsSeqSet(G->data); n++){
	    for(s=1; s<= (long) SqLenSeqSet(n,G->data); s++){
		if(prob[t][n][s] >= G->readcutoff){
			key = (keytyp) prob[t][n][s];
			insrtHeap(item,key,H);	
			item_t[item]=t;
			item_n[item]=n;
			item_s[item]=s;
			item++;
	        }
	    }
	}
    }
    for( ; (item=delminHeap(H)) != NULL; ){
	t = item_t[item];
	n = item_n[item];
	s = item_s[item];
	if(!OccupiedSite(t,n,s,S)) {
		AddSite(t,n,s,S);
		Add2FModel(SeqSeqSet(n,G->data),s,M[t]);
		AddBPrior(G->prior[t]);
	}
    }
    for(item=1,t = 1; t <= nTypeSites(G->sites); t++){ 
	for(n=1; n <= NSeqsSeqSet(G->data); n++){
	    for(s=1; s<= (long) SqLenSeqSet(n,G->data); s++){
		if((p=prob[t][n][s]) > 0.0){
			missing_info[t] += p*log(p);
		}
		if((p=prob[t][n][s]) < 1.0){
			missing_info[t] += (1.0-p)*log(1.0-p);
		}
	    }
	}
    }
    if(G->mfptr!=NULL){
        fprintf(G->mfptr,"//\nID   XXX\nAC   A00000\nDE   %s\n", G->name);
        fprintf(G->mfptr,"CC   comments\nNU   %d\n", nTypeSites(S));
    }
    for(t = 1; t <= nTypeSites(G->sites); t++){ 
	df = (double) nAlpha(G->A) - 1.0;
	df *= nColsFModel(M[t]);
	if(G->fragment) {
          NullSitesFModel(G->null, M[t]);
          PutSites(fptr,t,S,prob[t],G->null);
	  if(G->mfptr!=NULL) PutSitesMtfDBS(G->mfptr,t,S,prob[t],G->null);
	  if(G->sfptr!=NULL) PutScanSites(G->sfptr,t,S,G->null);
	} else {
	  PutSites(fptr,t,S,prob[t],NULL); 
	  if(G->mfptr!=NULL) PutSitesMtfDBS(G->mfptr,t,S,prob[t],NULL);
	  if(G->sfptr!=NULL) PutScanSites(G->sfptr,t,S,NULL);
	}
    	PutFModel(fptr,M[t]); 
    	if(G->test[4]) PutLogFModel(fptr,M[t]); 
	/** PutAlphaR(fptr,G->A); /***/
	fprintf(fptr,"\n MAP = %f\n\n",NetMapBGibbs1(t,M[t],G));
        if(G->verbose) {
	   minfo = LogLikeFModel(M[t]) - L0;
	   fprintf(fptr,"Model: %g bits\n", minfo);
	   binfo = LogLikeBPrior(G->prior[t]);
	   fprintf(fptr,"betaprior: %g bits\n", binfo);
	   missing_info[t] *= 1.4427;
	   fprintf(fptr,"Missing info: %g bits\n", missing_info[t]);
	   fprintf(fptr,"Total Info: %g bits\n", 
		minfo + binfo - missing_info[t]);
	   fprintf(fptr,"Information per parameter: %g bits\n", 
		(minfo + binfo - missing_info[t])/df);
	   fprintf(fptr,"AIC: %g bits\n", 
		2*(minfo + binfo - missing_info[t]) - (2.0*df*1.442695));
        }
	NilFModel(M[t]);
    }
    fprintf(fptr,"\n\n");
    NilSites(S); Nildheap(H);
    free(item_t); free(item_n); free(item_s);
    free(missing_info);
    free(M);
}

double	motif_sampler(gs_type G)
{
	long	w,i,b,n,s,s2,iter,k,t,N,ncycles=G->ncycles,end;
	long	number,ntyps,*endt,fcols;
	fm_type	*model;
	st_type	sites;
	ss_type  data = G->data;
	Boolean	quit=FALSE,lowtemp=FALSE;
	double	best,tprob,*prob,totlike=0.0,ratio,rand,nullmap;
	double	***rprob=G->readprob,*sumprob=G->sumprob;
	h_type	*H;

	InitGibbs(G);
	model = G->model; sites = G->sites; 
	N=NSeqsSeqSet(data);
	best = -DBL_MAX; ntyps = nTypeSites(sites);
	NEW(prob,ntyps+2,double); NEW(endt,ntyps+2,long);

	for(fcols=99999,t=1; t<=ntyps; t++) {
		fcols = MIN(long,fcols,nColsFModel(model[t]));
	}
	nullmap = NullMapBGibbs(G);
	fcols = MAX(long,fcols/3,2);
	for(number=0,iter=1; iter <= G->nconverge && !quit; iter++) {
	   if(!lowtemp && iter > 3  && iter % ncycles == 0) {
	      if(G->move && ntyps > 1 && G->fragment
				&& TransferColumn(G, ntyps, model)){
		nullmap = NullMapBGibbs(G);
		fprintf(stderr,"cycle %d\n", iter);
	      }
	      for(t = 1; t <= ntyps; t++) {
		if((G->fragment && MoveMultiCols(G,t,fcols))
		   || ContigFModel(model[t]) && Metropolis(G,t,model[t])){
			nullmap = NullMapBGibbs(G);
			fprintf(stderr,
			   " motif %c cycle %d AP %.1f (%d sites)\n",
			      t+'A'-1,iter,totlike,TotalSites(t,sites));
		}
	      } if(iter % 5 == 0) fprintf(stderr,"\r%d",iter);
	   }
	   for(n =1; n <= N && !quit; n++) {
		PosSites(n,G->pos,sites);
		for(w=end=0,t=1; t<=ntyps; t++) { 
		        if(w < SiteLen(t,sites)) w=SiteLen(t,sites);
		        endt[t] = SqLenSeqSet(n,data) - SiteLen(t,sites) + 1;
			end = MAX(long,end,endt[t]);
		}
		for(k=s=1,s2=G->pos[k]; s<=end; s++){
		   while(s2 > 0 && (s >= (s2 - w))){ 
			t=TypeSite(n,s2,sites);
			VacateSite(t,n,s2,sites);
			RmFModel(SeqSeqSet(n,data), s2, model[t]);
			k++; s2=G->pos[k];
			RmBPrior(G->prior[t]);
		   }
		   for(t=1; t<=ntyps; t++) {
		      if(s<=endt[t]){
		        if(iter > 3) G->p[t]= PostProbBPrior(G->prior[t]);
			prob[t]=ProbFModel(SeqSeqSet(n,data),s,
						G->p[t],model[t]);
		      } else prob[t]=0.0;
		   }
		   rand = (double) Random()/(double) LONG_MAX;
		   for(tprob=0.0, t = 1; t <= ntyps; t++) {
			tprob += prob[t];
		        if(prob[t] > 0.0 && rand <= tprob){
		 	    AddSite(t,n,s,sites);
			    Add2FModel(SeqSeqSet(n,data),s,model[t]);
			    AddBPrior(G->prior[t]);
		            s += SiteLen(t,sites) - 1;
			    break;
			}
		   }
		}
		if(iter > 3) {
        	 totlike=MapBGibbs(G) - nullmap;
		 if(totlike<=best) number++;
		 else {			
                   number=0;best=totlike;
		   for(t = 1; t <= ntyps; t++){
			G->p[t]= PostProbBPrior(G->prior[t]);
		   }
		   SaveBestGibbs(G);
		 }
		 if(number > G->stop) { /** stop fragmenting && quit **/
		     if(lowtemp) { quit=TRUE; break; }
		     else { lowtemp = TRUE; number=0; }
		 }
		} 
	    }
	}
	free(prob); free(endt); 
	return best;
}


double	***NearOptimumSampler(long niter, gs_type G)
/********************************************************************
 Returns probabilities p[t][n][s] for converged model; keep G->p[t] 
 constant during the run. - limits the search.
********************************************************************/
{
	long	w,i,b,n,s,s2,iter,k,t,N,end,ntyps,*endt;
	fm_type	*model;
	st_type	sites;
	ss_type  data = G->data;
	double	***rprob=G->readprob,tprob,*prob,rand,ap,map,tmp,p,cut;
	Boolean	read=FALSE, ***good;
	long	ngood,nmodel,nbad;
	h_type	H;

	model = G->model; sites = G->sites;
	N=NSeqsSeqSet(data); ntyps = nTypeSites(sites);
	NEW(prob,ntyps+2,double); NEW(endt,ntyps+2,long);
	for(t = 1; t <= ntyps; t++){
	   for(n = 1; n <= N; n++) {
		for(s=0;s<=(long)SqLenSeqSet(n,data);s++)rprob[t][n][s]=0.0;
	   }
	}
	cut = 0.0;
	cut = 0.01; 
    	map=NetMapBGibbs(G);
/**** find reasonable possibitities ****/
	NEWPP(good,ntyps+2,Boolean);
	for(t=0; t<=ntyps; t++) {
	   NEWP(good[t],N+2,Boolean);
	   for(n=1;n<=N;n++) NEW(good[t][n],SqLenSeqSet(n,data)+2,Boolean);
	}
	ngood=nmodel=nbad=0;
	for(n = 1; n <= N; n++) {
		PosSites(n,G->pos,sites);
		for(w=end=0,t=1; t<=ntyps; t++) { 
		        if(w < SiteLen(t,sites)) w=SiteLen(t,sites);
		        endt[t] = SqLenSeqSet(n,data) - SiteLen(t,sites) + 1;
			end = MAX(long,end,endt[t]);
		}
		/** get site locations & scan through sequences **/
		for(k=s=1,s2=G->pos[k]; s<=end; s++){
		   if(s==s2){	/** if there is a site here **/
		        t=TypeSite(n,s,sites);
			good[0][n][s] = TRUE;
			good[t][n][s] = TRUE;
			nmodel++;
			k++; s2=G->pos[k]; 
		   } else {
		   	for(t=1; t<=ntyps; t++) {
		          if(s<=endt[t]){
			    p=ProbFModel(SeqSeqSet(n,data),s,G->p[t],model[t]);
			    if(p > cut){
				 good[0][n][s] = TRUE;
				 good[t][n][s] = TRUE;
			    } 
		          } 
		   	}
			if(good[0][n][s]) ngood++;
			else nbad++;
		   }
		}
	}
	if(!G->test[2]) printf("sites: MAP = %d; maybe = %d; discard = %d\n",
			nmodel,ngood,nbad);
/**** end: find reasonable possibitities ****/
    	if(!G->test[2]) fprintf(stdout," Initial MAP = %f\n", map);
	for(iter =1; iter <= niter; iter++) {
	  if(iter % 5 == 0) fprintf(stderr,"\r%d",iter);
	  for(n = 1; n <= N; n++) {
		PosSites(n,G->pos,sites);
		for(w=end=0,t=1; t<=ntyps; t++) { 
		        if(w < SiteLen(t,sites)) w=SiteLen(t,sites);
		        endt[t] = SqLenSeqSet(n,data) - SiteLen(t,sites) + 1;
			end = MAX(long,end,endt[t]);
		}
		for(k=s=1,s2=G->pos[k]; s<=end; s++){
                 if(good[0][n][s]){
		   if(s==s2){
		        t=TypeSite(n,s2,sites);
			VacateSite(t,n,s2,sites);
			RmFModel(SeqSeqSet(n,data), s2, model[t]);
			RmBPrior(G->prior[t]);
			k++; s2=G->pos[k];
		   }
		   for(t=1; t<=ntyps; t++) {
		      if(s<=endt[t] && good[t][n][s] 
				&& !OccupiedSite(t,n,s,sites)){
			prob[t]=ProbFModel(SeqSeqSet(n,data),s,G->p[t],model[t]);
		      } else prob[t]=0.0;
		   }
		   rand = (double) Random()/(double) LONG_MAX;
		   for(tprob=0.0, t = 1; t <= ntyps; t++) {
			tprob += prob[t];
		        if(prob[t] > 0.0 && rand <= tprob){
		 	    AddSite(t,n,s,sites);
			    Add2FModel(SeqSeqSet(n,data),s,model[t]);
			    AddBPrior(G->prior[t]);
			    if(read) rprob[t][n][s] += 1.0;
		            s += SiteLen(t,sites) - 1;
			    break;
			}
		   }
		 }
		}
	   }
	   ap = NetMapBGibbs(G);
	   if(ap > map) { 
		map=ap; 
		if(!read) {
			fprintf(stderr," MAP = %f\n",map);
			iter=0; 
		}
	   }
	   if(!read && iter > 30) { iter=0; read = TRUE; }
	}
	for(t=1; t<=ntyps; t++) {
	   if(G->verbose) H = Histogram("site probabilities", 0, 1, 0.05);
	   for(n=1; n<=N; n++){
	     for(s=1;s<=(long)SqLenSeqSet(n,data); s++) {
		rprob[t][n][s]/=(double)niter;
		if(G->verbose && rprob[t][n][s] >= 0.01) {
			IncdHist(rprob[t][n][s],H);
		}
	     }
	     free(good[t][n]);
	   }
	   free(good[t]);
	   if(G->verbose) { PutHist(stdout,60,H); NilHist(H); }
	}
	for(n=1; n<=N; n++) free(good[0][n]); free(good[0]);
	free(good);
	free(endt); free(prob);
	if(!G->test[2]) fprintf(stdout," MAP = %f\n",map);
	return rprob;
}

