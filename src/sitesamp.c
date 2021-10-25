#include "gibbs.h"

long    RunGibbs(FILE *fptr, gs_type G)
{
        long time1;

        time1=time(NULL);
        if(G->prior != NULL) MotifSampler(fptr,G);
        else SiteSampler(fptr, G);
        if(!G->test[2])
          fprintf(stdout,"\ttime: %ld seconds (%0.2f minutes)\n",
                time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
}

/*************************************************************************
Gibbs Samplers for finding multiple sites in multiple sequences.
**************************************************************************/
void	SiteSampler(FILE *fptr, gs_type G)
{
	long	nruns=G->nruns;
	long	b,n,s,iter,k,newsite,t,N=NSeqsSeqSet(G->data),npos,hpsz=100;
	long	*df,site, ncycles=G->ncycles,end;
	long	nconverg = G->nconverge,ntyps = G->ntyps;
	long	*order=G->order,run,number,nipp;
	FILE	*ifptr=G->ifptr;
	fm_type	*model,*finalmodel=NULL;
	st_type	mapsites=NULL,sites;
	o_type	R=NULL;
	ss_type  data = G->data;
	a_type	A=SeqSetA(G->data);
	Boolean	flag;
	double	ipp,L0,*Like[3],*MissInfo[3],TotLike,*dTemp;
	double	pseudo = G->pseudo,qseudo=G->qseudo,totipp;
	double	best,oldbest,ratio,*pos_prob;
	char	str[50],*seq;
	mh_type	H;
	keytyp	key;
/*** NEW ****/
	double	best_prob;
	long	*typ_order,i;
/*** NEW ****/

   fprintf(stderr,"nconv = %d\n", nconverg);
   H = Mheap(hpsz,3); L0 = LogL0SeqSet(data); 
   for(t=0; t<= 2; t++) {
	NEW(Like[t],ntyps+1, double);
	NEW(MissInfo[t], ntyps+1, double);
   }
   NEW(df, ntyps+1, long);
   oldbest = 1.0; best = 0.0; nipp=0;
   for(flag = TRUE, run = 1; run <= nruns && oldbest != best; run++){
	if(oldbest < best) oldbest = best;
	best = 0.0;
	model = InitGibbs(G);
	sites = G->sites;
	if(run == 1){
	   for(npos=0,t=1; t<= ntyps; t++) {
	       	df[t] = (G->site_len[t]*(nAlpha(A) - 1));
		k=nSites(t,1,sites); npos += k;
		for(n = 2; n <= N; n++) 
			if(k!=nSites(t,n,sites)) flag = FALSE;
   	   }
	}
	if(ntyps > 1 && flag && G->use_order) {
		if(R != NULL) { PutOrder(stderr,R); NilOrder(R); }
		R=Order(ntyps,npos,qseudo*(double)N);
		InitOrder(R);
		for(n = 1; n <= NSeqsSeqSet(data); n++) {
			OrderSites(n,order,sites);
			Add2Order(order,R);
		}
	} else R = NULL;
	fprintf(stderr,"\r** %d **\n", run);
	for(number=0,iter =1; iter <= nconverg; iter++) {
	  if(iter > 2 && iter % ncycles == 0) {
            if(G->move && ntyps>1 && TransferColumn(G,ntyps,model)){
                fprintf(stderr,"cycle %d\n", iter);
            }
	    for(t = 1; t <= nTypeSites(sites); t++) {
		if(G->fragment && MoveColumn(G,t)){ 
		   fprintf(stderr," motif %c cycle %d\n",t+'A'-1,iter);
		} 
		if(ContigFModel(model[t]) && Metropolis(G,t,model[t])){
			fprintf(stderr," motif %c cycle %d\n",
				t+'A'-1,iter);
		}
		if(iter % 5 == 0) fprintf(stderr,"\r%d",iter);
	    }
	  }
	  for(n = 1; n <= N; n++) {
	      for(t = 1; t <= nTypeSites(sites); t++) {
		PosTSites(t,n,G->pos,sites);
		end = nSites(t,n,sites);
		for(k = 1; k <= end; k++){
		    if(R!=NULL) {
			OrderSites(n, order,sites);
			RmOrder(order,R);
		    }
		    s = G->pos[k];
		    VacateSite(t,n,s,sites);
		    RmFModel(SeqSeqSet(n,data), s, model[t]);
		    GetFreqProb(t,n,model[t],G,R);
		    s = ChooseSite(t,n,sites);
		    Add2FModel(SeqSeqSet(n,data),s,model[t]);
		    if(R!=NULL) {
			OrderSites(n,order,sites);
			Add2Order(order,R);
		    }
		}
	      }
	  }
	  for(TotLike=totipp=0.0, t = 1; t <= nTypeSites(sites); t++) {
       	 	MissInfo[0][t] = MissInfoSites(t,sites);
        	Like[0][t] = LogLikeFModel(model[t]);
		TotLike += Like[0][t] - L0;
		ipp = (Like[0][t] - MissInfo[0][t] - L0)/(double) df[t];
		totipp += ipp;
	  }
	  if(ifptr != NULL) { fprintf(ifptr,"%1.3g\n",totipp); nipp++; }
	  key = (keytyp) -TotLike;
	  InsertMheap(key, H);
	  if(TotLike > best){ 
		number = 0; best = TotLike; 
		SaveBestGibbs(G);
		dTemp = Like[1]; Like[1] = Like[0]; Like[0]= dTemp;
		dTemp = MissInfo[1]; MissInfo[1] = MissInfo[0]; 
		MissInfo[0]= dTemp;
	  } else if(++number > G->limit) break; 
	}
	if(oldbest < best) {
	   dTemp = Like[2]; Like[2] = Like[1]; Like[1]= dTemp;
	   dTemp = MissInfo[2]; MissInfo[2] = MissInfo[1]; 
	   MissInfo[1]= dTemp;
	   SaveFinalGibbs(G);
	} 
	for(t = 1; t <= ntyps; t++){ NilFModel(model[t]); }
	free(model); G->model = NULL; 
	NilSites(sites); G->sites = NULL; 
	if(oldbest == best) break;
    }
    fprintf(stderr,"best = %g; oldbest = %g\n",best,oldbest);
    finalmodel = InitMAPGibbs(G);
    mapsites = G->sites;
    if(R != NULL) {
	NilOrder(R);
	R=Order(nTypeSites(mapsites),npos,qseudo*(double)N);
	InitOrder(R);
	for(n = 1; n <= NSeqsSeqSet(data); n++) {
		OrderSites(n,order,mapsites);
		Add2Order(order,R);
	}
    }
    typ_order = NULL;
    if(G->mfptr!=NULL){
	if(R != NULL) typ_order = ConcensusOrder(R);
        fprintf(G->mfptr,"//\nID   XXX\nAC   A00000\nDE   %s\n", G->name);
        fprintf(G->mfptr,"CC   comments\nNU   %d\n", nTypeSites(mapsites));
    }
    if(typ_order == NULL) {
	  NEW(typ_order,nTypeSites(mapsites)+3,long);
	  for(t=1; t<=nTypeSites(mapsites); t++) typ_order[t] = t;
    }
    for(i=1; i <= nTypeSites(mapsites); i++) {
	t = typ_order[i];
	sprintf(str,"MOTIF %c",t + 'A' -1);
	fprintf(fptr,"\n\n%*s%*s%s",
		23,"",(SiteLen(t,mapsites)-7)/2,"",str);
	if(G->fragment) {
	  NullSitesFModel(G->null, finalmodel[t]);
	  PutSites(fptr,t,mapsites,NULL,G->null);
/****
	  if(G->mfptr!=NULL) PutSitesMtfDBS(G->mfptr,t,mapsites,NULL,G->null);
	  if(G->sfptr!=NULL) PutScanSites(G->sfptr,t,mapsites,G->null);
/****/
	} else { 
	  PutSites(fptr,t,mapsites,NULL,NULL); 
/****
	  if(G->mfptr!=NULL) PutSitesMtfDBS(G->mfptr,t,mapsites,NULL,NULL);
	  if(G->sfptr!=NULL) PutScanSites(G->sfptr,t,mapsites,NULL);
/****/
	}
	PutFModel(fptr,finalmodel[t]);
	fprintf(fptr,
		"\n\nComplete log-likelihood ratio  = %4d bits\n",
                       (long) (Like[2][t] - L0));
	fprintf(fptr,"Missing position information   = %4d bits\n",
                       (long) MissInfo[2][t]);
       	fprintf(fptr,"Log-likelihood ratio statistic = %4d bits\n",
                       (long)(Like[2][t] -  MissInfo[2][t] - L0));
       	fprintf(fptr,"Degrees of freedom             = %4d\n",df[t]);
       	fprintf(fptr,"Information per parameter      = %1.3g bits\n\n",
                       ((Like[2][t] -  MissInfo[2][t] - L0)/(double)df[t]));
/****** NEW: compute probabilites *****/
	for(best_prob=-DBL_MAX,n = 1; n <= N; n++) {
		end = SqLenSeqSet(n,G->data) - LenFModel(finalmodel[t]) + 1;
		pos_prob = PosProbSite(t,n,mapsites);
		seq = SeqSeqSet(n,G->data);
		for(pos_prob[0]=0.0, s= 1; s<= end; s++){
		   if(TypeSite(n,s,mapsites) == t){
		        RmFModel(seq, s, finalmodel[t]);
		   	pos_prob[s] = (double)
			   LikelihoodFModel(seq, s, finalmodel[t]);
		        Add2FModel(seq,s,finalmodel[t]);
		   } else pos_prob[s] = (double)
			   LikelihoodFModel(seq, s, finalmodel[t]);
		   pos_prob[0] += pos_prob[s];
		   best_prob = MAX(double,best_prob,pos_prob[s]);
		}
	}
	best_prob = log(best_prob);
	for(n = 1; n <= N; n++) {
		end = SqLenSeqSet(n,G->data) - LenFModel(finalmodel[t]) + 1;
		pos_prob = PosProbSite(t,n,mapsites);
		for(s = 1; s <= end; s++) {
			pos_prob[s] = log(pos_prob[s]);
			pos_prob[s] /= best_prob;
		}
	}
/***** OMIT *****
	for(n = 1; n <= N; n++) {
		pos_prob = PosProbSite(t,n,mapsites);
		end = SqLenSeqSet(n,G->data) - SiteLen(t,mapsites) + 1;
		for(s = 1; s <= end; s++){
		    if(TypeSite(n,s,mapsites) == t){
			if(pos_prob[s] < 0.25){
		    	   VacateSite(t,n,s,mapsites);
		    	   RmFModel(SeqSeqSet(n,data), s, finalmodel[t]);
			}
		    } else if(!OccupiedSite(t,n,s,mapsites)) {
			if(pos_prob[s] >= 0.25){
			   AddSite(t, n, s, mapsites);
		    	   Add2FModel(SeqSeqSet(n,data),s,finalmodel[t]);
			}
		    }
		}
	}
/***** OMIT *****/
	sprintf(str,"MOTIF %c",t + 'A' -1);
	fprintf(fptr,"\n\n%*s%*s%s",
		23,"",(SiteLen(t,mapsites)-7)/2,"",str);
	if(G->fragment) {
	  NullSitesFModel(G->null, finalmodel[t]);
	  PutSites(fptr,t,mapsites,ProbSite(t,mapsites),G->null);
	  if(G->mfptr!=NULL) PutSitesMtfDBS(G->mfptr,t,mapsites,
				ProbSite(t,mapsites),G->null);
	  if(G->sfptr!=NULL) PutScanSites(G->sfptr,t,mapsites,G->null);
	} else { 
	  PutSites(fptr,t,mapsites,ProbSite(t,mapsites),NULL); 
	  if(G->mfptr!=NULL) PutSitesMtfDBS(G->mfptr,t,mapsites,
				ProbSite(t,mapsites),G->null);
	  if(G->sfptr!=NULL) PutScanSites(G->sfptr,t,mapsites,NULL);
	}
	PutFModel(fptr,finalmodel[t]);
/****** NEW *****/
    }
    free(typ_order);
    fprintf(fptr,"\n    seed %ld",G->seed);
/****
    fprintf(stderr,"highest likelihood ratios:\n");
    for(t=0; !EmptyMheap(H); t++){
	fprintf(stderr,"%5.0f", -MinKeyMheap(H));
	if(t%10 == 9) fprintf(stderr,"\n");
	DelMinMheap(H);
    }
    fprintf(stderr,"\n\n");
/****/
    if(ifptr != NULL) { 
		fprintf(stderr,"\tnumber of ipp entries = %d\n",nipp);
    }
    if(R != NULL) { 
	PutOrder(fptr,R); 
	if(G->verbose) PutTypeSites(fptr, mapsites);
	NilOrder(R); 
    }
    NilMheap(H); free(df);
    for(t=0; t<= 2; t++) { free(Like[t]); free(MissInfo[t]); }
    if(ifptr != NULL) fclose(ifptr);
}

