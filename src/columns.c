#include "gibbs.h"

Boolean	MoveMultiCols(gs_type G, long t, long num)
{
	long	i,lemon;
	Boolean	result=FALSE;

	ClearMarksFModel(G->model[t]);
	for(i=1; i<=num; i++){
		lemon = ChoiceLemonFModel(G->model[t]); /****/
		/*** lemon = ChoiceOrangeFModel(G->model[t]); /****/
		if(move_column_gibbs(G, lemon, t)) {
			result=TRUE; fprintf(stderr,"\n");
		} else break;
	}
	return result;
}

Boolean	MoveColumn(gs_type G, long t)
{
	long	lemon;
	
	lemon = LemonFModel(G->model[t]); /** proportional to prob **/
	/*** lemon = OrangeFModel(G->model[t]); /*** randomly ***/
	return move_column_gibbs(G, lemon, t);
}

double	mv_col_weight(long c, long w, long n, long lemon, Boolean *null)
/********************************************************************

   If a fmodel has c on columns and a length of w.  Then there are 
   w-2 choose c-2  configurations for that model.

          n     1 lemon       w_old
          |     |  |           |
        ........*??*???????????*.............
                 |--- (w-2) --|
               .....
              -10234
      n = -inf..0 or +2..+inf.	 (i.e., n != 1)

       lemon = 1..w and n = -f..end.

 ********************************************************************/
{
	long start=1, end=w, w_new,w_old;
	double	b0,b2,P;
	long	i,s,e;

	w_old = w;
	if(lemon == 1) {     /** case 1: **/
	   if(n < 1){          /** case 1a: increase width by |n| **/
		w_new = w_old - n + 1;  /** note: n is <= 0 **/
	   } else if(n > 1){   /** case 1b: **/
		for(start=2; null[start]; start++) ;
		if(n < w_old){
		    start = MIN(long,start,n);
		    w_new = w_old - start + 1;
		} else if(n > w_old){
		    w_new = n - start + 1;   
		} else print_error("input error in mv_col_weight( )");
	   } else print_error("input error in mv_col_weight( )");
	} else if(lemon == w_old){	/** case 2: **/
	   for(end = w_old-1; null[end]; end--) ;
	   if(n < 1){
		w_new = end - n + 1;  /** note: n is <= 0 **/
	   } else if(n > 1){
		if(n < end) w_new = end;
		else if(n > end) w_new = n; 
		else print_error("input error in mv_col_weight( )");
	   } else print_error("input error in mv_col_weight( )");
	} else { 		/** case 3: 1 < lemon < w_old **/
	   if(n < 1){
		w_new = w_old - n + 1;  /** note: n is <= 0 **/
	   } else if(n > 1){
		if(n < w_old) w_new = w_old; 
		else if(n > w_old) w_new = n;
		else print_error("input error in mv_col_weight( )");
	   } else print_error("input error in mv_col_weight( )");
	}
/****
	b0 = bico(w_old-2,c-2);
	b2 = bico(w_new-2,c-2);
	P = b0/b2;
if(P > 1000000){
	fprintf(stderr,"w_old: ");
	s=MIN(long,n,start);
	e=MAX(long,n,w_old);
	for(i=s; i<=w_old; i++){
	   if(i >= 1 && i <= w_old){
		if(null[i]) fprintf(stderr,".");
		else fprintf(stderr,"*");
	   }else fprintf(stderr," ");
	}
	fprintf(stderr,"\nw_new: ");
	s=MIN(long,n,start);
	e=MAX(long,n,end);
	for(i=s; i<=e; i++){
	   if(i==lemon) fprintf(stderr,".");
	   else if(i==n) fprintf(stderr,"*");
	   else if(i >= s && i <= e){
		if(i >= start && i <= end && !null[i]) fprintf(stderr,"*");
		else fprintf(stderr,".");
	   }else fprintf(stderr," ");
	}
	s=MIN(long,n,start);
	e=MAX(long,n,end);
	fprintf(stderr," (%d-%d)",s,e);
	fprintf(stderr,"\nw_old = %d w_new = %d; c= %d\n",
			w_old,w_new,c);
	fprintf(stderr,"P = %g/%g = %g; n=%d; lemon = %d\n",b0,b2,P,n,lemon);
}
/****/
	return (bico(w_old-2,c-2)/bico(w_new-2,c-2));
}


Boolean	move_column_gibbs(gs_type G, long lemon, long t)
/*******************************************************************
  Sample a column to remove and then sample a column to replace it.
  if lemon = new site then don't bother to move ??? 
********************************************************************/
{
	ss_type	P=G->data;
	st_type S=G->sites;
	long	**site_freq;
	long	start,fend,end,i,j,o,d,n,k,flank,leng,newstart,s;
	double	*ratio,total,rand_no;
	fm_type M=G->model[t];
	char	c;
	a_type	A=SeqSetA(P);
	Boolean	left;

	double	weight;
	long	ncol;

	if(!G->test[0]){
	  ncol = nColsFModel(M);
	  NullSitesFModel(G->null, M);
	}
	flank = (G->maxlen[t] - LenFModel(M))/2;
	fend = end =  LenFModel(M) + flank;
	leng = LenFModel(M); site_freq = G->tmpfreq;
	ratio = G->tmpratio; 
	site_freq += flank; ratio += flank; site_freq[1] = NULL;
	for(start=1,total=1.0,i=0; i>= -flank; i--){
		site_freq[i] = GetSiteFreq(G,t,i-1);
		if(site_freq[i] == NULL) { start = i+1; break; }
		ratio[i] = RatioFModel(site_freq[i],lemon, M);
		if(!G->test[0]){
		  weight = mv_col_weight(ncol, leng, i, lemon, G->null);
		  ratio[i] *= weight;
		}
		total += ratio[i];
	}
	site_freq[1] = NULL;
	for(i=2; i<= end; i++){
	   if(NullSiteFModel(i,M)){
		site_freq[i] = GetSiteFreq(G,t,i-1);
		if(site_freq[i] == NULL) { end = i-1; break; }
		ratio[i] = RatioFModel(site_freq[i],lemon, M);
		if(!G->test[0]){
		  weight = mv_col_weight(ncol, leng, i, lemon, G->null);
		  ratio[i] *= weight;
		}
		total += ratio[i];
	   } else { site_freq[i] = NULL; }
	}
	rand_no = ((double)Random()/(double)LONG_MAX)*total;
	for(total=0.0, i=start; i <= end; i++){
	   if(NullSiteFModel(i,M)){
	      if((rand_no -= ratio[i]) <= 0.0) {
/****
  if(ratio[i] > 1e+100){
		NullSitesFModel(G->null, M);
		PutSites(stderr,t,S,NULL,G->null);
  }
/****/
		d = MvColumnFModel(site_freq[i],lemon,i,M);
		site_freq[i]= NULL;
		fprintf(stderr," (new=%d; old=%d; w=%d) %g",
			i,lemon,LenFModel(M), ratio[i]);
		if(LenFModel(M) < leng){	/* model shrunk */
			k = leng - LenFModel(M);
			for(j = 1; j <= k; j++) { ShrinkSites(t,S); }
			if(d != 0) ShiftSitesM(S, t, d);
		}else if(LenFModel(M) > leng){	/* model has grown */
			k = LenFModel(M) - leng;
			if(d != 0) ShiftSitesM(S, t, d);
			for(j = 1; j <= k; j++) { GrowSites(t,S); }
		} else if(d!=0) ShiftSitesM(S, t, d);
		for(j=-flank; j<=fend; j++){
			if(site_freq[j]!=NULL) {
				free(site_freq[j]); site_freq[j] = NULL;
			}
		}
		return TRUE; 
	      }
	   }
	}
	for(j=-flank; j<=fend; j++){
		if(site_freq[j]!=NULL) { 
			free(site_freq[j]); site_freq[j] = NULL;
		}
	}
	return FALSE;
}

Boolean	TransferColumn(gs_type G, long ntyp, fm_type *model)
{
	ss_type	P=G->data;
	st_type S=G->sites;
	long	**site_freq;
	double	*ratio,total,rand_no;
	long	t1,t2,start,leng1,leng2;
	long	fend,end,i,j,o,d,n,k,flank,lemon,newstart,s;
	char	c;
	fm_type	M1,M2;
	a_type	A=SeqSetA(P);
	bp_type	B1,B2;
	Boolean	left;

	
	do{	/* select a model for lemon */
		rand_no = (double)Random()/(double)LONG_MAX;
		rand_no *= (double) ntyp + 0.5;
		t1 = 1 + (long) rand_no;
	} while(t1 < 1 || t1 > ntyp);
	M1 = model[t1];
	if(G->prior != NULL) B1 = G->prior[t1];
	else B1 = NULL;

	do{	/* select a model for new column */
		rand_no = (double)Random()/(double)LONG_MAX;
		rand_no *= (double) ntyp + 0.5;
		t2 = 1 + (long) rand_no;
	} while(t2 < 1 || t2 > ntyp || t2 == t1);
	M2 = model[t2];
	if(G->prior != NULL) B2 = G->prior[t2];
	else B2 = NULL;

	leng1 = LenFModel(M1);
	lemon = LemonFModel(M1);
	if(lemon == 1 || lemon == leng1) return FALSE;
	/*** flank = (G->maxlen[t2] - LenFModel(M2))/2; /*** OLD ***/
	flank = 0; /*** need to weight columns if changes model size ***/
	fend = end =  LenFModel(M2) + flank;
	leng2 = LenFModel(M2);
	site_freq = G->tmpfreq; ratio = G->tmpratio;
	site_freq += flank; ratio += flank; site_freq[1] = NULL;
	for(start=1,total=1.0,i=0; i>= -flank; i--){
		site_freq[i] = GetSiteFreq(G,t2,i-1);
		if(site_freq[i] == NULL) { start = i+1; break; }
		ratio[i] = RatioFModel(site_freq[i],lemon, M1);
		if(G->prior != NULL) ratio[i] *= RatioBPrior(B2, B1);
		total += ratio[i];
	}
	site_freq[1] = NULL;
	for(i=2; i<= end; i++){
	   if(NullSiteFModel(i,M2)){
		site_freq[i] = GetSiteFreq(G,t2,i-1);
		if(site_freq[i] == NULL) { end = i-1; break; }
		ratio[i] = RatioFModel(site_freq[i],lemon, M1);
		if(G->prior != NULL) ratio[i] *= RatioBPrior(B2, B1);
		total += ratio[i];
	   } else { site_freq[i] = NULL; }
	}
	rand_no = (double)Random()/(double)LONG_MAX;
	rand_no *= total;
	for(total=0.0, i=start; i <= end; i++){
	   if(NullSiteFModel(i,M2)){
	      if((rand_no -= ratio[i]) <= 0.0) {
	/*****
		G->loglike[t] += log(RatioFModel(site_freq[i],lemon, M2));
	/*****/
	if(t1!=t2) {  /* t1 can't be == t2  - see above. */
	   /*** FIRST REMOVE OLD COLUMN ***/
		d = RmColumnFModel2(lemon, M1);
		G->ncol[t1]--;

		fprintf(stderr," (old=%d; length=%d) %g model %c\n",
			lemon,LenFModel(M1), ratio[i], t1+'A'-1);

		if(LenFModel(M1) < leng1){	/* model shrunk */
			k = leng1 - LenFModel(M1);
			for(j = 1; j <= k; j++) { ShrinkSites(t1,S); }
			if(d != 0) ShiftSitesM(S, t1, d);
		}else if(LenFModel(M1) > leng1){	/* model has grown */
			print_error("this should not happen");
		} else if(d!=0) ShiftSitesM(S, t1, d);

	  /*** THEN ADD NEW COLUMN ***/
		d = AddColumnFModel2(site_freq[i], i, M2);
		G->ncol[t2]++;
		site_freq[i]= NULL;

		fprintf(stderr," (new=%d; length=%d) %g model %c\n",
			i,LenFModel(M2), ratio[i],'A' + t2 - 1);

		if(LenFModel(M2) < leng2){	/* model shrunk */
			print_error("this should not happen");
		}else if(LenFModel(M2) > leng2){	/* model has grown */
			k = LenFModel(M2) - leng2;
			if(d != 0) ShiftSitesM(S, t2, d);
			for(j = 1; j <= k; j++) { GrowSites(t2,S); }
		} else if(d!=0) ShiftSitesM(S, t2, d);

	/*** FREE POTENTIAL COLUMNS ***/
		for(j=-flank; j<=fend; j++){
			if(site_freq[j]!=NULL) {
				free(site_freq[j]);
				site_freq[j] = NULL;
			}
		}
		return TRUE; 
	} else print_error("this should not happen!?");
	      }
	   }
	}
	for(j=-flank; j<=fend; j++){
		if(site_freq[j]!=NULL) { 
			free(site_freq[j]); site_freq[j] = NULL;
		}
	}
	return FALSE;
}


