#include "fmodel.h"

fm_type	MkFModel(Boolean *null, long length, double npseudo, long *counts, 
	a_type A)
/* create and return a fmodel with totsites = 0 and unfragmented columns. */
{
	fm_type	M;
	long	i,j,r;
	
	NEW(M,1,fmodel_type);
	M->A = A; M->totsites = 0; M->length = length; 
	if(null != NULL){
	   for(M->ncols=length,j=1; j<=length; j++) if(null[j]) M->ncols--;
	   if(M->ncols < 4) fmodel_error("model must have at least 3 columns");
	   M->maxlen = M->ncols * 5;
	   if(M->maxlen <= length) fmodel_error("model too sparce");
	} else { M->maxlen = length * 10; M->ncols = length; }
	M->start = MAX(long,1,(M->maxlen - length)/2);
	M->end = M->start + length - 1;
	M->npseudo = MAX(double,npseudo,0.0001);
	NEW(M->counts,nAlpha(A)+2,long);
	for(r=0; r<= nAlpha(A); r++) M->counts[r]= counts[r];
	NEW(M->freq,nAlpha(A)+2,double);
	NEW(M->observedS,nAlpha(M->A) +2, long);
	NEWP(M->observed, M->maxlen+1, long);
	NEWP(M->likelihood, M->maxlen+1, double);
	NEWP(M->target, M->maxlen+1, double);
	NEW(M->mark, M->maxlen+1, Boolean);
	NEW(M->sumlgamma, M->maxlen+1, double);
	NEW(M->tmp_val, M->maxlen+1, double);
	for(j = 0; j <= M->maxlen; j++){ M->observed[j] = NULL; }
	for(M->ncols = 0, i=1,j = M->start; j <= M->end; j++,i++){	
	  if(null == NULL || !null[i]){
	    NEW(M->observed[j],nAlpha(M->A) +2, long);
	    NEW(M->likelihood[j], nAlpha(M->A)+2, double);
	    NEW(M->target[j], nAlpha(M->A)+2, double);
	    M->ncols++;
	  }
	}
	NEW(M->Ps,nAlpha(M->A)+2,double);
	NEW(M->target[0],nAlpha(M->A)+2,double);
	NEW(M->observed[0],nAlpha(M->A) +2, long);
	InitFModel(M); 
	return M;
}

fm_type	CopyFModel(fm_type M)
/* return a fmodel with totsites = 0 and unfragmented columns. */
{
	fm_type	M2;
	Boolean	*null;
	long	i,j;
	
	NEW(null,M->length+2,Boolean);
	for(i=1,j = M->start; j <= M->end; j++,i++){	
		if(M->observed[j] == NULL) null[i] = TRUE;
		else null[i] = FALSE;
	}
	M2 = MkFModel(null, M->length, M->npseudo, M->counts, M->A);
	free(null);
	return M2;
}

void	InitFModel(fm_type M)
/* Initialize model M to no sites */
{
	long	j,b;
	double	d;

	M->totsites = 0; M->tot_cnts=0.0;
	for(b=1; b<= nAlpha(M->A); b++) M->tot_cnts += M->counts[b];
	for(b=1; b<= nAlpha(M->A); b++) {
		M->freq[b] = (double) M->counts[b]/(double) M->tot_cnts;
		M->Ps[b]= M->npseudo * M->freq[b]; 
	}
	for( j = M->start; j <= M->end; j++){	
	  if(M->observed[j]!=NULL){
	   for(b=0; b <= nAlpha(M->A); b++) M->observed[j][b] = 0;
	  }
	}
	ClearMarksFModel(M);
	M->recalc = M->update = TRUE;
}

fm_type NilFModel(fm_type M)
/* Destroy model M */
{
   long j;
   if(M!=NULL){ 
	free(M->freq); free(M->counts); free(M->Ps); 
	free(M->target[0]); free(M->observed[0]); free(M->mark);
	free(M->sumlgamma);
	for(j=M->start; j<=M->end; j++){
	  if(M->observed[j]!=NULL){
	    free(M->observed[j]); free(M->likelihood[j]);
	    free(M->target[j]);
	  }
	}
	free(M->observed); 
	free(M->observedS);
	free(M->target); free(M->likelihood); free(M->tmp_val);
	free(M);
   }
   return (fm_type) NULL;
}

double	LikelihoodFModel(register char *seq, register long pos, 
	register fm_type M)
/* Return the likelihood of site at pos in seq being in model M */
{
        register long j;
        register double p;
	if(M->update) update_fmodel(M);
        for(p=1.0,j=M->start; j<=M->end; j++,pos++){
                if(M->likelihood[j]!=NULL)p*=M->likelihood[j][(seq[pos])];
        }
        return p;
}

double	ProbFModel(register char *seq, register long pos, register double p,
	register fm_type M)
/* Return the relative probability of site at pos in seq being in model M */
{
        register long j;
        register double L;

	if(M->update) update_fmodel(M);
        for(L=1.0,j=M->start; j<=M->end; j++,pos++){
                if(M->likelihood[j]!=NULL) L*=M->likelihood[j][(seq[pos])];
        }
	return (L*p/((1.0 - p) + L*p));
}

void    RmFModel(char *seq, long site, fm_type M)
/* remove the segment at site in seq from model */
{
        long j;
        for(j=M->start; j<=M->end; j++,site++){
	    if(M->observed[j] != NULL) M->observed[j][seq[site]]--;
        }
	M->totsites--; M->recalc = M->update = TRUE;
}

void	Add2FModel(char *seq, long site, fm_type M)
/* Add the segment at site in seq to model */
{
	long j;
	for(j=M->start; j<=M->end; j++,site++){
	    if(M->observed[j] != NULL)
	    	M->observed[j][seq[site]]++;
	}
	M->totsites++; M->recalc = M->update = TRUE;
}

Boolean	NullSiteFModel(long s,fm_type M)
/* if site s == null in model M return TRUE; else return FALSE. */
{
	long	site;
	if(s < 1 || s > M->length) return TRUE;
	site = s + M->start - 1;
	if(M->observed[site] == NULL) return TRUE;
	else return FALSE;
}

long	NullSitesFModel(Boolean *null, fm_type M)
/* modifies null array so that null sites in model are TRUE and
   nonnull sites are FALSE; returns the length of model */
{
	long	i,j;

	for(i=1,j=M->start; j <= M->end; j++,i++){
		if(M->observed[j]==NULL) null[i] = TRUE;
		else null[i] = FALSE;
	}
	return M->length;
}

/************************ Column Move Routines ****************************/

void	ClearMarksFModel(fm_type M)
/** set all marks to FALSE **/
{
	long	j;

	for(j=0; j <= M->maxlen; j++) {
		M->mark[j] = FALSE;
	}
}

long	ChoiceLemonFModel(fm_type M)
/** Sample a column in model proportional to how BAD it is. */
{
	long	i,k;
	double	r,rand,total;
	long	*observed;

	for(total=0.0,k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL && !M->mark[i]){
		r = RatioFModel(M->observed[M->start],k,M);
		M->tmp_val[k] = r; total += r;
	     }
	}
	rand = ((double) Random()/(double) LONG_MAX)*total;
	for(k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL && !M->mark[i]){
		rand -= M->tmp_val[k];
		if(rand <= 0.0) { M->mark[i] = TRUE; return k; }
	     }
	}
	fmodel_error(" LemonFModel( )... this should not happen.");
}

long	ChoiceOrangeFModel(fm_type M)
/** Sample a column in model at random. */
{
	long	i,k;
	double	r,rand,total;
	long	*observed;

	for(total=0.0,k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL && !M->mark[i]) total += 1.0;
	}
	rand = ((double) Random()/(double) LONG_MAX)*total;
	for(k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL && !M->mark[i]){
		rand -= 1.0;
		if(rand <= 0.0) { M->mark[i] = TRUE; return k; }
	     }
	}
	fmodel_error(" LemonFModel( )... this should not happen.");
}

long	LemonFModel(fm_type M)
/** Sample a column in model proportional to how BAD it is. */
{
	double	r,rand,total;
	long	i,k,*observed;

	for(total=0.0,k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL){
		r = RatioFModel(M->observed[M->start],k,M);
		M->tmp_val[k] = r; total += r;
	     }
	}
	rand = ((double) Random()/(double) LONG_MAX)*total;
	for(k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL){
		rand -= M->tmp_val[k];
		if(rand <= 0.0) return k;
	     }
	}
	fmodel_error(" LemonFModel( )... this should not happen.");
}

long	OrangeFModel(fm_type M)
/** Sample a random column in model. */
{
	double	r,rand,total;
	long	i,k,*observed;

	for(total=0.0,k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL) total += 1.0;
	}
	/** total = (double) M->ncols; /***/
	rand = ((double) Random()/(double) LONG_MAX)*total;
	for(k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL) if((rand-= 1.0) <= 0.0) return k;
	}
	fmodel_error(" RandColFModel( )... this should not happen.");
}

double	RatioFModel(long *observed, long d, fm_type M)
/* Return the ratio of new to old observed.  */
{
	double	sum;
	long	b,j,*mobserved;

	if(observed == NULL || d < 1 || d > M->length) return 0.0;
	if(M->recalc){
	    for(j=0; j<=M->maxlen; j++) M->sumlgamma[j] = FMODEL_UNDEF;
	    M->recalc = FALSE;
	}
	j = d+M->start-1;
	if(M->sumlgamma[j] == FMODEL_UNDEF){
	  if((mobserved = M->observed[j])==NULL) return 0.0;
	  for(M->sumlgamma[j]=0.0,b=1;b<= nAlpha(M->A); b++) {
	     if(M->Ps[b] != 0.0){
		M->sumlgamma[j] += lgamma(((double)mobserved[b]+M->Ps[b]));
	     }
	  }
	} 
	for(sum=0.0,b=1;b<= nAlpha(M->A); b++) {
	   if(M->Ps[b]!=0.0) sum+=lgamma(((double)observed[b]+M->Ps[b]));
	}
	return exp(sum-M->sumlgamma[j]);
}

long	MvColumnFModel(long *observed, long lemon, long pos, fm_type M)
/*  Remove column at position lemon and add observed to pos in M */
{
	long	oms,ms,oldstart;
	
	oldstart = M->start;
	if(!RmColumnFModel(lemon, M)) fmodel_error("can't remove column!?");
	if(oldstart != M->start) pos += oldstart - M->start;
	oms = M->start; 	/* center_model may alter oms */
	ms = AddColumnFModel(observed, pos, M);
	if(ms!=oms) oldstart += ms - oms; /* re-adjust oldstart */
	return  M->start - oldstart;
}

long	RmColumnFModel2(long lemon, fm_type M)
/** remove column and return change in start site **/
{
	long	oldstart;
	oldstart = M->start;
	if(!RmColumnFModel(lemon, M)) fmodel_error("can't remove column!?");
	M->ncols--;
	return  M->start - oldstart;
}

Boolean RmColumnFModel(long pos, fm_type M)
/*  if possible position "pos" is removed from model and TRUE is returned;
    if position "pos" cannot be removed FALSE is returned */
{
	long	site;

	if(M->length == 1) fmodel_error("zero length model not allowed.");
	site = pos + M->start - 1;
	if(site < M->start || site > M->end)  return FALSE;
	else if(M->observed[site]!=NULL){
		free(M->observed[site]); free(M->likelihood[site]); 
		free(M->target[site]); 
		M->observed[site] = NULL; M->likelihood[site] = NULL; 
		M->target[site] = NULL; 
		if(site == M->start){		/* shrink from left */
			while(M->observed[site]==NULL) site++;
			M->start=site; 
			M->length = M->end - M->start + 1;
		}else if(site == M->end){	/* shrink from right */
			while(M->observed[site]==NULL) site--;
			M->end=site; 
			M->length = M->end - M->start + 1;
		}
		M->ncols--;
		M->update = TRUE;
		return TRUE;
	} else return FALSE;
}

long	AddColumnFModel2(long *observed, long pos, fm_type M)
/** remove column and return change in start site **/
{
	long	d,oldstart;
	oldstart=AddColumnFModel(observed, pos, M);
	return  M->start - oldstart;
}

long	AddColumnFModel(long *observed, long pos, fm_type M)
/********************************************************************
 Add the column observed at pos in model M.  The center_model( ) 
 function may change the value of M->start therefore AddColumnFModel( )
 returns the value of the (possibly changed) old M->start to the 
 calling environment which may need this value.
*********************************************************************/
{
	long	site,b,oldstart=M->start;

	site = pos + M->start - 1;
	if(site <= 0 || site >= M->maxlen){
		fprintf(stderr,"\ncentering model\n");	/*****/
		center_model(site, M);
		site = pos + M->start - 1;
		oldstart=M->start;
	}
	if(site < M->start){ 
		M->start = site; 
		M->length = M->end - M->start + 1;
	} else if(site > M->end){	
		M->end = site; 
		M->length = M->end - M->start + 1;
	} else if(M->observed[site]!=NULL) 
		fmodel_error("attempt to add column where one exists!?");
	NEW(M->likelihood[site],nAlpha(M->A) +2, double);
	NEW(M->target[site],nAlpha(M->A) +2, double);
	M->observed[site] = observed; M->ncols++;
	M->update = TRUE;
	return oldstart;
}

void	ShiftFModel(long *observed, Boolean left, fm_type M)
{
	if(left){
		RmColumnFModel(1, M);
		AddColumnFModel(observed, M->length+1, M);
	} else {
		RmColumnFModel(M->length, M);
		AddColumnFModel(observed, 0, M);
	}
	M->update = TRUE;
}

/****************************** Output ********************************/
void	PutFModel(FILE *fptr, fm_type M)
/* Report the current frequency model. */
{
	long	j,b,r,pos,v;
	double	i,total,info,p,q,maxinfo,ave;
	long	*observed;

	if(M->update) update_fmodel(M);
/*** freq model ****/
	fprintf(fptr,"Motif model (residue frequency x 100):\n");
	fprintf(fptr,"POS  ");
	for(pos=1,j=M->start; j<= M->end; j++,pos++){
	  if(M->observed[j]!=NULL){
	    for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
	      total += (double) M->observed[j][b] + M->Ps[b];
	      if(j==M->start) fprintf(fptr,"%3c", AlphaChar(b, M->A));
	    }
	    if(j==M->start) fprintf(fptr,"  Info\n");
	    fprintf(fptr,"%4d ",pos);
	    for(info=0.0,b = 1; b <= nAlpha(M->A); b++){
	    	p = ((double) M->observed[j][b]+M->Ps[b])/total;
		if(p > 0.0){
			q = M->freq[b];
			i = p*log(p/q)/log(2.0);
			info += i;
		}
		v = (long)floor(10*i+0.5);
		if(v<=0) fprintf(fptr,"  .");
	    	else fprintf(fptr,"%3d", (long)(100*p+0.5));
	    }
	    fprintf(fptr,"   %1.1f\n", info);
	  }
	}
	fprintf(fptr,"non-\nsite ");	
	for(total=0.0,b = 1; b <= nAlpha(M->A); b++){
		total += (double) M->observedS[b];
	}
	for(b = 1; b <= nAlpha(M->A); b++){
	    	p = (double) M->observedS[b]/total;
		if(p > 0.0){ q = M->freq[b]; i = p*log(p/q)/log(2.0); }
		else i=0.0;
		v = (long)floor(10*i+0.5);
		if(v<=0) fprintf(fptr,"  .");
		else {
			r = (long) (100.0 * M->freq[b]);
	    		fprintf(fptr,"%3d", r);
		}
	}
	fprintf(fptr,"\nsite ");	
	for(b = 1; b <= nAlpha(M->A); b++){
	    	p = (double) M->observedS[b]/total;
		if(p > 0.0){ q = M->freq[b]; i = p*log(p/q)/log(2.0); }
		else i=0.0;
		v = (long)floor(10*i+0.5);
		if(v<=0) fprintf(fptr,"  .");
		else {
		   r = (long) (100.0 * (double) M->observedS[b]/total);
	    	   fprintf(fptr,"%3d", r);
		}
	}
	fprintf(fptr,"\n%d columns\n\n",M->ncols);
/*** Information contributions ****/
	fprintf(fptr,
	  "Information (relative entropy) contribution in tenth bits:\n");
	fprintf(fptr,"POS  ");
	for(pos=1,j=M->start; j<= M->end; j++,pos++){
	  if(M->observed[j]!=NULL){
	    for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
	      total += (double) M->observed[j][b] + M->Ps[b];
	      if(j==M->start) fprintf(fptr,"%3c", AlphaChar(b, M->A));
	    }
	    if(j==M->start) fprintf(fptr,"  Info\n");
	    fprintf(fptr,"%4d ",pos);
	    for(info=0.0,b = 1; b <= nAlpha(M->A); b++){
	    	p = ((double) M->observed[j][b]+M->Ps[b])/total;
		if(p > 0.0){
			q = M->freq[b];
			i = p*log(p/q)/log(2.0);
			info += i;
		} else i = 0.0;
		v = (long)floor(10*i+0.5);
		if(v<=0) fprintf(fptr,"  .");
		else fprintf(fptr,"%3d", v);
	    }
	    v = (long)floor(10*info+0.5);
	    fprintf(fptr,"  %3d\n", v);
	  }
	}
	fprintf(fptr,"\nsite ");	
	for(total=0.0,b = 1; b <= nAlpha(M->A); b++){
		total += (double) M->observedS[b];
	}
	for(info=0.0,b = 1; b <= nAlpha(M->A); b++){
	    	p = (double) M->observedS[b]/total;
		if(p > 0.0){
			q = M->freq[b];
			i = p*log(p/q)/log(2.0);
			info += i;
		} else i=0.0;
		v = (long)floor(10*i+0.5);
		if(v<=0) fprintf(fptr,"  .");
		else fprintf(fptr,"%3d", v);
	}
	v = (long)floor(10*info+0.5);
	fprintf(fptr,"  %3d\n\n", v);
}

#include "blosum62.h"
void	PutLogFModel(FILE *fptr, fm_type M)
/* Report the current frequency model. */
{
	long	i,j,b,r,s,pos,max_r;
	long	info,end;
	double	total,p,q,ave,tot_info,val;
	long	*observed;

	char	**res;
	dh_type	H;
	keytyp  key;
	
	NEWP(res,nAlpha(M->A)+2,char);
	H = dheap(nAlpha(M->A),3);
	for(r=1; r<=nAlpha(M->A); r++){
	   NEW(res[r],nAlpha(M->A)+2,char);
	   for(s=1; s<=nAlpha(M->A); s++){
	     if(s != r){
		key = (keytyp) -blosum62[r][s];
		insrtHeap(s,key,H);
	     }
	   }
	   for(i = 0; TRUE; i++){ 
		if((s=delminHeap(H)) != NULL) res[r][i] = s;
		else break;
	   }
	}

	if(M->update) update_fmodel(M);
	fprintf(fptr,"POS  ");
	for(pos=1,j=M->start; j<= M->end; j++,pos++){
	  if(M->observed[j]!=NULL){
	    for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
	      /** total += (double) M->observed[j][b] + M->Ps[b]; /***/
	      total += (double) M->observed[j][b];
	      if(j==M->start) fprintf(fptr,"%3c", AlphaChar(b, M->A));
	    }
	    if(j==M->start) fprintf(fptr,"  Info\n");
	    fprintf(fptr,"%4d ",pos);
	    for(tot_info=0.0,b = 1; b <= nAlpha(M->A); b++){
	    	/** p = ((double) M->observed[j][b]+M->Ps[b])/total; /***/
		if(M->observed[j][b] > 0){
	    		p = ((double) M->observed[j][b])/total;
			q = M->freq[b];
			val = log(p/q)/log(2.0);
			if(val > 0.0){
				key = (keytyp) -val;
				insrtHeap(b,key,H);
			}
			info = (long) floor((2*val) + 0.5);
			tot_info += p*log(p/q)/log(2.0);
	    		fprintf(fptr,"%3d", info);
		} else fprintf(fptr,"   ");
	    }
	    fprintf(fptr,"   %1.1f\n", tot_info);
	    fprintf(fptr,"Column: ");
	    for(i = 0; TRUE; i++){ 
		if((r=delminHeap(H)) != NULL){
		  if(i!=0) fprintf(fptr," %c", AlphaChar(r,M->A));
		   else max_r = r;
		} else break;
	    }
	    end = MIN(long,i,nAlpha(M->A)-1);
	    fprintf(fptr," (major = %c)\nBlosum: ",AlphaChar(max_r,M->A));
	    for(i = 0; i < nAlpha(M->A)-1; i++){
		fprintf(fptr," %c",AlphaChar(res[max_r][i],M->A)); 
	    }
	    fprintf(fptr,"\n");
	  }
	}
	fprintf(fptr,"\n\n%d columns\n",M->ncols);
	for(r=1; r<=nAlpha(M->A); r++) free(res[r]);
	free(res);
	Nildheap(H);
}

/****************************** Statistics ********************************/
double	LnMapFModel(fm_type M)
/* return Log of Maximum a posteriori probability (MAP) for model */
{
	long	b,i,k;
	double	value,v,total;
	
	if(M->update) update_fmodel(M);
	k = M->length;
        /******************* column width weight ******************/
        value = -lnbico(M->length - 2, M->ncols - 2);
        /***********************************************************/
	for(i=M->start; i <= M->end; i++) {
	  if(M->observed[i] != NULL){
	    for(b=1;b<= nAlpha(M->A); b++) {
	      if(M->Ps[b] != 0.0){
		value += lgamma((double)M->observed[i][b]+M->Ps[b]);
	      }
	    }
	    value -= lgamma((double)(M->totsites + M->npseudo));
	  } else k--;
	}
	v = lgamma((double)M->npseudo);
	for(total = 0.0, b=1;b<= nAlpha(M->A); b++) {
	   if(M->Ps[b] != 0.0){
		value += lgamma(((double) M->observed[0][b]+M->Ps[b]));
		total += (double) M->observed[0][b];
		v -= lgamma((double)M->Ps[b]);
	   }
	}
	value -= lgamma((double)(total + M->npseudo));
	v *= ((double) k + 1.0);
	return (value+v);
}

double	LnMapFModelNull(fm_type M)
{
	long	b,i,k;
	double	value,v,total;
	
	if(M->update) update_fmodel(M);
	k = M->length;
	for(value=0.0,i=M->start; i <= M->end; i++) {
	  if(M->observed[i] != NULL){
	    for(b=1;b<= nAlpha(M->A); b++) {
	      if(M->Ps[b] != 0.0){
		value += lgamma(M->Ps[b]);
	      }
	    }
	    value -= lgamma((double)(M->npseudo));
	  } else k--;
	}
	v = lgamma((double)M->npseudo);
	for(total = 0.0, b=1;b<= nAlpha(M->A); b++) {
	   if(M->Ps[b] != 0.0){
		value += lgamma((M->counts[b]+M->Ps[b]));
		total += M->counts[b];
		v -= lgamma((double)M->Ps[b]);
	   }
	}
	value -= lgamma((double)(total + M->npseudo));
	v *= ((double) k + 1.0);
	return (value+v);
}

double	LogLikeFModel(fm_type  M)
/* log likelihood for non-converged (NC) model */
{
	double	term1,term1b,n;
	long	j,b;

	if(M->update) update_fmodel(M);
	for(term1=0.0, j=M->start; j<= M->end; j++){
	   if(M->observed[j]!= NULL){
		for(b=1; b<= nAlpha(M->A); b++){
		    if((n=(double)M->observed[j][b]) > 0.0){
			term1 += n * log(n/(M->totalS - M->npseudo));
		    }
		}
	   }
	} /*** null model ****/
        for(term1b=0.0, b=1; b<= nAlpha(M->A); b++){
            if((n=M->observed[0][b]) > 0.0){
                term1b += n * log(n/(M->total0 - M->npseudo));
	    }
	}
	return (1.4427*(term1 + term1b));
}

double	InfoColFModel(long *observed, fm_type M)
/* Report the information content of column "observed" relative to model M. */
{
	long	b;
	double	total,info,p,q;

	if(observed != NULL){
	    for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
		total += (double)observed[b] + M->Ps[b];
	    }
	    for(info=0.0,b = 1; b <= nAlpha(M->A); b++){
	    	p = ((double)observed[b]+M->Ps[b])/total;
		if(p > 0.0){
			q = M->freq[b];
			info += p*log(p/q); /* info in nats */
		}
	    }
	    return 1.442695*info;	/* info in bits */
	}
	return 0.0;
}

double	InfoFModel(fm_type M)
/* return the relative entropy of model M in bits. */
{
	long	j;
	char	b;
	double	E;
	
	if(M->update) update_fmodel(M); 
        for(E=0.0,j=M->start; j<=M->end; j++){
	    if(M->likelihood[j]!= NULL){
		for(b = 1; b <= nAlpha(M->A); b++){
		  if(M->Ps[b] > 0.0){
              	    E += M->target[j][b]*log(M->likelihood[j][b]);
		  } 
		}
	    }
	}
	return E*1.442695;
}

/**************************** Simulations *****************************/
Boolean	NonSiteSegFModel(char *seg, fm_type M)
/* Get a random segment drawn from the null model.  */
{
	long     i,j,c;
	double  r;

	if(M->update) update_fmodel(M);
        for(i=1,j=M->start; j<=M->end; j++,i++){
		r = (double) Random()/(double) LONG_MAX; /* 0 <= r <= 1 */
		for(c=1; c <= nAlpha(M->A); c++){
			if((r-=M->freq[c]) <= 0.0){ seg[i]=c; break; }
	    	}
	}
}

Boolean	GetSegFModel(char *seg, fm_type M)
/* Get a random segment drawn from the model M. Note: positions 
   corresponding to null columns are drawn from the non-site model */
{
	long     i,j,c;
	double  r;

	if(M->update) update_fmodel(M);
        for(i=1,j=M->start; j<=M->end; j++,i++){
	    r = (double) Random()/(double) LONG_MAX; /* 0 <= r <= 1 */
            if(M->observed[j]!=NULL){
		for(c=1; c <= nAlpha(M->A); c++){
		    if((r-=M->target[j][c]) <= 0.0){ seg[i]=c; break; }
	    	}
	    } else {
		for(c=1; c <= nAlpha(M->A); c++){
		    if((r-=M->freq[c]) <= 0.0){ seg[i]=c; break; }
	    	}
	    }
	}
}

void	SetPseudoFModel(double npseudo, fm_type M)
{ 	char b;
	double d;
	M->npseudo = npseudo;
	d = 1.0/(double)nAlpha(M->A);
	for(b=1; b<= nAlpha(M->A); b++) {
		/*** M->Ps[b]= M->npseudo * d;/****/
		M->Ps[b]= M->npseudo * M->freq[b]; /****/
	}
	M->update=TRUE;
}

long	*RandColFModel(fm_type M)
/* return a random column for M using nonsite frequencies */
{
	double	rand_no,sum,total;
	long	*observed;
	long	i,b,T;

	NEW(observed,nAlpha(M->A)+2,long);
	for(total=0.0,b=1; b<=nAlpha(M->A); b++){
		total += M->freq[b];
		observed[b] = 0;
	}
	for(T=i=0; i < M->totsites; i++){
	   rand_no  =  (double) Random()/(double) LONG_MAX;
	   rand_no *= total;
	   for(sum=0.0,b=1; b<=nAlpha(M->A); b++){
		if((sum+=M->freq[b]) >= rand_no){
			observed[b]++;
			T++; break;
		}
	   }
	   if(sum < rand_no) fmodel_error("this should not occur!?!\n");
	}
	/** fprintf(stderr,"total residues = %d; ",T); /****/
	return observed;
}

/********************************* PRIVATE ********************************/

void    fmodel_error(char *s) { fprintf(stderr,"fmodel: %s\n",s); exit(1); }

long	center_model(long site, fm_type M)
/************************************************************************
 Recenters the model so that don't run off end of array. site is location
 of site to be added that is beyond on of the ends. 
 CAUTION: This operation clears all marks and sumlgammas.
 ************************************************************************/
{
	long	length,flank,start,end,i,j,k;

	if(site <= 0){
		length = M->end - site +1;
		if(length >= M->maxlen-2) 
			fmodel_error("added site creates model > maxleng.");
	} else if(site >= M->maxlen) {
		length = site - M->start +1;
		if(length >= M->maxlen-2) 
			fmodel_error("added site creates model > maxleng.");
	} else return;
	ClearMarksFModel(M); M->recalc = TRUE;
	flank = M->maxlen - length;
	start = flank/2; end = start + M->length - 1;
	if(start <= 0) fmodel_error("centering error; start < 1.");
	if(start < M->start){		/*** :...new..old: ***/
	    for(i = start,j= M->start; i <= end; i++,j++){
		M->observed[i] = M->observed[j];
		M->likelihood[i] = M->likelihood[j];
		M->target[i] = M->target[j];
		M->observed[j] = NULL; M->likelihood[j] = NULL;
		M->target[j] = NULL;
	    }
	} else if(start > M->start){  /*** :old..new...: ***/
	    for(i = end,j= M->end; i >= start; i--,j--){
		M->observed[i] = M->observed[j];
		M->likelihood[i] = M->likelihood[j];
		M->target[i] = M->target[j];
		M->observed[j] = NULL; M->likelihood[j] = NULL;
		M->target[j] = NULL;
	    }
	}
	M->start = start; M->end = end;
}

void	update_Ps_fmodel(register fm_type M)
{
	register long	j,b;

	for(b=1; b <= nAlpha(M->A); b++) {
		M->Ps[b] = (double)M->totsites*M->npseudo*M->freq[b];
	}
}

void	update_fmodel(register fm_type M)
/* Normalize frequencies to avoid overflow */
{
	register long	j,b;
	register a_type	A=M->A;

	/*** DETERMINE REGULAR AND NORMALIZED NONSITE FREQUENCIES ***/
	for(b=1; b <= nAlpha(A); b++) M->observedS[b]=0;
	for(j=M->start; j<=M->end; j++){
	    if(M->observed[j] != NULL){
		for(b=1; b <= nAlpha(A); b++){
			M->observedS[b] += M->observed[j][b];
		}
	    }
	}
	for(M->total0=0.0,b=1; b <= nAlpha(A); b++) {
		M->observed[0][b] = M->counts[b] - M->observedS[b];
		M->total0 += (double) M->observed[0][b]+M->Ps[b];
	}
	for(b=1; b <= nAlpha(A); b++){	/** nonsite freqs **/
	   if(M->Ps[b] > 0.0){		/** if residue present **/
             M->target[0][b] = (M->observed[0][b]+M->Ps[b])/M->total0;
	   }
	}
	M->totalS = (double) (M->totsites + M->npseudo);
        for(j=M->start; j<=M->end; j++){
	    if(M->observed[j] != NULL){
               	M->likelihood[j][0] = 0.05; /*** 'X' residues ***/
        	for(b=1; b <= nAlpha(A); b++){
               	    M->target[j][b]=(M->observed[j][b]+M->Ps[b])/M->totalS;
               	    M->likelihood[j][b] = M->target[j][b]/M->target[0][b];
           	}
	    }
        }
	M->update = FALSE; 
}

