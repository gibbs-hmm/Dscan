/****
A Gibbs Sampler algorithm for finding multiple sites in multiple sequences 
****/
#include "gibbs.h"

Boolean	Metropolis(gs_type G, long t, fm_type M) 
{return ShiftGibbs(G,t,M);}

Boolean	ShiftGibbs(gs_type G, long t, fm_type M)
/* like Metropolis( ) only allows shifts into other elements
   of same type */
{
	ss_type	P=G->data;
	st_type S=G->sites;
	long	*site_freq[2][15],*stfreq,b,p,len;
	long	end[2],i,j,o,d,d2,n,k,shift,t2;
	double	ratio[2][15], total,rand_no;
	Boolean	left[2],*overlap;
	e_type	E;
	char	c;

	shift = MIN(long,LenFModel(M)/2,10);
	shift = MAX(long,shift,1);
	NEW(overlap,shift+3,Boolean);
	end[0] = shift; left[0]=FALSE;	/* note: right == 0 == FALSE */
	end[1] = shift; left[1]=TRUE;	/* note: left == 1 ==TRUE */
	len = SiteLen(t,S);
	for(total=1.0,i=0; i<2; i++){
	   ratio[i][0] = 1.0;
	   for(d=1; d<= end[i]; d++){
		/*** = GetSiteFreq ***/
		MEW(stfreq, nAlpha(SeqSetA(P))+2,long);
		for(b=0;b<= nAlpha(SeqSetA(P)); b++) stfreq[b] = 0;
		if(left[i]){	/*** [__site__] ===> ***/
		  d2 = d + len - 1;
		  for(n=1; n <=NSeqsSeqSet(P) && stfreq!= NULL; n++){
		    E = SeqSetE(n,P);
		    for(k=1; k<=nSites(t,n,S); k++){
			p = SitePos(t,n,k,S) + d2;
			t2 = TypeSite(n,p,S);
			overlap[d]=OpenPos(n,p,S);
			if(p<=(long)SqLenSeqSet(n,P) && 
			    (overlap[d] || t2 == t || t2==BlockedSite(t,S))){
			   b=XnuSeq(p,E); stfreq[b]++; 
			   /** b=ResSeq(p,E); stfreq[b]++;  /****/
			} else { free(stfreq); stfreq = NULL; break; }
		    }
		  }
		} else {	/*** <=== [__site__] ***/
		  d2 = -d;
		  for(n=1; n <=NSeqsSeqSet(P) && stfreq!= NULL; n++){
		    E = SeqSetE(n,P);
		    for(k=1; k<=nSites(t,n,S); k++){
			p = SitePos(t,n,k,S) + d2;
			t2 = TypeSite(n,p,S);
			overlap[d]=OpenPos(n,p,S);
			if(p>=1 && (overlap[d] || t2==BlockedSite(t,S))){
			   b=XnuSeq(p,E); stfreq[b]++; 
			   /** b=ResSeq(p,E); stfreq[b]++; /***/
			} else { free(stfreq); stfreq = NULL; break;}
		    }
		  }
		}
		/*** end GetSiteFreq ***/
		if(stfreq == NULL) { end[i] = d-1; break; }
		else {
		  site_freq[i][d] = stfreq;
		  if(left[i]) d2 = d;
		  else d2 = SiteLen(t,S) - d + 1;
		  ratio[i][d]=RatioFModel(site_freq[i][d],d2, M);
		  ratio[i][d] *= ratio[i][d-1];
		  total += ratio[i][d];
		}
	   }
	}
	rand_no = ((double)Random()/(double)LONG_MAX)*total;
	for(total=0.0, i=0; i < 2; i++){
	   if(left[i]) c= '+'; else c = '-';
	   for(d=1; d<= end[i]; d++){
	      if((total += ratio[i][d]) >= rand_no) {
		fprintf(stderr,"\r[");
		for(j=1; j<=d; j++){
		    fprintf(stderr,"%c", c);
		    ShiftFModel(site_freq[i][j], left[i], M); 
		    if(overlap[d]){
			ShrinkSites(t,S);
			ShiftSites(S,t,left[i]);;
			GrowSites(t, S);
		    } else { ShiftSites(S,t,left[i]);; }
		}
		fprintf(stderr,"] %g ",ratio[i][d]);
		for(j=d+1; j<=end[i]; j++){
			if(site_freq[i][j]!=NULL) free(site_freq[i][j]);
		}
		for(o=(i+1)%2,j=1; j<=end[o]; j++){
			if(site_freq[o][j]!=NULL) free(site_freq[o][j]);
		}
		free(overlap);
		return TRUE; 
	      }
	   }
	}
	for(i=0; i<2; i++){
	   for(d=1; d<= end[i]; d++){
		if(site_freq[i][d]!=NULL) free(site_freq[i][d]);
	   }
	}
	free(overlap);
	return FALSE;
}

