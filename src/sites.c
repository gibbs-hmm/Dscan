#include "sites.h"

st_type	MkSites(long ntyps, long *len_elem, ss_type data)
/*********************************************************************
 Create and return site structure containing no sites. 
							 nsites: A B C
  seq[1]    = .................................................* 0 0 0
  seq[2]    = ...........................................*	 0 0 0
     :
  seq[nseq] = ...............................................*   0 0 0
**********************************************************************/
{
    long		n,m,t,max;
    st_type	S;

    if(ntyps > MAX_NO_TYPESITES) sites_error("Too many types.");
    NEW(S,1,sites_type);
    S->data = data;
    S->ntyp = ntyps;
    S->nseq = NSeqsSeqSet(data);
    S->len_seq = LengthsSeqSet(data);
    NEWP(S->type, S->nseq+1, char);
    NEW(S->pos, S->nseq+1, ol_type);
    for(max=0,n = 1; n <= S->nseq; n++){
	max = MAX(long,max,S->len_seq[n]);
	NEW(S->type[n],(S->len_seq[n]+2),char);
	S->type[n][(S->len_seq[n]+1)] = ENDTYPESITE;
	S->pos[n] = Olist(S->len_seq[n]+1);
    }
    NEW(S->tmp,max+1,long);
    NEW(S->len_elem, ntyps+1, long);
    NEW(S->totsites, ntyps+1, long);
    NEWP(S->nsites,ntyps+1,long);
    NEWPP(S->pos_prob,ntyps+1,double);
    NEWPP(S->site_pos, ntyps+1,unsigned short);
    for(t = 1; t <= ntyps; t++) {
	S->totsites[t] = 0;
        S->len_elem[t] = len_elem[t];
        NEW(S->nsites[t], S->nseq+1,long);
        NEWP(S->pos_prob[t],S->nseq+1,double);
	NEWP(S->site_pos[t], NSeqsSeqSet(data)+1,unsigned short);
	for(n = 1; n <= NSeqsSeqSet(data); n++) {
		S->nsites[t][n] = 0;
        	NEW(S->pos_prob[t][n],(S->len_seq[n]+2),double);
                m = ((long) SqLenSeqSet(n,data)/(long)len_elem[t]) +1;
                NEW(S->site_pos[t][n],m+2,unsigned short);
	}
	if(t==1) S->maxinc = S->len_elem[t] - 1;
	else S->maxinc = MIN(long,S->maxinc,(S->len_elem[t]-1));
    }
    S->maxinc = MAX(long,1,S->maxinc);
    return S;
}

Boolean	ShuffleSites(st_type S)
/* shuffle sites retaining the same number of each type in 
   each sequences **/
{
	long	t,n,k,nsite,s;
	
    for(t = 1; t <= S->ntyp; t++){
	for(n = 1; n <= S->nseq; n++) {
	    nsite = nSites(t,n,S);
	    for(k= 1; k<= nsite; k++) S->tmp[k] = S->site_pos[t][n][k];
	    for(k= 1; k<= nsite; k++) {
		/** WARNING: VacateSite( ) modifies S->site_pos **/
		VacateSite(t, n, S->tmp[k], S);
	    }
	    for(k= 1; k<= nsite; k++) {
		if(AddRandomSite(t,n,100000,S)==0) return FALSE;
	    }
	}
    }
    return TRUE;
}

st_type	ReadRandomSites(long ntyps, long *len_elem, ss_type P)
/* Ask how many sites in each sequence, give option to all have 
   the same number & the length of the sites.  */
{
	st_type	S;
	char	c;
	long	t,n,m,i,nsite_all,s;

    S = MkSites(ntyps, len_elem, P);
    for(t = 1; t <= ntyps; t++){
	fprintf(stderr,"For type %ld sites...\n",t);
	GETCHAR("\tDo all seqs have the same # of sites (Y/N)?",&c);
	if(c == 'Y' || c == 'y'){
	   do {
	   	GETINT("\tHow many sites in each seq?", &nsite_all);
	   } while(nsite_all < 1);
	   fprintf(stderr,"%ld sites for each\n",nsite_all);
	}
	for(n = 1; n <= NSeqsSeqSet(P); n++) {
	   if(c == 'Y' || c == 'y'){
	        for(i= 1; i<= nsite_all; i++){
		   if(AddRandomSite(t,n,100000,S) == 0){
		      sites_error("too many sites; try a smaller number");
		   }
		}
	   } else {
	      do{
		fprintf(stderr, "\tsequence %ld:",n);
		GETINT(" how many sites?", &m);
	      } while(m < 1);
	      for(i= 1; i<= m; i++){
		if((s=AddRandomSite(t,n,100000,S)) == 0){
		    sites_error("too many sites; try a smaller number");
		}
	      }
	   }
	}
   }
   return S;
}

/****************************** Archive Sites **************************/
sti_typ	ArchiveSites(st_type S)
{
	sti_typ	X;
	long	i,t,n,s,m;

    NEW(X,1,sites_info_type);
    X->ntyp = S->ntyp;
    X->data = S->data;
    NEW(X->len_elem, X->ntyp+1, long);
    NEWP(X->nsites,X->ntyp+1,long);
    NEWPP(X->site_pos, X->ntyp+1,unsigned short);
    for(t = 1; t <= X->ntyp; t++) {
        X->len_elem[t] = S->len_elem[t];
        NEW(X->nsites[t], S->nseq+1,long);
	NEWP(X->site_pos[t], NSeqsSeqSet(X->data)+1,unsigned short);
	for(n = 1; n <= NSeqsSeqSet(X->data); n++) {
		X->nsites[t][n] = S->nsites[t][n];
                m = ((long) SqLenSeqSet(n,X->data)/(long)S->len_elem[t]) +1;
                NEW(X->site_pos[t][n],m+2,unsigned short);
		for(i = 1; i <= S->nsites[t][n]; i++) {
			X->site_pos[t][n][i] = S->site_pos[t][n][i];
		}
	}
    }
    return X;
}

Boolean	SameArchiveSites(sti_typ X1, sti_typ X2)
{
	long	i,t,n;

    if(X1->ntyp != X2->ntyp || X1->data != X2->data) return FALSE;
    for(t = 1; t <= X1->ntyp; t++) {
	if(X1->len_elem[t] != X2->len_elem[t]) return FALSE;
	for(n = 1; n <= NSeqsSeqSet(X1->data); n++) {
	   if(X1->nsites[t][n] != X2->nsites[t][n]) return FALSE;
	   for(i = 1; i <= X1->nsites[t][n]; i++) {
		if(X1->site_pos[t][n][i] != X2->site_pos[t][n][i]) return FALSE;
	   }
	}
    }
    return TRUE;
}

st_type	ExtractSites(sti_typ X)
{
	st_type	S;
	long	t,n,i;

    S = MkSites(X->ntyp, X->len_elem, X->data);
    for(t = 1; t <= X->ntyp; t++) {
	for(n = 1; n <= NSeqsSeqSet(X->data); n++) {
		for(i = 1; i <= X->nsites[t][n]; i++) {
			AddSite(t, n, X->site_pos[t][n][i], S);
		}
	}
    }
    return S;
}

void	NilArchiveSites(sti_typ X)
{
    long i,t,n;

    for(t = 1; t <= X->ntyp; t++) {
	for(n = 1; n <= NSeqsSeqSet(X->data); n++) free(X->site_pos[t][n]);
        free(X->nsites[t]); free(X->site_pos[t]);
    }
    free(X->len_elem); free(X->nsites); free(X->site_pos);
    free(X);
}

/**************************** end Archive Sites ************************/
void	InitSites(st_type S)
/* inititalize S to contain no sites; i.e., vacate sites. */
{
    long		n,t,i;

    for(n = 1; n <= S->nseq; n++) {
	ClearOlist(S->pos[n]);
	for(i = 0; i<= S->len_seq[n]; i++) S->type[n][i] = VACANT;
    }
    for(t = 1; t <= S->ntyp; t++) {
	for(n = 1; n <= S->nseq; n++) {
	   S->nsites[t][n] = 0;
	   S->site_pos[t][n][1] = 0;
	   S->pos_prob[t][n][0] = (double) S->len_seq[n];
	   for(i = 1; i<= S->len_seq[n]; i++) S->pos_prob[t][n][i] = 1.0;
	}
   }
}

st_type CopySites(st_type S)
/* Create a Null "copy" of S; i.e., S2 has the same number and lengths of 
   sequences and the same number and lengths of elements - 
   WARNING: NO SITE ARE ADDED. */
{ return MkSites(S->ntyp,S->len_elem, S->data); }

void	NilSites(st_type S)
/* destroy sites structure S. */
{
    long	t,n;

   free(S->len_elem); free(S->len_seq); free(S->tmp); free(S->totsites); 
   for(n = 1; n <= S->nseq; n++) { NilOlist(S->pos[n]); free(S->type[n]); }
   free(S->type); free(S->pos);
   for(t = 1; t <= S->ntyp; t++) {
   	for(n = 1; n <= S->nseq; n++){
		free(S->pos_prob[t][n]); free(S->site_pos[t][n]);
	}
   	free(S->pos_prob[t]);
	free(S->site_pos[t]);
	free(S->nsites[t]);
   }
   free(S->site_pos); free(S->pos_prob); free(S->nsites);
   free(S);
}

Boolean	NRandomSites(long t, long N, long max, st_type S)
{
	long	n,try, s;
	dh_type	H;

        H=dheap(S->nseq+2,3);
        for(try=0, s=0; s< N; ){
              for(n =1; n <= S->nseq; n++) insrtHeap(n,(keytyp)Random(),H);
              while((n =delminHeap(H)) != 0){
                if(AddRandomSite(t,n,100,S) == 0 && try++ > max){
			return FALSE;
                }
                if(++s >= N) break;
              }
        }
        Nildheap(H);
	return TRUE;
}

long	AddRandomSite(long t,long n, long ntries, st_type S)
/* Add and return a random available type t site in sequence n */
{
	long	newsite,seqlen,iter=0;
	double	r,ran;

	if(S->len_seq[n] < S->len_elem[t]) return 0;
	seqlen = S->len_seq[n] - S->len_elem[t] + 1;
	do {
	   do{
		  ran = (double) Random();
		  r = ran/(double) LONG_MAX;
       		  newsite=(long)(r*(double) seqlen)+1;
	   } while(newsite > seqlen);
	   if(iter++ > ntries) return 0;
      	} while(OccupiedSite(t,n,newsite,S));
	AddSite(t,n,newsite,S);
	return newsite;
}

void	ShiftSitesM(st_type S, long t, long d)
/* shift all type t sites d spaces to the left(?). */
{
	long	i,c;
	Boolean	left;

	if(d < 0){ c = '-'; d *= -1; left = FALSE;}
	else if(d > 0){ c = '+'; left = TRUE; }
	else return;
	fprintf(stderr,"[");
	for(i=1; i <=d; i++) {
	  fprintf(stderr,"%c", (char) c);
		ShiftSites(S, t, left);
	}
	fprintf(stderr,"]");
}

void	ShiftSites(st_type S, long t, Boolean left)
/*********************** shift left *******************
 Shift all type t sites one position to the left or right.
 	WARNING: assumes that new positions are available.
 SHIFT LEFT:

     site 		  type[n][site] = VACANT ('o')
      |			  type[n][site + 1] = t; ('A')
      A   a   a   a   ...   a   o       	<- sequence n
    |+0 |+1 |+2 |+3 | ... |w-1| w |
      o   A   a   a   ...   a   a       	
          |		  type[n][site + w] = BLOCKED (t)
        site + 1
	
 SHIFT RIGHT:

         site 		  type[n][site] = BLOCKED (t)
          |		  type[n][site - 1] = t; ('A')
      o   A   a   a   ...   a   a       	<- sequence n
    |-1 |+0 |+1 |+2 | ... |w-2|w-1|
      A   a   a   a   ...   a   o       	
      |			  type[n][site + w - 1] = VACANT ('o')
     site - 1
*******************************************************/
{
	long	n,k,site;

    if(left){		/* free 1; 1..w-1 = 2..w; w = new */
/**** fprintf(stderr,"\nmodel %d ********\n", t);
	PutTypeSites(stderr, S); PutSites(stderr, t, S, NULL,NULL);
/****/
	for(n=1; n<= S->nseq; n++){
	 if(S->nsites[t][n] > 0){
           GetOlist(S->tmp,S->pos[n]);
	   for(k=1;(site=S->tmp[k]) != 0; k++){
	   	site=S->tmp[k];
		if(S->type[n][site]==t){
		   RmOlist(site,S->pos[n]);	   /* move pattern left */
		   InsertOlist(site +1, S->pos[n]);
		}
	   }
	   bubble_sort_sites(S->site_pos[t][n]);
	   for(k=S->nsites[t][n]; k > 0; k--){
	   	site = S->site_pos[t][n][k];
		S->type[n][site] = VACANT; 
		S->type[n][site+1] = t;	 	/* move over +1 */
		if(S->type[n][site+S->len_elem[t]] != VACANT)
			sites_error("shift operation is blocked.");
		S->type[n][site+S->len_elem[t]] = BLOCKED(t);
		S->site_pos[t][n][k]++; /** move all pos. ahead one **/
	   }
	 }
	}
/**** PutTypeSites(stderr, S); PutSites(stderr, t, S, NULL,NULL); /****/
    } else {		/* free w; w..2 = w-1..1; 1 = new */
	for(n=1; n<= S->nseq; n++){
	 if(S->nsites[t][n] > 0){
           GetOlist(S->tmp,S->pos[n]);
	   for(k=1;(site=S->tmp[k]) != 0; k++){
		if(S->type[n][site]==t){
		   RmOlist(site,S->pos[n]);	   /* move pattern left */
		   InsertOlist(site-1, S->pos[n]);
		}
	   }
	   bubble_sort_sites(S->site_pos[t][n]);
	   for(k=1; k <= S->nsites[t][n]; k++){
	   	site = S->site_pos[t][n][k];
		if(S->type[n][site-1] != VACANT)
			sites_error("shift operation is blocked.");
		S->type[n][site-1] = t;	   /* move over -1 */
		S->type[n][site] = BLOCKED(t);
		S->type[n][site+S->len_elem[t]-1] = VACANT;
		S->site_pos[t][n][k]--;  /** move all pos. back one **/
	   }
	 }
	}
    }
}

void    GrowSites(long t, st_type S)
/*********************** grow right *************************
 Lengthen all tyep t sites one position to the right.
	 WARNING: assumes that new positions are available.
 GROW RIGHT:
     site 
      |
      A   a   a   a   ...   a   o    	<- sequence n
    |+0 |+1 |+2 |+3 | ... |w-1|	w |	 len_elem[t]++;
      A   a   a   a   ...   a   a       	
           		 type[n][site + w] = BLOCKED(t)
***************************************************************/
{
	long	w,n,k,site;

	w = S->len_elem[t]++;		/* w = old length */
	for(n=1; n<= S->nseq; n++){
           GetOlist(S->tmp,S->pos[n]);
	   for(k=1;(site=S->tmp[k]) != 0; k++){
		if(S->type[n][site]==t){
			if(S->type[n][site + w] != VACANT)
				sites_error("grow operation is blocked.");
			S->type[n][site + w] = BLOCKED(t);
		}
	   }
	}
	for(S->maxinc=w+1,t=1; t<=S->ntyp; t++) 
		S->maxinc = MIN(long, S->maxinc,(S->len_elem[t]-1));
}

void    ShrinkSites(long t, st_type S)
/*********************** shift left *******************
 Shortens all type t sites one position on the right.
 SHRINK RIGHT:
     site 
      |
      A   a   a   a   ...   a   a      	<- sequence n
    |+0 |+1 |+2 |+3 | ... |w-2|w-1|		 len_elem[t]--;
      A   a   a   a   ...   a   o           	
           		 type[n][site + w - 1] = VACANT 

*******************************************************/
{
	long	w,n,k,site;

	w = S->len_elem[t]--;		/* w = old length */
	if(w<3) sites_error("cannot shrink element to length < 3");
	for(n=1; n<= S->nseq; n++){
           GetOlist(S->tmp,S->pos[n]);
	   for(k=1;(site=S->tmp[k]) != 0; k++){
		if(S->type[n][site]==t) S->type[n][site + w-1] = VACANT;
	   }
	}
	S->maxinc = MIN(long, S->maxinc,(S->len_elem[t]-1));
}

void	VacateSite(long t, long n, long site, st_type S)
/*********************** vacate site *******************
 Remove site in sequence n of type t by vacating all positions. 
     site 
      |
      A   a   a   a   ...   a   a      	<- sequence n
    |+0 |+1 |+2 |+3 | ... |w-2|w-1|	
      o   o   o   o   ...   o   o           	
        for p = site ... site + w - 1 -> type[n][p] = VACANT 
*******************************************************/
{
	long	i,end;
	unsigned short  *site_pos=S->site_pos[t][n];

	if(S->type[n][site] != t){
		PutTypeSites(stderr, S);
		PutSites(stderr, t, S, NULL,NULL);
		fprintf(stderr,"ELEMENT %c; seq %ld; site %ld\n",
			(char) ('A' +t -1), n,site);
		sites_error("attempt to remove site where none exists.");
	}
	RmOlist(site,S->pos[n]);
	end = site + S->len_elem[t] -1; 
	for(i=site; i <= end ; i++) S->type[n][i] = VACANT; 
/****************************************************************
            i= 1   2   3   4    (nsites = 4)       
  remove(35, [12, 35, 78, 104, NULL]) 

	-> [12, 104, 78, 104, NULL] -> [12, 104, 78, NULL ] 

  remove(35, [35, NULL]) 

	-> [35, NULL] -> [NULL]

 ****************************************************************/
	for(i=1; i <= S->nsites[t][n]; i++){
		if(site_pos[i] == site){
			site_pos[i] = site_pos[S->nsites[t][n]];
			site_pos[S->nsites[t][n]] = 0;
			S->nsites[t][n]--;
			S->totsites[t]--;
			return ;
		}
	}
	sites_error("VacateSite( ) - this should not happen");
}

void	AddSite(long t, long n, long site, st_type S)
/*************************** add site ********************************
 Add type t site in sequence n.
     site 
      |
      o   o   o   o   ...   o   o           	
    |+0 |+1 |+2 |+3 | ... |w-2|w-1|	
      A   a   a   a   ...   a   a      	<- sequence n
	type[n][site] = t  ('A')
        for p = site + 1 ... site + w - 1 -> type[n][p] = BLOCKED (t) 
***********************************************************************/
{
	long	p,end;

	if(OccupiedSite(t, n, site, S)){
		PutTypeSites(stderr, S);
		PutSites(stderr, t, S, NULL,NULL);
		fprintf(stderr,
			"ELEMENT %c(length=%ld); seq %ld; site %ld\n",
			(char) ('A' +t -1), S->len_elem[t], n,site);
		sites_error("attempt to add site where one exists.");
	}
	InsertOlist(site, S->pos[n]);
	S->type[n][site] = t;
	end = site + S->len_elem[t] - 1; 
	for(p=site+1; p <= end ; p++) S->type[n][p] = BLOCKED(t); 
	S->nsites[t][n]++;
	S->totsites[t]++;
	S->site_pos[t][n][S->nsites[t][n]] = (unsigned short) site;
}

Boolean OccupiedSite(register long t, register long n, register long site, 
	register st_type S)
/* determine if a site is blocked. */
{
	register long	p,end;

	end = site + S->len_elem[t] - 1; 
	for(p=site; p<end ; p+=S->maxinc){if(S->type[n][p])return TRUE;}
	if(S->type[n][end]) return TRUE;
	return FALSE;
}

long     ChooseSite(long t, long n, st_type S)
/* sample a t site in sequence n of S. return site location 
   WARNING: Assumes that BLOCKED sites have zero probability. */
{
        double  rand_no, cum_prob;
        long     site,seqlen;

	seqlen = S->len_seq[n] - S->len_elem[t] + 1;
	do {
        	rand_no = (double) Random()/(double) LONG_MAX;
	} while(rand_no == 0.0);
        rand_no *= S->pos_prob[t][n][0];
        for(site=1,cum_prob = 0.0; site <= seqlen; site++){
           if((cum_prob += (double) S->pos_prob[t][n][site]) >= rand_no){
        	AddSite(t,n,site,S);
        	return site;
	   }
	}
	PutTypeSites(stderr, S);
	PutSites(stderr, t, S, NULL,NULL);
	fprintf(stderr,"total prob = %g; rand_no = %g\n",
		S->pos_prob[t][n][0],rand_no);
	sites_error("ChooseSite( ) - this should not happen!?");
}

double	MissInfoSites(long typ, st_type S)
/* Given a 2xarray of probabilities for sites occurs at each position 
   in each sequence return the missing position information in bits */
{
	double	term2,term3,d,**zeta,lambda;
	long	s,t,end;
	double  **pos_prob=S->pos_prob[typ];

	NEWP(zeta, S->nseq+1,double);
	for(s=1; s<=S->nseq; s++) {
	   NEW(zeta[s], S->len_seq[s],double);
	   end = S->len_seq[s] - SiteLen(typ,S) + 1;
	   for(d=0.0, t=1; t<= end; t++) d += pos_prob[s][t];
	   for(t=1; t<= end; t++)  zeta[s][t] = pos_prob[s][t]/d;
	}
	for(term2=0.0,s=1; s<=S->nseq; s++) {
		end = S->len_seq[s] - SiteLen(typ,S) + 1;
		lambda = 1.0/(double) end;
		for(t=1; t<= end; t++) 
			term2 += zeta[s][t]*log(lambda);
	}
	for(term3=0.0,s=1; s<=S->nseq; s++) {
	   end = S->len_seq[s] - SiteLen(typ,S) + 1;
	   for(t=1; t<= end; t++) 
		if(zeta[s][t]>0.0) term3 += zeta[s][t]*log(zeta[s][t]);
	}
	for(s=1; s<=S->nseq; s++) free(zeta[s]);
	free(zeta);
	return (-1.4427*(term2 - term3));
}

void	OrderSites(long n, long *order, st_type S)
/* modifies an array to give the order of the types of sites in
   sequence n.  WARNING: array is assumed to be long enough to 
   hold all sites in seq n. */
{
	long	s;
	GetOlist(S->tmp, S->pos[n]); 
	for(s=1; S->tmp[s] != 0; s++) order[s] = S->type[n][S->tmp[s]];
	order[s] = 0;
}

void	PosTSites(long t, long n, long *pos, st_type S)
/* modifies array pos to contain the positions of type t sites in seq. n */
{
	long	s,site,i;
	GetOlist(S->tmp, S->pos[n]); 
	for(i=s=1; (site=S->tmp[s]) != 0; s++) {
		if(S->type[n][site]==t) pos[i++] = S->tmp[s];
	}
	pos[i] = 0;
}

void	PosSites(long n, long *pos, st_type S)
/* modifies array pos to contain the positions of sites in sequence n */
{
	long	s;
	GetOlist(S->tmp, S->pos[n]); 
	for(s=1; S->tmp[s] != 0; s++) pos[s] = S->tmp[s];
	pos[s] = 0;
}

/********************************** output ******************************/

void    PutSitesF(FILE *fptr,long t,st_type S,double **site_prob, 
	Boolean *off)
/** PutSites to create a figure for publication **/
{
        ss_type P=S->data;
        long    i,n,s,leng,e,end,N,start,length;
        Boolean are_sites = FALSE;
        e_type  E;
        char    r,c;
        /** long        flank=5; /****/
        /** long        flank=10; /****/
        /**/ long       flank=0; /****/

        fprintf(fptr,"\n\n");
        length = SiteLen(t,S);
        for(N=0, n=1; n<= NSeqsSeqSet(P); n++){
           bubble_sort_sites(S->site_pos[t][n]);
           E=SeqSetE(n,P);
           if(nSites(t,n,S) > 0) N++;
           for(s=1; s<= nSites(t,n,S); s++){
                are_sites = TRUE;
                /*** fprintf(fptr,"%2d-%-2d ",n,s);/***/
                start= S->site_pos[t][n][s];
                fprintf(fptr,"%4ld  ",start);
                e = start + length - 1;
                end = e + flank;
                for(i=start-flank; i <= end; i++){
                   if(i < 1 || i > (long) LenSeq(E)) fprintf(fptr," ");
                   else{
                        r = ResSeq(i,E);
                        if(OpenPos(n,i,S)) c = AlphaCharLow(r,SeqSetA(P));
                        else c = AlphaChar(r,SeqSetA(P));
                        if(i == e) fprintf(fptr,"%c ", c);
                        else if(i == start) fprintf(fptr," %c", c);
                        else fprintf(fptr,"%c", c);
                   }
                }
                fprintf(fptr," %4ld",e);
                if(site_prob != NULL) {
                        fprintf(fptr," (%0.2f) ",site_prob[n][start]);
                        /***/ PutSeqID(fptr,E); /*** TEMP ***/
                }
                fprintf(fptr,"\n");
           }
        }
        if(are_sites){
           if(off != NULL){
                fprintf(fptr,"sites:%17c",' ');
                for(s=1; s <= S->len_elem[t]; s++){
                        if(off[s])fprintf(fptr," ");
                        else fprintf(fptr,"*");
                 }
           }
           fprintf(fptr,"\n%*s", 23, "");
           for(i=1; i<= SiteLen(t,S); i++) {
                if(i%5==0) fprintf(fptr,"%5ld",i);
           }
        }
        fprintf(fptr,"\n\t(%ld sites in %ld sequences)\n\n",S->totsites[t],N);
}

void	PutSites(FILE *fptr,long t,st_type S,double **site_prob, Boolean *off)
{
	ss_type	P=S->data;
	long	i,n,s,leng,e,end,N,start,length;
	Boolean	are_sites = FALSE;
	e_type	E;
	char	r,c;
   
	fprintf(fptr,"\n\n");
	length = SiteLen(t,S);
	for(N=0, n=1; n<= NSeqsSeqSet(P); n++){
	   bubble_sort_sites(S->site_pos[t][n]);
	   E=SeqSetE(n,P);
	   if(nSites(t,n,S) > 0) N++; 
	   for(s=1; s<= nSites(t,n,S); s++){
		are_sites = TRUE;
	   	fprintf(fptr,"%2ld-%-2ld ",n,s);
		start= S->site_pos[t][n][s]; 
		fprintf(fptr,"%4ld  ",start);
		e = start + length - 1;
		end = e + 10;
		for(i=start-10; i <= end; i++){
		   if(i < 1 || i > (long) LenSeq(E)) fprintf(fptr," ");
		   else{
			r = ResSeq(i,E);
			if(OpenPos(n,i,S)) c = AlphaCharLow(r,SeqSetA(P));
			else c = AlphaChar(r,SeqSetA(P));
			if(i == e) fprintf(fptr,"%c ", c);
			else if(i == start) fprintf(fptr," %c", c);
			else fprintf(fptr,"%c", c);
		   }
		}
		fprintf(fptr," %4ld",e);
		if(site_prob != NULL) {
			fprintf(fptr," (%0.2f)",site_prob[n][start]);
			/** PutSeqID(fptr,E); /*** TEMP ***/
		} 
		fprintf(fptr,"\n");
	   }
	}
	if(are_sites){
	   if(off != NULL){
	   	fprintf(fptr,"sites:%17c",' ');
	   	for(s=1; s <= S->len_elem[t]; s++){
			if(off[s])fprintf(fptr," "); 
			else fprintf(fptr,"*"); 
	  	 }
	   }
	   fprintf(fptr,"\n%*s", 23, "");
	   for(i=1; i<= SiteLen(t,S); i++) {
		if(i%5==0) fprintf(fptr,"%5ld",i);
	   }
	}
	fprintf(fptr,"\n\t(%ld sites in %ld sequences)\n\n",S->totsites[t],N);
}

void    PutTypeSites(FILE *fptr, st_type S)
/*  	e.g.    3: D(24)-A(35)-C(45)  */
{
	long i,n,t,tot;

	fprintf(fptr,"\n");
	for(n=1; n<=S->nseq; n++){
	   for(tot=0,t=1; t<=S->ntyp; t++) tot += S->nsites[t][n];
	   if(tot > 0){
	      fprintf(fptr,"%3ld: ",n);
	      for(i=1; i<=S->len_seq[n]; i++){
		if(S->type[n][i] > 0){
		   t = S->type[n][i];
	   	   if(tot == 0) fprintf(fptr, "%c %ld", '-',n); tot = 0;
		   if(S->ntyp <= 26){
		     fprintf(fptr,"%c(%ld)",(char) ('A'+t-1),i);
		   } else { fprintf(fptr,"%ld(%ld)",t,i); }
		}
	      }
	      fprintf(fptr,"\n");
	   }
	}
	fprintf(fptr,"\n");
}

/*** output for creating an entry in the motif database ***/

void    PutSitesMtfDBS(FILE *fptr,long t,st_type S, double **prob,
	Boolean *off)
/** PutSites to create a figure for publication **/
{
        ss_type P=S->data;
        long    i,n,s,leng,e,end,N,start,length;
        e_type  E;
        char    r,c;

        fprintf(fptr,"AL   ");
        if(off != NULL){
                for(s=1; s <= S->len_elem[t]; s++){
                        if(off[s])fprintf(fptr,".");
                        else fprintf(fptr,"*");
                 }
        } else for(s=1; s <= S->len_elem[t]; s++) fprintf(fptr,"*");
        fprintf(fptr,"\n");
        length = SiteLen(t,S);
        for(N=0, n=1; n<= NSeqsSeqSet(P); n++){
           bubble_sort_sites(S->site_pos[t][n]);
           E=SeqSetE(n,P);
           if(nSites(t,n,S) > 0) N++;
           for(s=1; s<= nSites(t,n,S); s++){
        	fprintf(fptr,"     ");
                start= S->site_pos[t][n][s];
                e = start + length - 1;
                end = e;
                for(i=start; i <= end; i++){
                   if(i < 1 || i > (long) LenSeq(E)) fprintf(fptr," ");
                   else{
                        r = ResSeq(i,E);
                        c = AlphaChar(r,SeqSetA(P));
                        fprintf(fptr,"%c", c);
                   }
                }
                fprintf(fptr," (%-ld ",start);
                PutSeqID(fptr,E); 
                if(prob != NULL) {
                        fprintf(fptr," %.2f",prob[n][start]);
                }
                fprintf(fptr,")\n");
           }
        }
}

void	PutScanSites(FILE *fptr, long t, st_type S, Boolean *off)
{
	ss_type	P=S->data;
	long	i,r,n,s,leng,end,length=SiteLen(t,S);
	Boolean	flag;
	e_type	E;
   
	fprintf(fptr,"\n\n");
	for(flag=TRUE, n=1; n<= NSeqsSeqSet(P); n++){
	   bubble_sort_sites(S->site_pos[t][n]);
	   E=SeqSetE(n,P);
	   for(s=1; s<= nSites(t,n,S); s++){
		if(flag){
	   	   if(off != NULL){
	   		for(i=1; i<=length; i++){
			   if(off[i])fprintf(fptr,"."); 
			   else fprintf(fptr,"*"); 
	  	 	}
			fprintf(fptr,"\n");
	   	   } else for(i=1; i<=length; i++) fprintf(fptr,"*"); 
		   flag = FALSE;
		}
		end = S->site_pos[t][n][s] + length - 1;
		for(i=S->site_pos[t][n][s]; i <= end; i++){
		   if(i < 1 || i > (long) LenSeq(E)) fprintf(fptr," ");
		   else{
			r = ResSeq(i,E);
			fprintf(fptr,"%c", AlphaChar(r,SeqSetA(P)));
		   }
		}
		fprintf(fptr,"\n");
	   }
	}
	fprintf(fptr,"\n");
}

/********************************** private ******************************/
long	print_sites(unsigned short *L)
{
	long	i;
	L++;
	fprintf(stderr,"\n");
	for(i=0; L[i]!=0; i++) fprintf(stderr," %d",L[i]);
	fprintf(stderr," (null)\n");
}

long     bubble_sort_sites(unsigned short *L)
{
        unsigned short i,j,t,n;
        unsigned short X[2000];
	long temp=0;

	L++;
        for(n=i=0; L[i]!=0; i++) n++;
/** DEBUG **
	for(i=0; L[i]!=NULL; i++) { X[i]=L[i]; } X[i]=NULL;
/** DEBUG **/
        for(i=1; i < n; i++){
                for(j=i;j > 0 && (L[j-1] > L[j]);j--){
                        t=L[j]; L[j]=L[j-1]; L[j-1]=t;
                }
        }
/** DEBUG **
	for(i=1; L[i]!=NULL; i++) { if(L[i]==L[i-1]) temp++; }
if(temp){
	fprintf(stderr,"\n");
	for(i=0; X[i]!=NULL; i++) fprintf(stderr," %d",X[i]);
	fprintf(stderr,"\n");
	for(i=0; L[i]!=NULL; i++) fprintf(stderr," %d",L[i]);
	fprintf(stderr,"\n\n");
}
/** DEBUG **/
	return temp;
}

void	sites_error(char *s){fprintf(stderr,"sites error: %s\n",s);exit(1);}

