/* smodel.c 
 *
 * last modified: 11-19-03 mjp
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/smodel.c,v 2.1 2004/06/07 19:03:43 thompson Exp $
 *
 *
 * 10-23-03 - added SetSModelFromProbCounts() which is needed to support
 *	      probability matrices being specified to specify a mtf/model
 *	      in the snfile. formerly, only a list of sites could be 
 *	      specified.
 * 
 */

#include "smodel.h"

sm_type	SModel(long length, double pseudo, long *counts, a_type A)
{ return MkSModel(NULL,length, pseudo, counts, A); }


/***************************************************************************
 * MkSModel()
 *
 * mjp: null   == boolean array; T for OFF positions in model
	      	  and F for ON positions. 
        length == site (model) length (or really width, ie,
         	  number of bps).
	pseudo == by default, == 0.1
	counts == nuc type counts in database file, [1..4]
       	A      == alphabet datatype; when called from ReadScan, it's F->A 
 *
 */


sm_type	MkSModel(Boolean *null, long length, double pseudo, long *counts, a_type A)
/* create and return a smodel with totsites = 0 and unfragmented columns. */
{
	sm_type	M;
	long j,r;
	
	NEW(M,1,smodel_type);
	M->A = A;
	M->length = length;
	M->method = 'g';  /** Gribskov's method **/
	NEW(M->null, length+2, Boolean);
	if(null != NULL){
	   for(j=1; j<=length; j++) {
		if(null[j]) M->null[j] = TRUE; else M->null[j]=FALSE;
	   }
	} else { for(j=1; j<=length; j++) M->null[j] = FALSE; }
	M->pseudo = MAX(double,pseudo,0.00000001); M->totsites = 0;
	NEW(M->freq,nAlpha(A)+2,double);
	NEW(M->counts,nAlpha(A)+2,long);

	/* in the model, store the counts of each alphabet letter that are
	 * present in the *database*. */
	for(M->tot_cnts=0.0, r=0; r<= nAlpha(A); r++){
		M->counts[r]= counts[r];
		if(r) M->tot_cnts += counts[r];
	} 
	NEW(M->N0,nAlpha(M->A)+2,double);
	NEW(M->temp,nAlpha(M->A)+2,double);
	/* record freq of each type in *database*. */
	for(r=1; r<= nAlpha(A); r++) {
		M->freq[r] = (double) counts[r]/(double) M->tot_cnts;
	}
	NEWP(M->site_freq, M->length+2, long);
	NEWP(M->likelihood, M->length+2, double);
	NEW(M->tmp_val, M->length+2, double);
	for(j=0; j <= M->length; j++){	
	    NEW(M->site_freq[j],nAlpha(M->A) +2, long);
	    NEW(M->likelihood[j], nAlpha(M->A)+2, double);
	}
	M->smx = NULL;
	M->update = TRUE;
	return M;
}

void	InitSModel(sm_type M)
/* Initialize model M to no sites */
{
	long	j,b;

	M->totsites = 0;
	for( j = 1; j <= M->length; j++){	
	   for(b=0; b <= nAlpha(M->A); b++){ M->site_freq[j][b] = 0; }
	}
	M->update = TRUE;
}

sm_type NilSModel(sm_type M)
/* Destroy model M */
{
   long j;
   if(M!=NULL){ 
	free(M->freq); free(M->counts); free(M->N0); free(M->temp);
	for(j = 0; j <= M->length; j++){
	    free(M->site_freq[j]); 
	    free(M->likelihood[j]);
	}
	free(M->site_freq); 
	free(M->likelihood); free(M->tmp_val);
	if(M->smx != NULL) NilSMatrix(M->smx);
	free(M->null);
	free(M);
   }
   return (sm_type) NULL;
}

char	*MaxSegSModel(sm_type M)
{
	if(M->update || M->smx == NULL) get_smx_smodel(M);
	return MaxSegSMatrix(M->smx);
}

smx_typ	GetSMatrixSModel(sm_type M)
/* Return the score matrix for model M; score is in half bits */
{
	if(M->update || M->smx == NULL) get_smx_smodel(M);
	return M->smx;
}

long	ScoreSModel(register char *seq, register long pos, register sm_type M)
/* Return the log likelihood score for site at pos in seq using model M */
/* score is in half bits */
{
	if(M->update || M->smx == NULL) get_smx_smodel(M);
	return ScoreSMatrix(seq,pos,M->smx);
}

double	PvalSModel(register char *seq, register long pos, register sm_type M)
/* Return the p-value for log likelihood score of site at pos in seq using 
   model M score is in half bits. */
{
        register long	score;
	/** printf("method = %c\n",M->method);/***/
	if(M->update || M->smx == NULL) get_smx_smodel(M);
	score = ScoreSMatrix(seq,pos,M->smx);
        return SMatrixProb(score, M->smx);
}

void	get_smx_smodel(register sm_type M)
{
        register long i,score,c,ave,x;
	register double	s,factor; /*factor only needed for regular model*/
	int      debug = 0;

	if(M->update) update_smodel_freqN(M);
	if(M->smx == NULL) {
	   if(M->method == 'r' || M->method == 'a')
		factor=2.8853901/(M->maxscore/SMODEL_MAX_SCORE);
	   /*** fprintf(stderr,"DEBUG: maxscore = %g; factor = %g\n",
		M->maxscore, factor); /****/
	   M->smx = MkSMatrix(2.0,M->length,M->freq,M->A);
           for(i=1; i<=M->length; i++){
             if(!M->null[i]){
               for(ave=0,c=1; c<=nAlpha(M->A); c++){
		if(M->N0[c] != 0){
		 switch (M->method) {
		  case 'r': /****** Regular Model *******/
		  case 'a': /****** asymptotic Model *******/
		    score=(long)floor(factor*log(M->likelihood[i][c])+0.5); 
		    /* score=(long)floor(100.0*log(M->likelihood[i][c])); */  /* BT 05/27/04 */ /* temp for testing */
		    break;
		  case 'c': /****** CHIP'S METHOD??? *******/
		    for(s = 0.0,x = 1; x<=nAlpha(M->A); x++){
			s += ((double)M->site_freq[i][x]/(double)M->totsites)
					* (double) blosum62L[c][x];
		    }
                    score = (long) floor((10.0*log(s) + 0.5));
		    break;

		  case 'm': /****** PAM1 DNA METHOD??? *******/
		    for(s = 0.0,x = 1; x<=nAlpha(M->A); x++){
		      s += ((double)M->site_freq[i][x]/(double)M->totsites)
			* pam1ScoreL[c][x]; 
		    }
		    score = (long) floor((10.0*log(s) + 0.5)); 
		    break;
	
		  case 'i': /****** IDENTITY METHOD??? *******/
		    for(s = 0.0,x = 1; x<=nAlpha(M->A); x++){
		      s += ((double)M->site_freq[i][x]/(double)M->totsites)
			* (double) iden[c][x]; 
		    }
		    score = (long) floor((10.0*log(s) + 0.5)); 
		    break;
	
		  default: /**** GRIBSKOV'S PROFILE METHOD ****/
		    for(s = 0.0,x = 1; x<=nAlpha(M->A); x++){
			  s += (double)M->site_freq[i][x]*blosum62[c][x];
		    }
                    score = (long) floor(4.0*s + 0.5);
		 }
		 /**	 printf( "Debug: i,c,score %d %d %ld\n", i, c, score );/***/  /* DEBUG */

		 SetSMatrix(c,i,score,M->smx);
		 ave += score;
		}
	       }
	       SetSMatrix(0,i,0,M->smx);
	       /** SetSMatrix(0,i,(ave/nAlpha(M->A)),M->smx);/****/
	     }
           }
	} else smodel_error("second update not yet implemented");

	if( debug )
	  {
	    for( i = 1; i <= M->length; i++ )
	      {
		printf( "row %ld ", i );
		for( x = 1; x <=nAlpha(M->A); x++ )
		  {
		    printf( "%ld ", M->smx->score[i][x] );
		  }
		printf( "\n" );
	      }
	  }
}

double	**RealScoresSModel(sm_type M)
/** return a 2 dimensional array of real values log-odds scores **/
{
        long	i,c,x;
	double	s,factor; /*factor only needed for regular model*/
	double	ave,score,**scores;

	if(M->update) update_smodel_freqN(M);
	/** printf("method = %c\n",M->method);/***/
	if(M->method == 'r' || M->method == 'a') {
		factor=2.8853901/(M->maxscore/SMODEL_MAX_SCORE);
	}
	NEWP(scores,nAlpha(M->A)+2,double);
        for(c=0; c<=nAlpha(M->A); c++) NEW(scores[c],M->length+2,double);
        for(i=1; i<=M->length; i++){
             if(!M->null[i]){
               for(ave=0,c=1; c<=nAlpha(M->A); c++){
		if(M->N0[c] != 0){
		 switch (M->method) {
		  case 'r': /****** Regular Model *******/
		  case 'a': /****** asymptotic Model *******/
		     score=factor*log(M->likelihood[i][c]);
		    break;
		  case 'c': /****** CHIP'S METHOD??? *******/
		    for(s = 0.0,x = 1; x<=nAlpha(M->A); x++){
			s += ((double)M->site_freq[i][x]/(double)M->totsites)
					* (double) blosum62L[c][x];
		    }
                    score = log(s);
		    break;
		  default: /**** GRIBSKOV'S PROFILE METHOD ****/
		    for(s = 0.0,x = 1; x<=nAlpha(M->A); x++){
			  s += (double)M->site_freq[i][x]*blosum62[c][x];
		    }
                    score = s;
		 }
		 scores[c][i] = score;
		 ave += score;
		} else scores[c][i] = 0.0;
	       }
	       scores[0][i] = 0.0; /****/
	       /*** scores[0][i] = (ave/(double) nAlpha(M->A)); /****/
	     } else for(c=0; c<=nAlpha(M->A); c++) scores[c][i] = 0.0;
        }
	return scores;
}

long	DrawSeqSModel(char *seq, sm_type M)
/** draw and return a random sequence from the model **/
{
	long	c,j,b;
	double	r,total;

	if(M->update) update_smodel_freqN(M);
        for(j=1; j<=M->length; j++){
         if(!M->null[j]){
	   for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
	      total += M->site_freq[j][b] + M->N0[b];
	   }
	   r = (double) Random()/(double) LONG_MAX; /** 0 <= r <= 1 **/
	   for(c=1; c <= nAlpha(M->A); c++){
		r -= ((double) M->site_freq[j][c]+M->N0[c])/total;
		if(r <= 0.0){ seq[j] = c; break; }
	   }
	 } else seq[j] = 0;
        }
	return M->length;
}

long	CellScoreSModel(long r, long pos, sm_type M)
/* r = residue; pos = position */
{
	if(M->update || M->smx == NULL) get_smx_smodel(M);
        if(pos > 0 && pos <= M->length && r >= 0 && r <= nAlpha(M->A)){
		return ValSMatrix(pos,r,M->smx); 
	} else smodel_error("Referenced cell is out of bounds");
}

double	RelProbSModel(register char *seq, register long pos, register sm_type M)
/* Return the relative probability of site at pos in seq being in model M */
{
        register long j;
        register double p;
	if(M->update) update_smodel_freqN(M);
        for(p=1.0,j=1; j<=M->length; j++,pos++){
                if(!M->null[j]) p*=M->likelihood[j][(seq[pos])];
        }
        return p;
}

double	ProbSModel(register char *seq, register long pos, register double p,
	register sm_type M)
/* Return the relative probability of site at pos in seq being in model M */
{
        register long j;
        register double L;
	if(M->update) update_smodel_freqN(M);
        for(L=1.0,j=1; j<=M->length; j++,pos++){
                if(!M->null[j]) L*=M->likelihood[j][(seq[pos])];
        }
	return (L*p/((1.0 - p) + L*p));
}

void    RmSModel(char *seq, long site, sm_type M)
/* remove the segment at site in seq from model */
{
        long j;
        for(j=1; j<=M->length; j++,site++){
	    M->site_freq[j][seq[site]]--;
        }
	M->totsites--; M->update = TRUE;
}

void	Add2SModel(char *seq, long site, sm_type M)
/* Add the segment at site in seq to model */
{
	long j;
	for(j=1; j<=M->length; j++,site++){
	    M->site_freq[j][seq[site]]++;
	}
	M->totsites++; M->update = TRUE;
}

void	update_smodel_freqN(register sm_type M)
/* Normalize frequencies to avoid overflow */
{
	register long	j,b,score,max;
	register double factor;
	register a_type	A=M->A;

	if(M->smx != NULL) { NilSMatrix(M->smx); M->smx = NULL; }
	/*** DETERMINE REGULAR AND NORMALIZED NONSITE FREQUENCIES ***/
	M->npseudo = M->pseudo * sqrt((double)M->totsites);
	for(b=1; b <= nAlpha(A); b++) {
		M->N0[b]= M->npseudo * M->freq[b];
		M->site_freq[0][b] = M->counts[b];
	}
	for(factor=0,b=1; b <= nAlpha(A); b++) {
		factor += M->site_freq[0][b]+M->N0[b];
	}
	for(b=1; b <= nAlpha(A); b++){
	   if(M->N0[b] > 0.0){
             M->likelihood[0][b] = (M->site_freq[0][b]+M->N0[b])/factor;
	   }
	}
	/*** DETERMINE NORMALIZED SITE FREQUENCIES ***/
	factor = (double) (M->totsites + M->npseudo);
        for(j=1; j<=M->length; j++){
	    if(!M->null[j]){
               	M->likelihood[j][0] = 1.0;
        	for(b=1; b <= nAlpha(A); b++){
		  if(M->N0[b]==0.0){
		    M->likelihood[j][b] = 0.0;
		  } else {
               	    M->likelihood[j][b] = (M->site_freq[j][b]+M->N0[b])/factor;
               	    M->likelihood[j][b] /= M->likelihood[0][b];
		    /*	  fprintf(stderr,"DEBUG: likelihood[%d][%c] = %f\n",
			   j,AlphaChar(b,M->A), M->likelihood[j][b]); */ 
		    /* DEBUG */
		  }
           	}
	    }
        }  /* for j <= M->length */
        for(M->maxscore=0.0,j=1; j<=M->length; j++){
             if(!M->null[j]){
                for(max=0,b=1; b<=nAlpha(M->A); b++){
		  if(M->N0[b]!=0.0){
		    score = 
		     (long) floor(2.8853901*log(M->likelihood[j][b])+0.5);
		    max = MAX(long,score,max);
		  }
		}
		M->maxscore += max;
	     }
	}
	M->update = FALSE;
}

void    smodel_error(char *s) { fprintf(stderr,"smodel: %s\n",s); exit(1); }

Boolean	NullSiteSModel(long s,sm_type M)
/* if site s == null in model M return TRUE; else return FALSE. */
{ if(s < 1 || s > M->length) return TRUE; else return M->null[s]; }

double	PutSModel(FILE *fptr, sm_type M)
/* Report the current frequency model. */
{
	long	j,b,r;
	double	total,info,tot_info,p,q;
	Boolean flag = TRUE;

	fprintf(fptr,"PutSModel:\n");
	if(M->update) update_smodel_freqN(M);
	fprintf(fptr,"POS  ");
	for(tot_info=0.0, j=1; j<= M->length; j++){
	  if(!M->null[j]){
	    for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
	      if(M->method == 'a') total += M->site_freq[j][b];
	      else total += M->site_freq[j][b] + M->N0[b];
	      if(flag) fprintf(fptr,"%3c", AlphaChar(b, M->A));
	    }
	    if(flag) { fprintf(fptr,"  Info (bits)\n"); flag = FALSE; }
	    fprintf(fptr,"%4ld ",j);
	    for(info=0.0,b = 1; b <= nAlpha(M->A); b++) {
	    	if(M->method == 'a') p = (M->site_freq[j][b])/total;
	    	else p = (M->site_freq[j][b]+M->N0[b])/total;
		if(p > 0.0){
			q = M->freq[b];
			info += p*log(p/q)/log(2.0);
		}
#if 0
	    	fprintf(fptr,"%3d (%ld/%.4f) ", (long)(100*p+0.5), 
			M->site_freq[j][b],M->N0[b]);
#else
	    	fprintf(fptr,"%3ld", (long)(100*p+0.5));
#endif
	    }
	    fprintf(fptr,"   %1.3lf\n", info);
	    tot_info += info;
	  }
	}
	fprintf(fptr,"non-\nsite ");	
	for(b = 1; b <= nAlpha(M->A); b++){
		r = (long) (100.0 * M->N0[b]/M->npseudo);
	    	fprintf(fptr,"%3ld", r);
	}
	fprintf(fptr,"\n\tinformation = %g\n",tot_info);
	fprintf(fptr,"end PuTSModel.\n");
	return tot_info;
}

double	VarianceInfoSModel(double *average, long N, sm_type M)
/* determine the variance of the information content of model M using
   N independent simulated models. */
{
	long	i,j,k,n,s,c;
	long	*counts;
	double	d,var,I0,*I,r,ave;
	double	**target,p,*q,total;

	if(M->update) update_smodel_freqN(M);
	I0 = InfoSModel(M);
	q = M->freq;
	NEWP(target,M->length+2,double);
	total = (double) M->totsites;
        for(i=1; i<=M->length; i++) {
	     NEW(target[i],nAlpha(M->A)+2,double);
	     for(c=1; c <= nAlpha(M->A); c++){
		target[i][c] = (double) M->site_freq[i][c]/total;
	     }
	}
	NEW(I,N+2,double);
        NEW(counts,nAlpha(M->A)+2,long);
	for(ave=var=0.0, n=1; n <= N; n++) {
           for(I[n]=0.0, i=1; i<=M->length; i++) {
	      if(!M->null[i]){
                for(c=1; c <= nAlpha(M->A); c++) counts[c] = 0;
		for(s=1; s<=M->totsites; s++) {
		   r = (double) Random()/(double) LONG_MAX; /** 0 <= r <= 1 **/
                   for(c=1; c <= nAlpha(M->A); c++){
			if((r-=target[i][c]) <= 0.0){ counts[c]++; break; }
                   }
		}	
		for(c=1; c<=nAlpha(M->A); c++){
		   p = (double)counts[c]/total;
		   if(p > 0.0) I[n] += p*log(p/q[c]);
/****
		   printf("%5d   %5d  (%g nats)\n", 
				M->site_freq[i][c],counts[c],I[n]); /****/
		}
	      }
	   }
	   I[n] /= log(2.0);
	   printf("\t\tinfo[%ld] = %g\n",n,I[n]);
	   ave += I[n];
	}
	ave /= N;
	for(var=0.0, n=1; n <= N; n++) {
	   var += (I[n]-ave)*(I[n]-ave);
	}
	var /= (N-1);
        for(i=1; i<=M->length; i++) free(target[i]);
	free(target); free(counts); free(I);
	*average = ave;
	/*** printf("\t\taverage information = %g\n",ave); /****/
	return var;
}

double	InfoSModel2(sm_type M)
/* Return the information content of model M (using pseudocounts).*/
{
	long	j,b,r;
	double	total,info,tot_info,p,q;

	if(M->update) update_smodel_freqN(M);
	for(tot_info=0.0, j=1; j<= M->length; j++){
	  if(!M->null[j]){
	    for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
	      total += M->site_freq[j][b] + M->N0[b];
	    }
	    for(info=0.0,b = 1; b <= nAlpha(M->A); b++){
	    	p = (M->site_freq[j][b]+M->N0[b])/total;
		if(p > 0.0){
			q = M->freq[b];
			info += p*log(p/q);
		}
	    }
	    tot_info += info;
	  }
	}
	return tot_info/log(2.0);
}

double	InfoSiteSModel(long site, sm_type M)
/** return the information content of site in M **/
{
	long	b;
	double	info,p,q;
	
	if(site < 1 || site > M->length) return 0.0;
	if(M->null[site]) return 0.0;
	for(info=0, b = 1; b <= nAlpha(M->A); b++){
	    p = (M->site_freq[site][b])/(double)M->totsites;
	    if(p > 0.0){
		q = M->freq[b];
		info += p*log(p/q);
	    }
	}
	return info/log(2.0);
}

double	InfoSModel(sm_type M)
/* Return the information content of frequency model M. */
{
	long	j,b,r;
	double	info,p,q;

	for(info=0.0, j=1; j<= M->length; j++){
	  if(!M->null[j]){
	    for(b = 1; b <= nAlpha(M->A); b++){
	    	p = (M->site_freq[j][b])/(double)M->totsites;
		if(p > 0.0){
			q = M->freq[b];
			info += p*log(p/q);
		}
	    }
	  }
	}
	return info/log(2.0);
}

/***************************************************************************
 * SetSModelFromProbCounts()
 *
 *	it is assumed a matrix of normalized doubles are passed in. it's
 *	assumed the row sums are 1.0 (or close). 
 *	
 * 	for nucleotide data: 
 *	  the column order of modelProbMatrix[][] is ATCG, while dscan wants
 *	  CGAT. the columns will be stored in the proper order here. 
 *
 * 	for protein data data: 
 *	  the column order of modelProbMatrix[][] is AA's is alphabetic order.
 *	  dscan wants: CGASTNDEQKRHWYFVILMP
 *	  the columns will be stored in the proper order here. 
 *
 *	nAlpha() is used to determine if we're operating on nucleotide or
 *	protein data.
 *
 */

void SetSModelFromProbCounts(sm_type M, double **modelProbMatrix,
			     long modelNumSites)
{
  int i, j, *orderIndices,
      verbose = 0,
      dscanNucOrder[4] = {3,4,1,2},
		      /* A C D E F  G H  I  K  L  M  N P  Q R  S T V  W  Y
			 C G A S T  N D  E  Q  K  R  H W  Y F  V I L  M  P */
      dscanAAOrder[]   = { 3,1,7,8,15,2,12,17,10,18,19,6,20,9,11,4,5,16,13,14};


  if ( verbose )
    printf("SetSModelFromProbCounts()\n");

  /* a freq or counts matrix was read in and from that position counts
   * summing to a counts value that is either user specified or computed are  
   * are computed and stored. */
  M->totsites = modelNumSites;
  M->update = TRUE;
  
  if ( nAlpha(M->A) == 4 ) 
    orderIndices = dscanNucOrder;
  else
    orderIndices = dscanAAOrder;

  for ( i = 1; i <= M->length; i++ ) {
    for ( j = 0; j < nAlpha(M->A); j++ ) {
      M->site_freq[i][orderIndices[j]] = (long) 
                        ((modelProbMatrix[i-1][j] + 0.005) * modelNumSites);

      if ( verbose )
	printf(" %ld ",M->site_freq[i][orderIndices[j]]);
    }
    if ( verbose )
      printf("\n");
  } /* for i */

#if 0
  for ( i = 1; i <= M->length; i++ )
    for ( j = 0; j < nAlpha(M->A); j++ )
      M->site_freq[i][dscanNucOrder[j]] = (long) 
                        ((modelProbMatrix[i-1][j] + 0.005) * modelNumSites);
#endif

#if 0
  /* this doesn't take care of order swapping; note 'j' starts at 1. */
    for ( j = 1; j <= nAlpha(M->A); j++ )
      M->site_freq[i][j] = (long)((modelProbMatrix[i-1][j-1] + 0.005) * modelNumSites);
#endif


} /* SetSModelFromProbCounts */
