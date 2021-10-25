/* scan.c - product multinomial scan program
 *
 * last modified:  11-19-03 mjp
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/scan.c,v 2.3 2004/06/07 19:03:41 thompson Exp $
 *
 */


#include "scan.h"

scn_typ MkScan(double pseudo, char *snfile, long *counts, a_type A, long  hpsz,
	       char  method, short nFlag, int nCount, int nNoOverlap,
	       int   adjustN,            /* BT 11/06/2001 */
	       long *freqModelNumSites,  /* mjp, 11-03-03 */
	       int   numNumSitesArgs)   /* mjp, 11-03-03 */
{
	long	i;
	scn_typ F;

	NEW(F,1,scan_type);
	NEW(F->M,MAX_NUM_MODELS,sm_type);
	F->snfile = String(snfile);
	F->neg_mask = FALSE;
	F->seqfile = FALSE; F->repeats=FALSE; F->shuffle=FALSE;
	F->new_way = FALSE;
	F->min_segs = 2;
	F->max_segs = 10;
	F->maxEval = 0.01; F->logmaxEval = log10(0.01);
	F->singlePval = 0.0;
	F->A = A; F->HG = NULL; 
	F->method = method; 
	ReadScan(F,pseudo,F->snfile,counts,freqModelNumSites,numNumSitesArgs);
	NEWP(F->pos,F->N+2,long);
	for(i=1;i<=F->N; i++) { NEW(F->pos[i],hpsz+3,long); }
	NEW(F->E,hpsz+3,e_type);
	for(i=1;i<=hpsz; i++) F->E[i] = NULL;
	F->H = Mheap(hpsz,3);
	F->hpsz = hpsz;
	F->HG = Histogram("-log10(E-values)", -100, 100,1.0);
	F->nCountFlag = nFlag;      /* BT 10/17/97 */
	F->nCount = nCount;
        F->nNoOverlap = nNoOverlap;  /* BT 10/4/2000 */
	F->adjustN = adjustN;        /* BT 11/06/2001 */
	NEW( F->numSites, F->N+2, int);
	return F;
} /* MkScan */

long	SetMaxEvalScan(double maxEval, scn_typ F)
{
	if(maxEval >= 0.0) {
		F->maxEval = maxEval;
		F->logmaxEval = log10(maxEval);
	} else print_error("maximum E-value must be >= 0.0");
}

void	NilScan(scn_typ F)
{
	long	i;
	NilHist(F->HG); NilMheap(F->H);
	for(i=1;i<=F->N; i++) { NilSModel(F->M[i]); }
	free(F->M); free(F->freq); 
	for(i=1;i<=F->hpsz; i++){ 
		if(F->E[i]!=NULL) { NilSeq(F->E[i]); F->E[i] = NULL; }
	}
	if(F->repeats){
		for(i=1;i<=F->hpsz; i++) { free(F->pos[i]); }
	} else {
		for(i=1;i<=F->N; i++) { free(F->pos[i]); }
	}
	free(F->E); free(F->pos); free(F->snfile);
	free(F);
}
/************************************************************************
 * ReadScan()
 *
 * 10,11/03 this routine was modified to accept a freq and/or counts matrix,
 * in addition to what it has always accepted - a list of sites. the tricky 
 * part of accepting a matrix is that internally, dscan maintains a counts
 * matrix. when a probability freq matrix is present a counts matrix needs
 * to be created based on the supplied probabilities, but how many total 
 * counts should be used? when counts matrix is supplied, the number of counts
 * can be figured out, but in the matrix specification, we're not requiring
 * that position sums be equal, so what are the counts to use? 
 * the solution, ie, the implementation is such:
 * <yet to be written>
 */

long	ReadScan(scn_typ F, double pseudo, char *snfile, long *counts,

		 /* these 2 args are optional cmd line args used to specify
		  * the number of sites to use when a freq model/matrix has
		  * been specified in the snfile. an array is uses to support
		  * multiple models, order in freqModelNumSites[] corresponds
		  * to the order the models are specified in - just counting
		  * freq models, though, not lists of sites. */
		 long *freqModelNumSites,   /* mjp, 11-03-03 */
		 int   numNumSitesArgs)     /* mjp, 11-03-03 */
{
	FILE	*fptr;
	long	m,i,j,length,number,max,min,pos;
	char	**seq,c = '\n';
	long	total;
	Boolean	*null,
	         verbose=FALSE; 
	a_type	A = F->A;
	
	int      numFreqModels = 0;
	double   **modelProbMatrix;

	if ( verbose ) {
	  fprintf(stderr,"nAlpha: %ld\n",nAlpha(A));


	  printf("numNumSitesArgs = %d\n",numNumSitesArgs);
	  for ( i = 0; i < numNumSitesArgs; i++ )
	    printf("  %ld",freqModelNumSites[i]);
	  printf("\n");
	}


   NEW(null,MAX_BLOCK_LENGTH+2,Boolean);
   NEW(F->freq,nAlpha(A)+2,double);
   /* 'seq' is used when a list of sites describe the mtf model. this is the
    * historical default. */
   NEWP(seq,MAX_BLOCK_SIZE+2,char);
   for(i=0; i<=MAX_BLOCK_SIZE; i++) NEW(seq[i],MAX_BLOCK_LENGTH+2,char);

   /* mjp: 10-22-03 modelProbMatrix[][] is used when a probability matrix is
    * present to describe the mtf model. */
   NEWP(modelProbMatrix, MAX_BLOCK_LENGTH, double);
   for (i=0; i < MAX_BLOCK_LENGTH; i++) 
     NEW(modelProbMatrix[i],nAlpha(A), double);

   /*   for(total=i=0; i<=nAlpha(A); i++) total += counts[i];
   for(i=0; i<=nAlpha(A); i++) {
	F->freq[i]= (double)counts[i]/(double)total;
	if(verbose) fprintf(stderr,"freq[%d] = %d/%d = %g\n",
		i,total,counts[i],F->freq[i]);
		} */

   /* BT 11/05/2001 Changed so that we only count nucleotides/amino acid, 
    * not N's for DNA or X's */  
   for(total=0,i=1; i<=nAlpha(A); i++) total += counts[i];
   for(i=1; i<=nAlpha(A); i++) {
	F->freq[i]= (double)counts[i]/(double)total;
	if(verbose) fprintf(stderr,"freq[%ld] = %ld/%ld = %g\n",
		i,counts[i],total,F->freq[i]);
   }

   F->total = (double) total;
   if((fptr = fopen(snfile,"r")) == NULL) {
        fprintf(stderr,"Could not open file \"%s\"\n",snfile);
        print_error("File does not exist!");
   }

   for(c=' ',length = 0,F->N=0,m=1;c!=EOF; m++) {
        if(verbose) fprintf(stderr,"\n"); /****/
        for(; isspace(c); c=fgetc(fptr)) ; /** GO TO FIRST CHARACTER **/
	if(c == '*'){			/** fragmentation row **/
		length = 0;
		do {
                   if(c=='*') {
			length++; null[length] = FALSE;
			if(verbose) fprintf(stderr,"*");
		   } else if(c=='.') {
			length++; null[length] = TRUE;
			if(verbose) fprintf(stderr,".");
		   } else if(!isspace(c)) {
                        fprintf(stderr,"illegal character -> %c",c);
                        fprintf(stderr,"\n");
			print_error("fatal error.");
                   }
		   if(length >= MAX_BLOCK_LENGTH) 
                   		print_error("reset MAX_BLOCK_LENGTH");
                } while((c=fgetc(fptr))!='\n' && c!=EOF);
	} else print_error("model file input error(1): expected mask line");

	if(verbose) fprintf(stderr," length = %ld\n",length);
	if(length == 0) print_error("input error");
        for(; isspace(c); c=fgetc(fptr)) ; /** GO TO FIRST CHARACTER **/

	if ( isalpha(c) ) {
	  /* list of sites is present */
	  for(i=1,pos=0;c!=EOF && c!='*';i++,pos=0) {
		do{
                   if(c=='*') { break; }
		   if(isalpha(c)) {
			if(islower(c)) c=toupper(c);
                        if(verbose) fprintf(stderr,"%c",c);/***/
		        if(++pos >= MAX_BLOCK_LENGTH) 
                   		print_error("reset MAX_BLOCK_LENGTH");
			seq[i][pos] = AlphaCode(c,A);
			/* if(verbose) fprintf(stderr,"'%d'",seq[i][pos]); **/

		   } else if(!isspace(c)) {
                        fprintf(stderr,"seq %ld: illegal character -> %c",i,c);
                        fprintf(stderr,"\n");
			print_error("fatal error.");
                   }
                } while((c=fgetc(fptr))!='\n' && c!=EOF);
		if(verbose) {
		    if(c == '\n' && pos > 0) fprintf(stderr," #%ld\n",i);
		}
                if(i >= MAX_BLOCK_SIZE)
                   print_error("aligned segments > MAX_BLOCK_SIZE");
		else if(pos == 0)	i--;
		else if(pos != length) 
                   print_error("input segments have inconsistent lengths.");
	    }

	  /* here, 'i' is the number of sites in the model that was just read
	     in. */
	  i--; number = i; max = min = 0;
	  if (i < 1){
	    if(F->N==0) print_error("zero segments in snfile.");
	    else m--;	/** no segments in this pass **/
	  } 
	  else if (i > 0) {
		F->N++;
		/* mjp: null   == boolean array; T for OFF positions in model
		      		  and F for ON positions. 
		        length == site (model) length (or really width, ie,
		         	  number of bps).
			pseudo == by default, == 0.1
			counts == nuc type counts in database file, [1..4]
			A      == alphabet datatype; namely F->A 
		*/

		/* 'm' is an index into the 'models' array. */
		F->M[m] = MkSModel(null,length,pseudo,counts,A);
		SetMethodSModel(F->method,F->M[m]);
		for(i=1; i<=number; i++) {
		  Add2SModel(seq[i],1,F->M[m]);
		  /*
		    for(j=1; j<=length; j++){
		      printf("%c",AlphaChar(seq[i][j],A));
		      if (i%50 == 49) 
		      printf("\n");
		    }
		    printf("\n\n");
		  */
		}

		if(verbose){
		   fprintf(stderr,"\n\n");
		   PutSModel(stderr,F->M[m]); 
		   fprintf(stderr,"\n\n"); 
		}
	  } /* if (i > 0) */
	}
	else if ( isdigit(c) || c == '.' ) {
	  /* 10-22-03 mjp: probability or count matrix present; each row will
	   * be normalized. rows are model positions and, for DNA data, cols 
	   * come in as A T C G and alphabetic single letter code order
	   * when protein. For nuc 
	   * data, both bmc and Gibbs produce the data is this required 
	   * format. dscan want nucleotide cols as CGAT, so
	   * they'll be re-arranged later. the number of rows must match the
	   * 'length' variable set above when the fragmentation string was 
	   * read. might try to get fancy in the checking and issue warning
	   * messages if row sums differ; not sure yet. 
	   */
	  int      k, l, numOnPos = 0;
	  long     modelNumSites;

	  double   rowSum, dblVal, 
	           avePosSum,
	           sumOfONRows =  0.0, 
	    	   savedRowSum  = -1.0;


	  /* need to push the char just read back onto the stream. then do
	   * the reads. */
	  ungetc(c,fptr);

	  /*fprintf(stderr,"'c' is: %c\n",c); */

#if 0
	  /* old: */
	  /* need to sum each row and divide each position's sum. still 
	   * probably want to set the number of sites to 100 and do this:
	   * (long) ((dblVal + 0.005) * 100.0)
	   */
#endif
	  if ( verbose )
	    fprintf(stderr,"reading matrix. all rows are normalized.\n");

	  for ( l = 0; l < length; l++ )   {
	    rowSum = 0.0;
	    for ( k = 0; k < nAlpha(A); k++ ) {
	      if ( fscanf(fptr,"%lf",&dblVal) != 1 ) {
		fprintf(stderr,
		 "error reading prob matrix: model: %ld; line: %d; file %s\n",
			m,l,snfile);
		print_error("fatal error.");
	      }
	      else {
		rowSum += dblVal;
		/* store value for now. normalization of the values occurs
		 * below, after the complete row has been read in. */
		modelProbMatrix[l][k] = dblVal;

		/* old method before normalizing row values. */
		/* modelProbMatrix[l][k] =(long) ((dblVal + 0.005) * 100.0); */
		if ( verbose ) {
#if 1
		  fprintf(stderr,"%8.4f ",modelProbMatrix[l][k]);
#else
		  fprintf(stderr,"d: %.4lf ",dblVal);
#endif
		}
	      }
	    }    /* for k < nAlpha(A) */ 
	    
	    /* in an effort to issue warning messages, save the 1st row
	     * sum encountered as one to check all the others against. */
	    if ( savedRowSum == -1.0 )
	      savedRowSum = rowSum;

	    /* for curr row, all values have been read in. do some error 
	     * checking on the line and then normalize values. */
	    /* we expect to see '\n'; 1st read past any whitespace. */
	    /* can't use isspace() as a test here cause it will gobble up
	     * the newline; so allow spaces and tabs. */
	    for(c=fgetc(fptr); ((c==' ') || (c=='\t')) && c!= EOF; 
		c=fgetc(fptr) )  
	      ; 

	    if ( feof(fptr) ) {
	      fprintf(stderr,"unexpected eof in %s while reading prob matrix",
		      snfile);
	      fprintf(stderr,"line: %d;\nmtf len is probably bad.\n",l);
	      print_error("fatal error.");
	    }

	    if ( c != '\n' ) {
	      fprintf(stderr,
		    "\nexpected newline in prob matrix line: %d in %s.\n",
		      l,snfile);
	      fprintf(stderr,"found: '%c'; mtf len is probably bad.\n",c);
	      print_error("model file input error: 2");
	    }

	    /* a complete row has been successfully read in. normalize the 
	     * values - even though they might already be normalized. */
	    if ( verbose ) 
	      fprintf(stderr,"\trowSum: %.4f\n",rowSum);

	    if ( rowSum == 0.0 ) {
	      if ( null[l+1] ) {
		/* position is OFF, 0.0 sum is ok, but issue warning. */
		fprintf(stderr,"WARNING: position sum == 0.0; pos = %d\n",l+1);
	      }
	      else {
		/* position is supposed to be ON, yet row sum == 0.0; bad */
		fprintf(stderr,"FATAL ERROR: position sum == 0.0; pos = %d\n",
			l+1);
		print_error("model file input error.");
	      }
	    }
	    else if ( fabs(savedRowSum - rowSum) > 0.1 ) {
		fprintf(stderr,
	    "WARNING: FYI, row sum (%.4f) at pos %d != 1st row sum (%.4f).\n",
			rowSum, l+1, savedRowSum);
	    }

	    if ( ! null[l+1] ) {
	      /* if this position is ON, add in this row sum to the sum of
	       * sums */
	      sumOfONRows += rowSum;
	      ++numOnPos;
	    }

	    /* normalize based on the current row's (ie, position's) sum. */
	    for ( k = 0; k < nAlpha(A); k++ )
	      if ( rowSum == 0.0 ) 
		modelProbMatrix[l][k] = 0.0;
	      else
		modelProbMatrix[l][k] /= rowSum;

	  }  /* l < length */

	  /* to be compatible w/ how the code worked when only a list of
	   * sites could be specified, sync back upto the next model, or
	   * EOF, whatever's there. */
	  for(; isspace(c); c=fgetc(fptr)) ; /** GO TO FIRST CHARACTER **/

	  /* now make a 'model' from the prob matrix just read in. */

	  F->N++;
	  /* mjp: null   == boolean array; T for OFF positions in model
		      		  and F for ON positions. 
		  length == site (model) length (or really width, ie,
		         	  number of bps).
		  pseudo == by default, == 0.1
		  counts == nuc type counts in database file, [1..4]
		  A      == alphabet datatype; namely F->A 
	  */
	  /* 'm' is an index into the 'models' array. */

	  F->M[m] = MkSModel(null,length,pseudo,counts,A);

	  SetMethodSModel(F->method,F->M[m]);

	  /* figure out how many 'counts' or num_sites to use for this model.
	   * if a prob freq matrix was supplied (which is difficult to
	   * determine), hopefully the user suppiled counts to use. if not
	   * we'll used the default value (100). if a counts matrix has been
	   * supplied, also difficult to determine, compute the 'average'
	   * across all rows; but also need to check if user supplied a
	   * counts values on the cmd line. since a user specified counts
	   * takes precedence, approach it from that angle - but try to issue
	   * a warning message if things are way off. 
	   */

	  /* need this in a few places, compute it now. */
	  avePosSum = sumOfONRows / numOnPos;

	  if ( freqModelNumSites == NULL ) {
	    /* printf("HERE I AM, freqModelNumSites == NULL\n"); */
	    /* user did not specify counts for model. if counts matrix was
	     * supplied here, use average counts across all positions as the
	     * number of counts to use, if prob freq matrix was supplied,
	     * use the default num counts (100). 
	     * it's hard to tell if a prob freq matrix or a counts matrix was
	     * supplied. look at row sum, or ave row sum, if it's above ~1.0,
	     * assume a counts matrix was supplied and the aver row sum,
	     * other wise, use the default (100).
	     */

	    /* is 1.2 a good choice??? */
	    if ( avePosSum > 1.2 ) {
	      /* assume counts were specified so use 'avePosSum' as the 
	       * number of counts/sites. */
	      modelNumSites = (long) floor(avePosSum + 0.5);
	    }
	    else
	      modelNumSites = (long) DEFAULT_FREQ_MODEL_NUM_SITES;
	  }
	  else if ( freqModelNumSites[numFreqModels] != 0.0 ) {
	    /* printf("HERE I AM, freqModelNumSites[i] != 0.0\n"); */

	    /* user specified a 'counts' for this model, we must use it. */
	    modelNumSites = freqModelNumSites[numFreqModels];

	    /* if a counts matrix was supplied, as opposed to a prob freq
	     * matrix, issue a warning if the value is considerably
	     * different from the ave counts. use the same test to determine
	     * if a counts matrix or a freq matrix was supplied. */
	    /* is 1.2 a good choice??? */
	    if ( avePosSum > 1.2 ) {
	      /* assume counts matrix was specified. if ave is really 
	       * different from user specified, issue warning. */
	      if ( fabs(avePosSum - (double)modelNumSites) > 2.0 ) { 
		fprintf(stderr,
	    "WARNING: FYI, for freq model %d, avePosSum (%f) is fairly different from user specified value (%ld); using user specified value.\n",
			numFreqModels+1,avePosSum, modelNumSites);
	      }
	    }
	    /* else, assume freq matrix was specified, so no warning. */
	  }
	  else {
	    /* printf("HERE I AM, freqModelNumSites[i] == 0.0\n"); */

	    /* user specified counts were supplied, but NOT for this matrix
	     * model. determine whether freq matrix or count matrix was
	     * specified and use either ave pos sum or default counts 
	     * value (100). */
	    /* is 1.2 a good choice??? */
	    if ( avePosSum > 1.2 ) {
	      /* assume counts were specified so use 'avePosSum' as the 
	       * number of counts/sites. */
	      modelNumSites = (long) floor(avePosSum + 0.5);
	    }
	    else
	      modelNumSites = (long) DEFAULT_FREQ_MODEL_NUM_SITES;
	  }
		     
	  printf("For freq matrix model %d, using position counts of %ld.\n",
		 numFreqModels+1,modelNumSites);

	  /* now, from the prob matrix just read in, place the counts (which
	   * were created from the probs) into the model's site_freq matrix. 
	   * remember, that here, the column order of modelProbMatrix[][]
	   * is different from what dscan wants (regardless of whether it's 
	   * DNA data or protein data). SetSModelFromProbCounts() takes
	   * care of the order switching.
	   * above, when a list of sites has beed specified for the model, 
	   * Add2SModel() is used to do this step.
	   */
	  SetSModelFromProbCounts(F->M[m],modelProbMatrix,modelNumSites);

	  if ( verbose) {
	    fprintf(stderr,"\n\n");
	    PutSModel(stderr,F->M[m]); 
	    fprintf(stderr,"\n\n"); 
	  }

	  /* incre the index to the user specified list of counts. */
	  ++numFreqModels;
	}
	else 
	  print_error("model file input error: 3");

   }  /* for ( c = ' ', length = 0, F->N=0, m=1; c!=EOF; m++ ) */

   if(verbose) fprintf(stderr,"\n\t %ld motif models\n",F->N);

   for(i=0; i<=MAX_BLOCK_SIZE; i++) free(seq[i]);
   free(seq); 

   for (j=0; j < MAX_BLOCK_LENGTH; j++) free(modelProbMatrix[j]);
   free(modelProbMatrix);

   free(null);

}  /* ReadScan */


long	MaskScan(FILE *fptr,long number, long *nsize, scn_typ F)   /* BT 9/26/97 */
{
	e_type	E;
	sm_type	*M=F->M;
	a_type	A = F->A;
	char	**mseq,*seq,c,*info;
	long	best,bestscore,**pos,s,score,i,j,item,*p,m,*w,tot_len;
	double	pval,k,factor;
	Boolean	okay;
	long	start,end,*p1=NULL,r,leng1;
	smx_typ	smx,*sM;

	if(F->repeats) print_error("Illegal sequence of operations");
	NEW(p,F->N+2,long); NEW(w,F->N+2,long);
	NEW(sM,F->N+2,smx_typ); NEWP(mseq,F->N+2,char);
	for(tot_len=0,m=1; m<=F->N; m++) { 
		w[m] = LenSModel(M[m]); tot_len += w[m];
	}
	smx = MkSMatrix(2.0,tot_len,F->freq,A);
	for(i=1,m=1; m<=F->N; m++) { 
	    for(j=1; j<=w[m]; j++,i++) { 
		for(c=0; c<=nAlpha(A); c++) { 
		   if(IsColumnSModel(j,M[m])){
			score = CellScoreSModel(c, j, M[m]); 
			SetSMatrix(c,i,score,smx);
		   }
		}
	    }
	    sM[m] = GetSMatrixSModel(M[m]);
	}
	/*** PutSMatrix(stderr,smx); /******/
	pos = F->pos;
        for(s=1; s<=number;s++) {
		E = ReadSeq(fptr,s,nsize[s],A);
		if(F->shuffle) ShuffleSeq(E);
		seq = SeqPtr(E);
		for(bestscore=0,factor= 1.0,okay=TRUE,m = 1; m<= F->N; m++){
		   best = INT_MIN;
		   mseq[m] = seq;
		   end = LenSeq(E) - w[m] +1;
		   if(end < 1) { okay=FALSE; break; }
		   factor *= (double) end;
		   for(i = 1; i<= end; i++){
		      score = ScoreSMatrix(seq, i, sM[m]);
		      if(score > best){ best = score; p[m]=i; }
		   }
		   bestscore += best;
		}
		if(!okay) NilSeq(E);
		else {
		   factor *= (double) number;
		   score = SplitScoreSMatrix(mseq, F->N, p, w, smx);
		   pval = SMatrixProb(score, smx);
		   pval = log10(pval*factor); 
		   /** fprintf(stderr,"score = %d = %d, pval = %f\n", 
			score, bestscore, pval); /****/
		   if(pval <= F->logmaxEval){
		     if(F->neg_mask){
/**** NEGATIVE MASKING *****/
			for(m = 1; m<= F->N; m++){
				if(m < F->N){ /** break if order not right **/
				  if((p[m+1]-p[m]) < w[m]) break;
				}
			}
		        if(s == 1){
			   if(m <= F->N) print_error("fatal: ordering error");
			   NEW(p1,F->N+2,long); 
			   leng1 = LenSeq(E);
			   for(m = 1; m<= F->N; m++) p1[m] = p[m];
			}
		        if(p1 == NULL) print_error("fatal: ordering error");
			if(m > F->N) {
			  printf(">");
			  PutSeqID(stdout,E);
			  printf(" (conserved regions)\n");
			  for(j=1,i=1,m = 1; m<= F->N; i=p1[m]+w[m],m++){
				for( ; i< p1[m]; i++,j++) {
					printf("x");
					if(j%50==0)printf("\n");
				}
				for(i=p[m]; i< p[m]+w[m]; i++,j++) {
					r = ResSeq(i,E);
					printf("%c",AlphaChar(r,A));
					if(j%50==0)printf("\n");
				}
			  }
			  for(i=p1[F->N]+w[F->N]; i<= leng1; i++,j++) {
				printf("x");
				if(j%50==0)printf("\n");
			  }
			  printf("\n");
			}
/**** NEGATIVE MASKING *****/
		     } else {
			for(m = 1; m<= F->N; m++){
				MaskSeq(p[m],(p[m]+w[m]-1), E);
			}
		     }
		   } 
		   if(!F->neg_mask) PutSeq(stdout,E,A);
		   NilSeq(E);
		}
	}
	if(F->neg_mask && p1 != NULL) { free(p1); }
	free(p); free(w); free(sM); free(mseq); NilSMatrix(smx);
}

long	ScanScan(FILE *fptr,long number, long *nsize, scn_typ F)
{
	e_type	E;
	sm_type	*M=F->M;
	mh_type H = F->H;
	a_type	A = F->A;
	h_type	HG = F->HG;
	char	**mseq,*seq,c,*info;
	long	best,bestscore,**pos,s,score,i,j,end,item,*p,m,*w,tot_len;
	double	pval,k,factor,ave;
	Boolean	okay;
	smx_typ	smx,*sM;
	double  p1;
	int     totalLen;

	if(F->repeats) print_error("Illegal sequence of operations");
	NEW(p,F->N+2,long);
	NEW(w,F->N+2,long);
	NEW(sM,F->N+2,smx_typ);
	NEWP(mseq,F->N+2,char);
	for(tot_len=0,m=1; m<=F->N; m++) { 
		w[m] = LenSModel(M[m]); tot_len += w[m];
	}
	smx = MkSMatrix(2.0,tot_len,F->freq,A);
	for(i=1,m=1; m<=F->N; m++) { 
	    for(j=1; j<=w[m]; j++,i++) { 
		for(c=0; c<=nAlpha(A); c++) { 
		   if(IsColumnSModel(j,M[m])){
			score = CellScoreSModel(c, j, M[m]); 
			SetSMatrix(c,i,score,smx);
		   }
		}
	    }
	    sM[m] = GetSMatrixSModel(M[m]);
	}

	totalLen = CountNumSites( fptr, number,  nsize, F, w, F->numSites );	

	/*** PutSMatrix(stderr,smx); /******/
	pos = F->pos;
	fprintf(stderr,"searching database of %ld sequences\n",number);
        for(ave=0,s=1; s<=number;s++) {
	  E = ReadSeq(fptr,s,nsize[s],A);
	  /******
		 printf("#%");
		 PutSeqID(stdout,E);
		 printf("{%d}\n",LenSeq(E));
		 /******/
	  ave += LenSeq(E);
	  /******/
	  if(F->shuffle) ShuffleSeq(E);
	  seq = SeqPtr(E);
	  for(bestscore=0,factor=1.0,okay=FALSE,m = 1; m<= F->N; m++){
	    best = INT_MIN;
	    mseq[m] = seq;
	    end = LenSeq(E) - w[m] +1;
	    if(end < 1) { okay=FALSE; break; }
	    if( ! F->adjustN )
	      factor *= (double) end;  /* BT 11/05/2001 */
	    else
	      factor *= (double) F->numSites[m];   /* BT 11/05/2001 */
	    for(i = 1; i<= end; i++){
	      score = ScoreSMatrix(seq, i, sM[m]);
	      if(score > best)
		{ 
		  okay=TRUE;    /* BT 01/30/04 */
		  if( ! F->nNoOverlap )    /* BT 10/4/2000 */
		    {
		      best = score; 
		      p[m]=i;
		    }
		  else
		    {
		      if( !Overlap( m, i, p, sM ) )
			{				
			  best = score; 
			  p[m]=i;
			}
		    }
		}
	    }
	    bestscore += best;
	  }
	  if(!okay) NilSeq(E);
	  else {
	    if( ! F->adjustN )
	      factor *= (double) number;   /* BT 11/05/2001 */
	    score = SplitScoreSMatrix(mseq, F->N, p, w, smx);
	    pval = SMatrixProb(score, smx);
	    p1 = pval;
	    pval = log10(pval*factor); 
	    /* fprintf(stderr,"seq = %d score = %ld best = %d, pval = %f factor = %f p1 = %g\n", 
	       s, score, bestscore, pval, factor, p1); /****/            /* DEBUG */
	    IncdHist(-pval,HG);
	    if(pval <= F->logmaxEval && 
	       (item = InsertMheap((keytyp)pval,H)) != 0){
	      for(m = 1; m<= F->N; m++){
		pos[m][item] = p[m];
	      }
	      if(F->E[item] != NULL) NilSeq(F->E[item]);
	      F->E[item] = E;
	    } else NilSeq(E);
	  }
	}
	printf("average length = %.1lf\n",ave/(double)number);
	free(p); free(w); free(sM); free(mseq); NilSMatrix(smx);
}

long	PutScan(FILE *fptr, scn_typ F)
{
  double	score,s;
  long	i,j,m,end;
  FILE	*ofptr;
  Boolean	flag;
  
  if(F->repeats) print_error("Illegal sequence of operations");
  PutHist(stdout,55,F->HG);
  if(F->seqfile) ofptr = open_file(F->snfile,".seq","w");
  while(TRUE)
    {
      score = -(double) MinKeyMheap(F->H);
      if((i=DelMinMheap(F->H))==0) break;
      /*** NEW: determine if score okay. ****/
      for(flag = TRUE, m = 1; m<= F->N; m++)
	{
	  end = LenSeq(F->E[i]) - LenSModel(F->M[m]) + 1;
	  s = PvalSModel(SeqPtr(F->E[i]),F->pos[m][i],F->M[m]);
	  if( F->nCountFlag )    /*BT 10/21/97 */
	    {
	      if( F->nCount <= 0 )
		{
		  flag = FALSE;
		  break;
		}		    
	    }
	  else
	    {
	      if( (-log10(s*end)) <= F->singlePval) 
		{ 
		  flag=FALSE; 
		  break;
		}
	    }
	}
      /*** NEW ****/
      if(flag)
	{
	  F->nCount--;
	  fprintf(fptr,"[%3.2f] ", score);
	  PutSeqInfo(fptr,F->E[i]);
	  for(m = 1; m<= F->N; m++)
	    {
	      end = LenSeq(F->E[i]) - LenSModel(F->M[m]) + 1;
	      score = PvalSModel(SeqPtr(F->E[i]),
				 F->pos[m][i], F->M[m]);
	      fprintf(fptr," %4.1f (%8.5e): ", -log10(score*end), score); /* BT 06/04/04 */
	      /*** OFF ****
		   fprintf(fptr,"(%4d) ",
		   ScoreSModel(SeqPtr(F->E[i]),
		   F->pos[m][i], F->M[m]));
		   /*** ****/
	      PutSeqRegion(fptr,F->pos[m][i],
			   LenSModel(F->M[m]),F->E[i],F->A);
	      fprintf(fptr,"\n");
	    }
	  fprintf(fptr,"\n\n");
	  if(F->seqfile) PutSeq(ofptr,F->E[i],F->A);
	}
    }
  if(F->seqfile)  fclose(ofptr);
}


long	OrderScanScan(FILE *fptr,long number, long *nsize, scn_typ F)
{
	e_type	E;
	sm_type	*M=F->M;
	mh_type H = F->H;
	a_type	A = F->A;
	h_type	HG = F->HG;
	char	**mseq,*seq,c,*info;
	long	best,bestscore,**pos,s,score,i,j,end,item,*p,m,*w,tot_len,n;
	long	*tmp;
	double	pval,k,factor;
	Boolean	okay;
	smx_typ	smx,*sM;
	mh_type	*mhp;
	long	**site,**scores,numsave=20,**Items;

	if(F->repeats) print_error("Illegal sequence of operations");
	NEW(p,F->N+2,long);
	NEW(w,F->N+2,long);
	NEW(sM,F->N+2,smx_typ);
	NEWP(mseq,F->N+2,char);
	for(tot_len=0,m=1; m<=F->N; m++) { 
		w[m] = LenSModel(M[m]); tot_len += w[m];
	}
	NEW(mhp,F->N+2,mh_type);
	NEWP(site,F->N+2,long);
	NEWP(scores,F->N+2,long);
	NEWP(Items,F->N+2,long);
	for(m=1; m<=F->N; m++) { 
		mhp[m]=Mheap(numsave,3);
		NEW(site[m],numsave+2,long);
		NEW(scores[m],numsave+2,long);
		NEW(Items[m],numsave+2,long);
	}
	smx = MkSMatrix(2.0,tot_len,F->freq,A);
	for(i=1,m=1; m<=F->N; m++) { 
	    for(j=1; j<=w[m]; j++,i++) { 
		for(c=0; c<=nAlpha(A); c++) { 
		   if(IsColumnSModel(j,M[m])){
			score = CellScoreSModel(c, j, M[m]); 
			SetSMatrix(c,i,score,smx);
		   }
		}
	    }
	    sM[m] = GetSMatrixSModel(M[m]);
	}
	/*** PutSMatrix(stderr,smx); /******/
	pos = F->pos;
	fprintf(stderr,"searching database\n");
        for(s=1; s<=number;s++) {
	  E = ReadSeq(fptr,s,nsize[s],A);
	  if(F->shuffle) ShuffleSeq(E);
	  seq = SeqPtr(E);
	  for(bestscore=0,okay=TRUE,m = 1; m<= F->N; m++){
		   best = INT_MIN;
		   mseq[m] = seq;
		   end = LenSeq(E) - w[m] +1;
		   if(end < 1) { okay=FALSE; break; }
		   for(i = 1; i<= end; i++){
		     if( seq[i] == '\0' )       /* BT 12/19/97 */
		       score = -LONG_MAX;
		     else
		       score = ScoreSMatrix(seq, i, sM[m]);
		     if(score > best){ best = score; p[m]=i; }
		   }
		   bestscore += best;
	  }
/*** fprintf(stderr,"DEBUG: down - %d\n",s); /****/
	  for(i=0, m = 1; m<= F->N; m++){ i+= w[m]; }
	  if(i > (long) LenSeq(E) || !okay) NilSeq(E);
	  else {
	   n = LenSeq(E) - tot_len; m = F->N;
	   factor = bico(n+m,m) * (double) number;
	   score = SplitScoreSMatrix(mseq, F->N, p, w, smx);
	   pval = SMatrixProb(score, smx);
	   pval = log10(pval*factor); 
	   /** fprintf(stderr,"score = %d = %d, pval = %f\n", 
			score, bestscore, pval); /****/
	   /** if best is within limit then do rigorous check **/
	   if(pval > F->logmaxEval) NilSeq(E);
	   else {
	     /*** get best regions ***/
	     for(m = 1; m<= F->N; m++){		
		end = LenSeq(E) - w[m] +1;
		for(i = 1; i<= end; i++){
		    score = ScoreSMatrix(seq, i, sM[m]);
		    if((item=InsertMheap(-(keytyp)score, mhp[m])) !=0){
			scores[m][item]=score;
			site[m][item]=i;
		    }
		}
	     }
	     /*** find best ordered motifs ***/
	     for(m = 1; m<= F->N; m++){		
	        for(i=0;(item=DelMinMheap(mhp[m]))!=0; i++){
			Items[m][i] = item;
		}
		Items[m][i] = 0;
	     }
	     NEW(tmp, F->N+3,long);
	     score = INT_MIN;
	     foo_scan(1,1,F->N,Items,scores,site,p,w,tmp,0,&score);
	     free(tmp);
	     /*** calc pval for ordered motifs ***/
	     pval = SMatrixProb(score, smx);
	     pval = log10(pval*factor); 
	     if(pval <= F->logmaxEval && 
		(item = InsertMheap((keytyp)pval,H)) != 0){
	   	   IncdHist(-pval,HG); 
		   for(m = 1; m<= F->N; m++){
				pos[m][item] = p[m];
		   }
		   if(F->E[item] != NULL) NilSeq(F->E[item]);
		   F->E[item] = E;
	     } else NilSeq(E);
	   }
	  }
/*** fprintf(stderr,"DEBUG: up\n"); if(s>=20) exit(1); /****/
	}
	free(p); free(w); free(sM); free(mseq); NilSMatrix(smx);
	for(m=1; m<=F->N; m++) { 
		NilMheap(mhp[m]);
		free(site[m]); free(scores[m]); free(Items[m]);
	}
	free(site); free(scores); free(Items); free(mhp);
}

long	foo_scan(long start, long m, long N, long **Items, long **scores, 
		long **site, long *path, long *w, long *tmp, long score, long *best)
{
	long	i,item,s;

	if(m > N) {
		if(score > *best){
			*best = score;
			for(i=1;i<=N;i++) path[i]=tmp[i];
		} return 0;
	}
	for(i=0; (item = Items[m][i]) !=0; i++){
	    if((s = site[m][item]) >= start){
		tmp[m] = s;
		foo_scan(s+w[m],m+1,N,Items,scores,site,path,w,tmp,
			(score + scores[m][item]), best);
	    }
	}
	return 0;
}

long	ScanRScan(FILE *fptr,long number, long *nsize, scn_typ F)
{
	e_type	E;
	sm_type	M=F->M[1];
	mh_type sH,H = F->H;
	a_type	A = F->A;
	h_type	HG = F->HG;
	char	*seq,c,*info;
	long	**pos,s,i,j,k,end,item,p,m,w,*tpos;
	long	bestnsegs,nsegs,*scores,score,*sites,item2,obs,bestobs;
	double	pval,bestpval,factor,lambda,weight;
	double	sum_p,Vi,d_fact,N,R,qi,qim1;
	long	okay;
	smx_typ	sM;
	st_type	S=NULL;
	long	len_mtf[3];
	ss_type P;

	if(F->N > 1) print_error("only one block allowed for repeats");
	F->repeats=TRUE;
	w = LenSModel(M);
	sM = GetSMatrixSModel(M);
	/*	PutSMatrix(stderr,sM); /******/
	fprintf(stderr,"min_segs = %ld; max_segs = %ld\n\n",
		F->min_segs, F->max_segs);
	for(i=1;i<=F->N; i++) free(F->pos[i]); free(F->pos); 
	NEWP(F->pos,F->hpsz+3,long);
	for(i=1;i<=F->hpsz; i++) { NEW(F->pos[i],F->max_segs+3,long); }
	pos = F->pos;
	NEW(scores,F->max_segs+2,long);
	NEW(sites,F->max_segs+2,long);
	NEW(tpos,F->max_segs+2,long);
	fprintf(stderr,"searching database (%ld sequences)\n",number);
        for(s=1; s<=number;s++) {
		E = ReadSeq(fptr,s,nsize[s],A); seq = SeqPtr(E);
		if(F->shuffle) ShuffleSeq(E);
		okay=0; sH = NULL;
		end = LenSeq(E) - w +1;
		if(end > 0) {
		  for(i=1; i<= end; i++){

		    if( seq[i] == '\0' )       /* BT 12/19/97 */
		      score = -LONG_MAX;
		    else
		      score = ScoreSMatrix(seq, i, sM);

		      if(score > 0) {
			if(sH==NULL) {
			   sH=Mheap(F->max_segs,3); /** heap for best scores **/
			}
			if((item2 = InsertMheap(-(keytyp)score,sH)) != 0){
			   	okay++;
				scores[item2] = score;
				sites[item2] = i;
			}
		      }
		  }
		  factor = (double) end;
		}
		if(okay >= F->min_segs){
		   weight = 2.0; bestpval = 99999; 
		   /*** sites datatype stuff ***/
		   len_mtf[1]= w; 
		   P = SeqSet1("none", E,A);
                   S = MkSites(1,len_mtf,P);
		   /*** sites datatype stuff ***/
		   bestobs = 0;  obs=1;
		   sum_p = 0.0; Vi = 1.0; d_fact = 1.0; N = factor;
		   qim1 = 1.0;
		   while((item2=DelMinMheap(sH))!=0) {
		     k = sites[item2];
                     if(!OccupiedSite(1,1,k,S)){
		   	score = scores[item2];
		        lambda = SMatrixProb(score, sM);
			/*************** New Order statistics *************/
			R = (double) obs;
			qi = (1.0 - lambda);
			Vi = pow(qi/qim1,N-R+1);
			/** if get same score twice... *****/
			if(Vi == 1.0) {
				printf("qi/qim1 = %g/%g = %g\n",
					qi,qim1,qi/qim1);
				printf("score = %ld; Vi = %g; N = %g; R = %g\n",
					score,Vi, N, R); 
			}
			printf("1 - lambda = %g; qi = %g; qim1 = %g\n",
				1.0 - lambda,qi,qim1);
	/******/
			qim1 = qi;
			sum_p += (1.0 - Vi);
			d_fact *= obs;
			pval = weight * pow(sum_p,R)/d_fact;
			/******/
			if(sum_p >= 1.0) printf("sum_p = %g\n",sum_p);
			printf("score = %ld; Vi = %g (i = %ld) pval = %g\n\n",
				score,Vi,obs,pval);
	/******/
			if(obs >= F->min_segs) {
			    if(obs == F->max_segs) pval /= 2.0;
			    if(pval < bestpval)
			      {
				AddSite(1,1,k,S);
				bestpval = pval; bestobs = obs;
				tpos[obs] = sites[item2];
				obs++;
				weight *= 2.0;
			      } 
			    /******/    /* BT 10/21/97 */
			    else 
			      {
				AddSite(1,1,k,S);
				tpos[obs] = sites[item2];
				obs++;
				weight *= 2.0;
			      }  /*** look for better one at higher obs **/
			    /*** else weight *= 2.0; /*** keep going ***/
			    /** else break; /*** don't use ***/   /* BT 10/21/97 */
			} else {
                                AddSite(1,1,k,S);
				tpos[obs] = sites[item2];
				obs++;
				weight = 2.0;
			} 
		      }
		   }
		   NilSites(S); NilSeqSet(P);
/****
pval = bestpval;
/****/
		   bestpval *= (double) number; /* adjust for # of sequences */
		   bestpval = log10(bestpval); 
		   IncdHist(-bestpval,HG);
		   if(bestpval <= F->logmaxEval && 
			(item = InsertMheap((keytyp)bestpval,H)) != 0){
/****
PutSeqID(stdout,E);
printf(": %g (%d * %g)\n",bestpval,number,pval);
/****/
			if(F->E[item] != NULL) NilSeq(F->E[item]);
		   	F->E[item] = E;
			/*			for(i=1; i<= bestobs; i++){  */
			for(i=1; i < obs; i++){    /* BT 10/21/97 */
/****
pval = PvalSModel(SeqPtr(E),tpos[i], F->M[1]);
fprintf(stdout," %4.1f: ", -log10(pval*factor));
PutSeqRegion(stdout,tpos[i],LenSModel(F->M[1]),E,F->A);
printf("\n");
/****/

				pos[item][i] = tpos[i]; 
			}
			pos[item][i] = 0;
		   } else NilSeq(E);
		   if(sH != NULL) { NilMheap(sH);  sH = NULL; }
		} else {
		   if(sH != NULL) { NilMheap(sH);  sH = NULL; }
		   NilSeq(E);
		}
	}
	free(tpos); free(scores); free(sites);
}

long     PutRScan(FILE *fptr, scn_typ F)
{
  double  score;
  long     i,j,m,end;
  FILE    *ofptr;
  short   flag = TRUE;   /* BT 10/21/97 */
  
  if(!F->repeats) print_error("Illegal sequence of operations");
  PutHist(stdout,55,F->HG);
  if(F->seqfile) ofptr = open_file(F->snfile,".seq","w");
  while(flag)
    {
      score = -(double) MinKeyMheap(F->H);
      if((i=DelMinMheap(F->H))==0) break;
      else 
	{
	  if( F->nCountFlag )    /*BT 10/21/97 */
	    {
	      if( F->nCount <= 0 )
		{
		  flag = FALSE;
		  break;
		}		    
	    }
	  F->nCount--;
	  fprintf(fptr,"[%3.2f] ", score);
	  PutSeqInfo(fptr,F->E[i]);
	  end = LenSeq(F->E[i]) - LenSModel(F->M[1]) + 1;
	  for(j=1; F->pos[i][j] != 0; j++)
	    {
	      score = PvalSModel(SeqPtr(F->E[i]),F->pos[i][j], F->M[1]);
	      if( (! F->nCountFlag) && (-log10(score*(double)end)) <= F->singlePval) /* BT 04/29/2002 */
		{ 
 		  break;
		}
	      fprintf(fptr," %4.1f (%8.5e): ", -log10(score*(double)end), score); /* BT 06/04/04 */
	      PutSeqRegion(fptr,F->pos[i][j],
			   LenSModel(F->M[1]),F->E[i],F->A);
	      fprintf(fptr,"\n");
	    }
	  fprintf(fptr,"\n\n");
	  if(F->seqfile) PutSeq(ofptr,F->E[i],F->A);
	}
    }
  if(F->seqfile)  fclose(ofptr);
}

long	OrderScanMtf(char *name, scn_typ F)
{
	e_type	E;
	sm_type	*M=F->M;
	mh_type H = F->H;
	a_type	A = F->A;
	h_type	HG = F->HG;
	char	**mseq,*seq,c,*info;
	long	best,bestscore,**pos,s,score,i,j,end,item,*p,m,*w,tot_len,n;
	long	*tmp;
	double	pval,k,factor;
	Boolean	okay;
	smx_typ	smx,*sM;
	mh_type	*mhp;
	long	**site,**scores,numsave=20,**Items;

	/*** NEW ****/
	long	number,*len_elem;
	double	**prob,d;
	ss_type	P;
	st_type S;
	Boolean	*null;
	/*** NEW ****/

	if(F->repeats) print_error("Illegal sequence of operations");

	NEW(p,F->N+2,long);
	NEW(w,F->N+2,long);
	NEW(sM,F->N+2,smx_typ);
	NEWP(mseq,F->N+2,char);
	for(tot_len=0,m=1; m<=F->N; m++) { 
		w[m] = LenSModel(M[m]); tot_len += w[m];
	}
/*** NEW ****/
	P = SeqSet(name,A);
	number = NSeqsSeqSet(P);
	NEW(len_elem,F->N+2,long);
	NEWP(prob,number+2,double);
	for(i=1; i<=number; i++) {
		NEW(prob[i], SqLenSeqSet(i,P) +2,double);
	}
	for(m=1; m<=F->N; m++) len_elem[m] = LenSModel(M[m]);
	S = MkSites(F->N, len_elem, P);
/*** NEW ****/
	NEW(mhp,F->N+2,mh_type);
	NEWP(site,F->N+2,long);
	NEWP(scores,F->N+2,long);
	NEWP(Items,F->N+2,long);
	for(m=1; m<=F->N; m++) { 
		mhp[m]=Mheap(numsave,3);
		NEW(site[m],numsave+2,long);
		NEW(scores[m],numsave+2,long);
		NEW(Items[m],numsave+2,long);
	}
	smx = MkSMatrix(2.0,tot_len,F->freq,A);
	for(i=1,m=1; m<=F->N; m++) { 
	    for(j=1; j<=w[m]; j++,i++) { 
		for(c=0; c<=nAlpha(A); c++) { 
		   if(IsColumnSModel(j,M[m])){
			score = CellScoreSModel(c, j, M[m]); 
			SetSMatrix(c,i,score,smx);
		   }
		}
	    }
	    sM[m] = GetSMatrixSModel(M[m]);
	}
	/*** PutSMatrix(stderr,smx); /******/
	pos = F->pos;
	fprintf(stderr,"searching database\n");
        for(s=1; s<=number;s++) {
	  E = SeqSetE(s,P);
	  if(F->shuffle) ShuffleSeq(E);
	  seq = SeqPtr(E);
	  for(bestscore=0,okay=TRUE,m = 1; m<= F->N; m++){
		   best = INT_MIN;
		   mseq[m] = seq;
		   end = LenSeq(E) - w[m] +1;
		   if(end < 1) { okay=FALSE; break; }
		   for(i = 1; i<= end; i++){
		      score = ScoreSMatrix(seq, i, sM[m]);
		      if(score > best){ best = score; p[m]=i; }
		   }
		   bestscore += best;
	  }
/*** fprintf(stderr,"DEBUG: down - %d\n",s); /****/
	  for(i=0, m = 1; m<= F->N; m++){ i+= w[m]; }
	  if(i > (long) LenSeq(E) || !okay) NilSeq(E);
	  else {
	   n = LenSeq(E) - tot_len; m = F->N;
	   factor = bico(n+m,m) * (double) number;
	   score = SplitScoreSMatrix(mseq, F->N, p, w, smx);
	   pval = SMatrixProb(score, smx);
	   pval = log10(pval*factor); 
	   /** fprintf(stderr,"score = %d = %d, pval = %f\n", 
			score, bestscore, pval); /****/
	   /** if best is within limit then do rigorous check **/
	   if(pval > F->logmaxEval) NilSeq(E);
	   else {
	     /*** get best regions ***/
	     for(m = 1; m<= F->N; m++){		
		end = LenSeq(E) - w[m] +1;
		for(i = 1; i<= end; i++){
		    score = ScoreSMatrix(seq, i, sM[m]);
		    if((item=InsertMheap(-(keytyp)score, mhp[m])) !=0){
			scores[m][item]=score;
			site[m][item]=i;
		    }
		}
	     }
	     /*** find best ordered motifs ***/
	     for(m = 1; m<= F->N; m++){		
	        for(i=0;(item=DelMinMheap(mhp[m]))!=0; i++){
			Items[m][i] = item;
		}
		Items[m][i] = 0;
	     }
	     NEW(tmp, F->N+3,long);
	     score = INT_MIN;
	     foo_scan(1,1,F->N,Items,scores,site,p,w,tmp,0,&score);
	     free(tmp);
	     /*** calc pval for ordered motifs ***/
	     pval = SMatrixProb(score, smx);
	     pval = log10(pval*factor); 
	     if(pval <= F->logmaxEval){
/*** NEW ***/
		for(m = 1; m<= F->N; m++){
			end = LenSeq(E) - w[m] +1;
			AddSite(m, s, p[m], S);
			d = PvalSModel(SeqPtr(E),p[m], M[m]);
                    	prob[s][p[m]] = -log10(d*(double)end);
		}
/*** NEW ***/
		
		if((item = InsertMheap((keytyp)pval,H)) != 0){
	   	   IncdHist(-pval,HG); 
		   for(m = 1; m<= F->N; m++){
				pos[m][item] = p[m];
		   }
		   if(F->E[item] != NULL) NilSeq(F->E[item]);
		   F->E[item] = CopySeq(E);
		}
	     }
	   }
	  }
/*** fprintf(stderr,"DEBUG: up\n"); if(s>=20) exit(1); /****/
	}
/*** NEW ****/
	for(m = 1; m<= F->N; m++){
		NEW(null,LenSModel(M[m])+3,Boolean);
		for(i=1; i <= LenSModel(M[m]); i++){
			null[i] = NullSiteSModel(i,M[m]);
		}
		PutSitesMtfDBS(stdout,m,S,prob,null);
	}
	NilSites(S);
	NilSeqSet(P);
/*** NEW ****/
	free(p); free(w); free(sM); free(mseq); NilSMatrix(smx);
	for(m=1; m<=F->N; m++) { 
		NilMheap(mhp[m]);
		free(site[m]); free(scores[m]); free(Items[m]);
	}
	for(i=1; i<=number; i++) free(prob[i]);
	free(prob);
	free(site); free(scores); free(Items); free(mhp);
}


int Overlap( int m, int i, long *p, smx_typ *sM)
{
  int     overlap = FALSE;
  int     m1;
  int     i1;
 
  for( m1 = 1; m1 < m && !overlap; m1++ )
    {
      for( i1 = i; i1 < i + sM[m]->K; i1++ )
	{
	  if( i1 == p[m1] )
	    {
	      overlap = TRUE;
	      break;
	    }
	}
    
      for( i1 = p[m1]; i1 < p[m1] + sM[m1]->K && !overlap; i1++ )
	{
	  if( i1 == i )
	    {
	      overlap = TRUE;
	      break;
	    }
	}
    }

  return(overlap);
}


int CountNumSites( FILE *fptr, long number,  long *nsize, scn_typ F, long *w, int *numSites )
{
  int    s;
  int    m;
  int    n;
  int    j;
  int    totalLen = 0;
  e_type E;
  a_type A = F->A;
  int    siteOK;

  for( m = 1; m<= F->N; m++)
    {
      numSites[m] = 0;
    }

  for( s = 1; s <= number; s++ )
    {
      totalLen += nsize[s];
      E = ReadSeq(fptr, s, nsize[s], A);
      for( m = 1; m<= F->N; m++)
	{
	  for( n = 1; n <= nsize[s] - w[m] + 1; n++ )
	    {
	      if( E->S[n] )
		{
		  siteOK = TRUE;
		  for( j = 0; j < w[m]; j++ )
		    {
		      if( ! E->S[n+j] ) 
			{
			  siteOK = FALSE;
			  break;
			}			
		    }
		}
	      else
		siteOK = FALSE;
	      if( siteOK )
		numSites[m]++;
	    }
	}
    }
  fseek( fptr, 0, SEEK_SET );

  fprintf( stdout, "Total database length: %d\n", totalLen );
  for( m = 1; m<= F->N; m++)
    {
      fprintf( stdout, "Effective size for model %d: %d\n", m, numSites[m] );
    }
  fprintf( stdout, "\n" );

  return( totalLen );
}
