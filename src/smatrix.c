/* smatrix.c 
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/smatrix.c,v 2.2 2004/06/07 19:03:42 thompson Exp $
 *
 */

#include "smatrix.h"

smx_typ	MkSMatrix(double nsd, long K, double *freq, a_type A)
/* create a scoring matrix with all zero scores. */
{
	smx_typ	M;
	long	i,c;
	
	NEW(M,1,smatrix_type);
	M->A = A;
	M->nlet = nAlpha(A);
	M->K = K;
	M->nsd = nsd;
	M->neginf = (LONG_MIN/K) + 1;	/*** NEW ***/
	M->calc_prob = TRUE;
	M->calc_stats = TRUE;
	M->changed = TRUE;
	M->f0 = NULL;
	NEWP(M->score,K+2,long);
	NEW(M->maxseg,K+3,char);
	for(i=1; i<=K;i++) {
		NEW(M->score[i],M->nlet+2,long);
		for(c=0; c<=M->nlet; c++) M->score[i][c]=0;
	}
	NEW(M->max,K+2,long); NEW(M->min,K+2,long); 
	NEW(M->cmax,K+2,long); NEW(M->cmin,K+2,long);
	NEW(M->freq,M->nlet +2,double);
	for(c=0; c<=M->nlet; c++) M->freq[c]=freq[c];
	return M;
}

void	NilSMatrix(smx_typ M)
{
	long	i;

	for(i=1; i <= M->K;i++) { free(M->score[i]); }
	free(M->score); 
	if(M->f0 != NULL) free(M->f0);
	free(M->cmin); free(M->cmax);
	free(M->min); free(M->max);
	free(M->freq); free(M->maxseg);
	free(M);
}

long	SetSMatrix(long r, long row , long score, smx_typ M)
/* set the matrix for residue r in row equal to score */
{
	if(row < 1 || row > M->K) smatrix_error("input error");
	M->score[row][r]= MAX(long,score,M->neginf); /*** NEW ***/
	M->changed = TRUE;
}

void	PutSMatrix(FILE *fptr, smx_typ M)
{
	long	j,d;
	a_type	A = M->A;

	fprintf(fptr,"\n    ");
	for(d=0; d<=nAlpha(A); d++){
		fprintf(fptr,"  %c ", AlphaChar(d,A));
	}
	fprintf(fptr,"\n");
	for(j=1; j <= M->K; j++){
		fprintf(fptr,"%2ld: ",j);
		for(d=0; d<=nAlpha(A); d++){
			fprintf(fptr,"%3ld ", M->score[j][d]);
		}
		fprintf(fptr,"\n");
	}
	fprintf(fptr,"\n");
}

long	SplitScoreSMatrix(char **seq, long n, long *start, long *leng, smx_typ M)
/* determine score by split regions. WARNING: assumes lengths sum to M->K */
{
	register long	m,j,c,s,score;

	for(score=0,m=1,c=1; m <= n; m++){
	    for(s=start[m],j=0; j < leng[m]; s++,j++,c++){
		score += M->score[c][seq[m][s]];
	    }
	}
	return score;
}

long	ScoreSMatrix(register char *seq, register long start, register smx_typ M)
{
  register long	j,score;

  for(seq += start-1, score=0,j = 1; j <= M->K; j++)
    {
      if( seq[j] == '\0' )
	{
	  score = LONG_MIN;  /* BT 01/30/04 */ 
	  return score;
	}
      score += M->score[j][seq[j]];
    }
  return score;
}

double	SMatrixProbFast(long score, smx_typ M)
{
	double	f0[3],n;

	if(M->changed){M->changed=FALSE;M->calc_prob=M->calc_stats=TRUE;} 
	if(M->calc_stats) stats_smatrix(M);
	if(score <= M->cmin[M->K]) return 1.0;
	if(M->calc_prob){
		if(M->f0!= NULL){ free(M->f0); M->f0=NULL; }
		f0[0] = 1.0; min_max_smatrix(M);
		M->fx = NULL;
		smatrix_prob(0, f0, M);
		M->calc_prob = FALSE; 
	}
	if(score < M->cmin[M->K]) return 1.0;
	else return M->f[score];
}

double	SMatrixProb(long score, smx_typ M)
{
	double	f0[3],n;

	if(M->changed){
		stats_smatrix(M);
		if(M->f0!= NULL){ free(M->f0); M->f0=NULL; }
		f0[0] = 1.0; min_max_smatrix(M);
		M->fx = NULL;
		smatrix_prob(0, f0, M);
		M->changed = M->calc_prob = FALSE; 
	}
	if(score < M->cmin[M->K]) return 1.0;
	else return M->f[score];
}

char	*MaxSegSMatrix(smx_typ M)
{
	long	s,i,r,K=M->K;
	double	f0[3];

	if(M->changed){
		stats_smatrix(M);
		if(M->f0!= NULL){ free(M->f0); M->f0=NULL; }
		f0[0] = 1.0; min_max_smatrix(M);
		M->fx = NULL;
		smatrix_prob(0, f0, M);
		M->changed = M->calc_prob = FALSE; 
	}
	for(i=1; i<=K; i++){
	   for(r = 0; r <= M->nlet; r++){
		s = M->score[i][r];
	    	if(s == M->max[i]){
			M->maxseg[i]=r; break;
		}
	   }
	}
	M->maxseg[i]='\0';
	return M->maxseg;
}

long	MaxScoreSMatrix(smx_typ M)
{
	double	f0[3];

	if(M->changed){
		stats_smatrix(M);
		if(M->f0!= NULL){ free(M->f0); M->f0=NULL; }
		f0[0] = 1.0; min_max_smatrix(M);
		M->fx = NULL;
		smatrix_prob(0, f0, M);
		M->changed = M->calc_prob = FALSE; 
	}
	return M->cmax[M->K];
}

/****************************** private **********************************/
long	stats_smatrix(smx_typ M)
{
	long	i,r;
	double	E,mean,variance,V;
	
	for(mean=variance=0.0,i=1; i <= M->K; i++){
	    for(E=0.0,r=0; r <= M->nlet; r++){
		E += (double) M->score[i][r] * M->freq[r];
	    }
	    mean += E;
	    for(V=0.0,r=0; r <= M->nlet; r++){
		V += pow((E - (double) M->score[i][r]),2.0);
	    }
	    variance += V/(M->nlet - 1);
	}
	M->mean = mean; M->var = variance; M->sd = sqrt(variance);
	M->calc_stats = FALSE;
}

long	min_max_smatrix(smx_typ M)
/********************************************************************
      cmin0[K]                       n sd              cmax[K] 
	|-----------------------------|-----------------|

	cmin0[K-1]                cmin[K-1]        cmax[K-1]
	   |------------------------|-----------------|
		:		:		:
	        cmin0[1] cmin[1]     cmax[1]
	             |----|------------|

	                cmin[0]=cmax[0]=0
			       |

	cmax[0]   = 0
	cmax[i]   = cmax[i-1] + max[i];

	cmin[K]   = mean + n*sd;  where n = nsd
	cmin[i-1] = cmin[i] - max[i];
**********************************************************************/
{
	long	s,i,r,min,K=M->K;
	double	minscore;

	if(M->calc_stats) stats_smatrix(M);
	M->cmax[0] = 0;
	for(M->max[0]=M->min[0]=0,i=1; i<=K; i++){
	   M->max[i]=LONG_MIN; M->min[i]=LONG_MAX;
	   for(r = 0; r <= M->nlet; r++){
		s = M->score[i][r];
	   	M->max[i]=MAX(long,M->max[i],s);
	   	M->min[i]=MIN(long,M->min[i],s);
	   }
	   M->cmax[i] = M->cmax[i-1] + M->max[i];
	}
	minscore = (M->mean + M->nsd*M->sd);
	if( minscore >  M->cmax[K] )
          M->cmin[K] = 0;
        else if(minscore  >  -99.0){
                M->cmin[K] = (long)minscore - 1;
        } else M->cmin[K] = 0;  /** beware of negative infinity socres **/
	for(i=K-1; i>=1; i--){ M->cmin[i] = M->cmin[i+1] - M->max[i+1]; }
	M->cmin[0] = 0;
	for(min = 0,i=1; i<=K; i++){
	    min += M->min[i]; M->cmin[i] = MAX(long,M->cmin[i],min);
	}
	/** fprintf(stderr,"cmax=%d; cmin=%d; mean = %f\n",
		M->cmax[K],M->cmin[K],M->mean); /****/
	/** M->cmin[i] = min; /*** DEBUG ***/
	return M->cmax[K];
}

Boolean	smatrix_prob(register long k, register double *f, smx_typ M)
/*************************************************************************

  Call: smatrix_prob(0, f[0]=1.0, M)
  This multiplies the frequencies to get the coefficients of the generationg funtion
  for the Staden algorithm
 *************************************************************************/
{
	register double	*f2,*freq = M->freq;
	register long	s,s2,r,min,**S=M->score;
	double 	*f0;

	if(k == M->K){ 
		for(s=M->cmax[k], f[s+1]=0.0; s>=M->cmin[k]; s--){
			f[s] = f[s] + f[s+1];
		}
	} else {
	    k++;
	    NEW(f0,(M->cmax[k] - M->cmin[k] + 4),double); 
	    f2 = f0 - M->cmin[k]; min = M->cmin[k];
	    /*****
fprintf(stderr,"k = %d; cmin = %d; cmax = %d\n",k,M->cmin[k],M->cmax[k]);
/*****/
	    for(s=M->cmin[k-1]; s <= M->cmax[k-1]; s++){
		if(f[s] > 0.0){
		    for(r = M->nlet; r >= 0; r--){
			s2 = s + S[k][r];
			if(s2 >= min) f2[s2] += freq[r]*f[s];
			/*****
if(s2 > M->cmax[k])
	fprintf(stderr,"s2 = %d > max: s = %d; S[k][r] =%d; r=%d\n",
		s2,s,S[k][r],r);
/*****/
		    }
		}
	    }
	    if(M->fx != NULL) free(M->fx);
	    M->fx = f0;
	    smatrix_prob(k, f2, M);
	    if(k==M->K){M->f=f2; M->f0=f0; }
	} 
}

void	smatrix_error(char *s){fprintf(stderr,"SMatrix: %s\n",s);exit(1);}

