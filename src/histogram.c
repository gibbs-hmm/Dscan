/* histogram.c 
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/histogram.c,v 2.0 2003/11/20 15:13:13 palumbo Exp $
 *
 */

#include "histogram.h"

h_type	Histogram(char *id,long start,long end,double inc)
/* bin0[0] = underflow; bin0[1->nbins]= values; bin0[nbins+1]=overflow */
{
	h_type	H;

	NEW(H,1,histogram_type);
	NEW(H->id,(strlen(id)+1),char); strcpy(H->id,id);
	if(end <= start || inc <= 0.0){
		fprintf(stderr,"Histogram input error\n"); exit(1); 
	}
	H->min = (double) start;
	H->inc = inc;
	H->nbins = (long) ceil((double)(end - start)/inc);
	if(H->nbins < 1) {
		fprintf(stderr,"nbins = %ld\n",H->nbins);
		fprintf(stderr,"Histogram input error\n"); exit(1); 
	}
	H->max = (double) start + (inc * (double)(H->nbins));
	H->maxval = - DBL_MAX; H->minval = DBL_MAX;
	H->n = 0; H->total = 0.0; H->total_sq = 0.0;
	NEW(H->bin0,H->nbins+3,long);
	H->bin=(H->bin0+1);
	return(H);
}

void	NilHist(h_type H){ free(H->id); free(H->bin0); free(H); }

void	PutHist(FILE *fptr,long line_leng,h_type H)
/************************************************************************
	note: Var(X) = E(X*X) - E(X)*E(X).
 ************************************************************************/
{
	long	start,total,sum,i,j,n,m, max = 0,mult_factor;
	double	val;
	
	total = 0;
	for(i=0;i<=H->nbins+1;i++) {
		total += H->bin0[i];
		max = H->bin0[i] > max ? H->bin0[i] : max;
	}
	if(max <= line_leng) mult_factor = 1;
	else mult_factor = (max/line_leng) + 1;
	fprintf(fptr,"Distribution of %s:\n",H->id);
	if(total == 0) { fprintf(fptr,"empty.\n"); return; }

	if(mult_factor==1) fprintf(fptr, " '=' is %ld count.\n", mult_factor);
	else fprintf(fptr, "  ( '=' is %ld counts. )\n", mult_factor);

	sum = 0;
	if((n=H->bin0[0]) > 0){
		sum = n;
		fprintf(fptr,"\n<%7.2f : %-8ld |",H->min,n);
		m = n/mult_factor; for(j=0;j<m;j++) fprintf(fptr,"=");
		start = 1;
		val=H->min;
	} else {
	   for(start=0,val=H->min; H->bin0[start]==0 ;start++){
	      val+=H->inc;
	      if(start >= H->nbins) {
		fprintf(stderr,"Histogram: this should not happen\n");
		exit(1); 
	      }
	   }
	   val-=H->inc;
	}

	for(i=start; i<=H->nbins; val+= H->inc, sum+=H->bin0[i],i++) {
	   n = H->bin0[i];
	   if(sum==total) { 
		fprintf(fptr,"\n%8.2f : %-8ld |",val,n); 
		m = n/mult_factor; for(j=0;j<m;j++) fprintf(fptr,"=");
		break; 
	   } else fprintf(fptr,"\n%8.2f : %-8ld |",val,n);
	   m = n/mult_factor;  for(j=0;j<m;j++) fprintf(fptr,"=");
	}
	if((n=H->bin0[H->nbins+1]) > 0){
		fprintf(fptr,"\n>=%6.2f : %-8ld |",H->max,n); 
		m = n/mult_factor; for(j=0;j<m;j++) fprintf(fptr,"=");
	}
	fprintf(fptr,"\n%8s : %-8ld\n\n","total",total);
   if(H->n > 1){
	fprintf(fptr,"    mean = %g\n",MeanHist(H));
	fprintf(fptr,"   stdev = %g\n",sqrt(VarianceHist(H)));
	fprintf(fptr,"   range = %g .. %g\n\n",H->minval,H->maxval);
   }
}

double	MeanHist(h_type H) { return H->total/(double)H->n; }

double	VarianceHist(h_type H)
{
	double	m;
	m = H->total/(double)H->n;
	return ((H->total_sq/(double)H->n) - m*m); 
}

void	IncdMHist(double x, long number, h_type H)
{
	long	i;
	if(x < H->min) H->bin0[0] += number;
	else if(x < H->max) { 
		i= (long)floor((x-H->min)/H->inc); 
		H->bin[i] += number; 
	} else H->bin[H->nbins]+=number;
	H->total += number*x; H->n += number;
	H->total_sq += number*x*x;
	H->minval = MIN(double,x,H->minval);
	H->maxval = MAX(double,x,H->maxval);
}

void	IncdHist(double x,h_type H)
{
	long	i;

	H->total += x; H->n++;
	H->total_sq += x*x;
	H->minval = MIN(double,x,H->minval);
	H->maxval = MAX(double,x,H->maxval);
	if(x < H->min){ H->bin0[0]++; }	/*** underflow ***/
	else if(x < H->max) { 		/*** min ... < max ***/
		i= (long)floor((x-H->min)/H->inc); 
		H->bin[i]++; 
	} else /** x>=H->max **/ { H->bin[H->nbins]++; } 
}

long	IncdHistBin(double x,h_type H)
/** Increment histogram and return the number of the bin incremented **/
{
	long	i;

	H->total += x; H->n++;
	H->total_sq += x*x;
	H->minval = MIN(double,x,H->minval);
	H->maxval = MAX(double,x,H->maxval);
	if(x < H->min){ H->bin0[0]++; return 0; }
	else if(x < H->max) { 
		i= (long)floor((x-H->min)/H->inc); 
		/* printf("x=%g; H->bin[%d]=%d\n",x,i,H->bin[i]);/**/
		H->bin[i]++; 
		return (i+1);
	} else { H->bin[H->nbins]++; return (H->nbins+1); }
}

