/* histogram.h - histogram data type */

#if !defined(HIST)
#define HIST
#include "stdinc.h"

/**********************************************************************

          min                                             max
         start   o----> <inc>            o---->o---->   end 
long |_____|_____|_____|_____|___ ... ___|_____|_____|___:_|_____|
                                      
        0     1     2     3                           nbins  over
     (under)
                                                    round up end to max
                                               last bin includes overflow

/*************************** histogram type ***************************/
typedef struct {
	char	*id;		/* identifier string for histogram */
	long	*bin0;		/* array for histogram bin0ribution */
	long	*bin;		/* bins for histogram bin0ribution */
	long	nbins;		/* number of bin0inct values in histogram */
	double	maxval,minval;	/* maximum and minimum values */
	double	max;		/* maximum value */
	double	total;		/* total for mean */
	double	total_sq;	/* total X^2 for variance */
	double	min;		/* minimum value */
	double	inc;		/* bin0ance between values in histogram */
	long	n;		/* number of variables in bin0ribution */
} histogram_type;

typedef	histogram_type	*h_type;

/**************************** Public **********************************/
h_type  Histogram(char *id,long start,long end,double inc);
long    IncdHistBin(double x,h_type H);
void	IncdHist(double x,h_type H);
void    IncdMHist(double x, long number, h_type H);
void	PutHist(FILE *fptr,long line_leng,h_type H);
void    NilHist(h_type H);
double  MeanHist(h_type H);
double  VarianceHist(h_type H);

#define indexHist(x,H)		((long)floor(((double)(x)-(H)->min)/(H)->inc))
#define IncfHist(x,H)		IncdHist((double)(x),H)
#define MinValHist(H)		((H)->minval)
#define MaxValHist(H)		((H)->maxval)
#define NumBinsHist(H)		((H)->nbins+2)

#endif

