/* pam.c - calculate PAM1 scoaring matrix for DNA 
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/pam.c,v 2.0 2003/11/20 15:13:39 palumbo Exp $
 *
 */



#include "pam.h"

void CalcPam1ScoreMatrix( a_type A, long *counts )
{
  double total;
  int    i;
  int    j;
  double *freq;
  double log2 = log(2.0);

  NEW( freq, nAlpha(A) + 1, double );
  
  for(total=0, i=1; i<=nAlpha(A); i++) 
    total += counts[i];

  for(i=1; i<=nAlpha(A); i++)
    {
      freq[i] = ((double) counts[i]) / total;
    }

  for(i=1; i<=nAlpha(A); i++)
    {
      for(j=1; j<=nAlpha(A); j++)
	{
	  pam1Score[i][j] = log( (pam1[i][j] * pam1[j][i]) / (freq[i]*freq[j]) ) / (2 * log2);
	  pam1ScoreL[i][j] = pow( 2.0, pam1Score[i][j] );
	}
    }

  free( freq );
}


