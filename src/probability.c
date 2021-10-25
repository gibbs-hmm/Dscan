/* probability.c 
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/probability.c,v 2.0 2003/11/20 15:13:42 palumbo Exp $
 *
 */

#include "probability.h"

double  bico(long N,long k)
{ return floor(0.5+exp(lnfact(N)-lnfact(k)-lnfact(N-k))); }

double  lnbico(register long N, register long k)
{ return lnfact(N)-lnfact(k)-lnfact(N-k); }

double	lnfact(long n)
/* static variables are guaranteed to be initialized to zero */
{
	static double lnft[101];

	if (n <= 1) return 0.0;
	if (n <= 100) return lnft[n] ? lnft[n] : (lnft[n]=lgamma(n+1.0));
	else return lgamma(n+1.0);
}


