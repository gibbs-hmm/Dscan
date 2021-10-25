#if !defined(PROBABILITY)
#define PROBABILITY
#include <math.h>
#include "stdinc.h"

double lgamma(double x);
double  lnfact(long n);
double  bico(long N,long k);
double  lnbico(register long N, register long k);

/*********************************************************/

#define NATURAL_LOG_FACTOR 2.3025850929940459
#define LnBeta(a,b)             (lgamma((double)(a))+lgamma((double)(b))\
                                        - lgamma((double)((a)+(b))))

#endif

