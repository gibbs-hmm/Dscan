#if !defined (PAM)
#define PAM

#include <math.h>
#include "alphabet.h"
#include "blosum62.h"

void CalcPam1ScoreMatrix( a_type A, long *counts );

#endif

