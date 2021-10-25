/* purge.h - codes and constants for purge program. */
#if !defined (PURGE)
#define PURGE
#include <time.h>
#include "afnio.h"
#include "block.h"
#include "residues.h"
#include "seqset.h"
#include "pairaln.h"
#include "dheap.h"
#include "gblast.h"

Boolean	RmHomologs(long cutoff, char method, long minimum, Boolean query, 
	ss_type P);

        /********* N   A   C   T   G *******/
#define DNA_MTRX "-4  -4  -4  -4  -4 \
                  -4   5  -4  -4  -4 \
                  -4  -4   5  -4  -4 \
                  -4  -4  -4   5  -4 \
                  -4  -4  -4  -4   5 "

#define PURGE_USAGE	"\nusage: purge file score <options>\n\
   options:\n\
     [-b]    - use blast heuristic method (default)\n\
     [-e]    - use an exhaustive method\n\
     [-q]    - keep first sequence in the set\n\
     [-x]    - don't use xnu to mask low complexity regions\n\
\n"

#endif

