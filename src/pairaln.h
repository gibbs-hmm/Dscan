/* pairaln.h - pairwise alignment methods . */
#if !defined(PAIRALN)
#define PAIRALN
#include "stdinc.h"
#include "alphabet.h"
#include "sequence.h"
#include "mheap.h"

/******************************** private **********************************/
Boolean relate_seq(register char *seq1, register char *seq2,
        register char **R, register long n, register long cutoff);
long     diagonal_score_seq(register char *seq1, register char *seq2,
        register char **R, register long n);
long     get_diagonal_ends_seq(char *seq1, char *seq2, char **R, long n,
        long *begin, long *end);
long     repeat_score_seq(register char *seq, register char **R,
                register long o, register long n);
/******************************** PUBLIC **********************************/
long     PutDiagonalSeq(FILE *fptr, long offset, e_type E1, e_type E2, a_type A);
Boolean RelatedSeqs(long cutoff, e_type E1, e_type E2, a_type A);
long     AlignSeqFastp(e_type E1, e_type E2, a_type A);
Boolean RelateSeqFastp(e_type E1, e_type E2, a_type A,long score);
long     RepeatScoreSeq(e_type E, a_type A, long *offset);

#endif
