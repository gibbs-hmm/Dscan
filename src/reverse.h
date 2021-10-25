#ifndef _REV_
#define _REV_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "stdinc.h"

#define SEQ_ALLOC_SIZE 65536

/* #define TRUE	1 */
/* #define FALSE	0 */
/* #define	NEW(x,n,t)      (( (x=(t*) calloc(n,sizeof(t)))==NULL) ? \ */
/* 			 (t*) (intptr_t) (fprintf(stderr,"Out of Memory."),exit(1),0):x) */
#define print_error(str) for(fprintf(stderr,"%s\n",str); TRUE; exit(1))
#define Boolean	char

void rev_compl();

#endif
