#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include "scan.h"
#include "random.h"
#include "reverse.h"
#include "pam.h"

/*  smain.c 
 *
 * last modified:  01-30-04 bt
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/smain.c,v 2.3 2004/05/05 13:12:04 thompson Exp $
 *
 */

/* mjp: 10-22-03: nucleotide order is defined as (in residue.h):  XCGAT  */


char   revStr[] = "$Revision: 2.3 $";

/* #define	USAGE_START	"dscan 1.0.11 04/29/2002\n\ */
#define	USAGE_START	"dscan $Revision: 2.3 $\n\
   USAGE: dscan database snfile [options]\n\
   database = sequences to be searched, in fasta format\n\
   snfile = file with aligned segments generated using Gibbs\n\
   options:\n\
     -c         - create sequence file for output\n\
     -C         - create motif listing\n\
     -e<float>  - maximum E-value detected\n\
     -E<float>  - minimum -log10(P-value) required for each motif\n\
     -h<int>    - size of heap for saving sequences\n\
     -M<int>    - maximum # repeats per sequence (use with -r)\n\
     -m<int>    - minimum # repeats per sequence (use with -r)\n\
     -n         - use DNA alphabet and Identity matrix\n\
     -N<int>    - print top <int> values. Ignore -e.\n\
     -o         - require motifs to be in correct order\n\
     -P         - use product multinomial model instead of Gribskov method\n\
     -g         - use alternative Gribskov method\n\
     -p<float>  - pseudo counts for product multinomial model\n\
     -r         - scan for repeats \n\
     -S         - shuffle input sequences\n\
     -s<int>    - seed for random number generator\n\
     -X         - mask out sequence regions NOT matching motif(s)\n\
     -x         - mask out sequence regions matching motif(s)\n\
     -R         - reverse complement search of database file \n\
                  (a default file 'seqXXXXXXX' is created in '/tmp')\n\
     -a         - use PAM1 DNA matrix\n\
     -d<int>[,<int>] - num 'sites' for each *freq* model in snfile; comma\n\
		  separated list (100)\n\
     -q<int>    - maximum number of input sequences\n\
     -v         - don't allow overlaping motifs when using multiple models\n\
     -B         - alternate Bonferroni adjustment for the size of the database (single model only)\n\
     -h         - this message\n\n"


/**************************** Global Variables ******************************/
int	main(long argc,char *argv[])
{ 
	long	i,j,k,length,arg,number,hpsz=2000;
	long	min_segs=1,max_segs=10,*counts,time1,total;
	long	MAX_IN_SEQS=300000;
	/*	unsigned short	*nsize; */
	long	*nsize;          /* BT 9/26/97 */
	char	c,method='g';	/** default: Gribskov's method **/
	double	expect=0.01,pseudo=0.1,singlePval=0.0;
	a_type	A;
        FILE    *fptr;
	e_type	E;
	scn_typ F;
	Boolean seqfile=FALSE,order=FALSE,repeats=FALSE,reverse=FALSE;
	Boolean	mask=FALSE,shuffle=FALSE,neg_mask=FALSE,create_mtf=FALSE;
	unsigned long seed=18364592;
	char   *alphabet = AMINO_ACIDS;
	char   *Rmat = PROT_BLOSUM62;
	short  nFlag = FALSE;     /* BT 10/17/97 */
	int    nCount = 0;
	/*char   *seqrev="/tmp/seq.tmp";*/  /* default*/
	char   seqrev[] = "/tmp/seqXXXXXX";
        int    nNoOverlap = FALSE;	
	int    adjustN = FALSE;

	int    numNumSitesArgs      = 0;      /* mjp: 10-31-03 */
	long   *freqModelNumSites = NULL;   /* mjp: 10-31-03 */

	time1=time(NULL);
	if(argc < 3) print_error(USAGE_START);
	sRandom((unsigned long) time(NULL)/2);
	/*** A=MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);  ****/
	/*** A=MkAlpha(AMINO_ACIDS,PROT_PAM40);  /****/
	for(arg = 3; arg < argc; arg++)
	  {
	    if(argv[arg][0] != '-') print_error(USAGE_START);
	    switch(argv[arg][1]) 
	      {
	      case 'C': create_mtf = TRUE; break;
	      case 'c': seqfile = TRUE; break;
	      case 'h': 
		if(sscanf(argv[arg],"-h%ld",&hpsz)!=1)print_error(USAGE_START);
		break;
	      case 'E':
		if(sscanf(argv[arg],"-E%lf", &singlePval)!= 1)   /* BT 6/23/2000 */
		  print_error(USAGE_START);
		/*		singlePval = -singlePval;
				if(singlePval > 0)
				{
				print_error("input error: E must be <= 0");
				} */
		if(singlePval < 0)          /* BT 10/14/97 */
		  {
		    print_error("input error: E must be >= 0");
		  }
		break;
	      case 'e':
		if(sscanf(argv[arg],"-e%lf", &expect)!= 1)   /* BT 6/23/2000 */
		  print_error(USAGE_START);
		break;
	      case 'M': 
		if(sscanf(argv[arg],"-M%ld",&max_segs)!=1)
		  print_error(USAGE_START);
		break;
	      case 'm': 
		if(sscanf(argv[arg],"-m%ld",&min_segs)!=1)
		  print_error(USAGE_START);
		break;
	      case 'P': method = 'r'; break;
	      case 'g': method = 'c'; break;          /* alternative Gribsov method */
	      case 'i': method = 'i'; break;          /* identity method */
	      case 'p':
		if(sscanf(argv[arg],"-p%lf",&pseudo)!=1)print_error(USAGE_START); /* BT 04/05/04 */
		break;
	      
	      case 'o': 
		order = TRUE; 
		if( adjustN )
		  {
		    fprintf( stderr, "****Warning -o and -B options are incompatable. -B ignored\n" );
		    adjustN = FALSE;
		  }
		break;
		
	      case 'r': 
		repeats = TRUE; 
		if( adjustN )
		  {
		    fprintf( stderr, "****Warning -r and -B options are incompatable. -B ignored\n" );
		    adjustN = FALSE;
		  }
		break;

	      case 'S': shuffle = TRUE; break;
	      case 's': 
		if(sscanf(argv[arg],"-s%ld",&seed)!=1)
		  print_error(USAGE_START);
		break;
	      case 'x': mask = TRUE; neg_mask = FALSE; break;
	      case 'X': mask = TRUE; neg_mask = TRUE; break;
		
	      case 'n':                    /* BT 9/19/97 */
		alphabet = NUCLEOTIDES;
		Rmat = ID_MATRIX;
		if( method != 'r' )
		  method = 'i';               /* BT 04/12/04 */
		break;

	      case 'a':                    /* BT 8/3/2000 */
		alphabet = NUCLEOTIDES;
		Rmat = ID_MATRIX;
		if( method != 'r' )
		  method = 'm';               /* BT 04/12/04 */
		break;


	      case 'd':                    /* mjp 10/31/2003 */
		/* read possibly multiple 'num sites'; 1 for each matrix
		 * model in snfile */
		/* get mem for array that will hold num sites values.
		 * NEW() init's elements to 0. 
		 */
		{ 
		  char   argBuf[255], *tokP;
		  int    i;

		  /* NEW() uses calloc, so elements are init'ed to 0. */
		  NEW(freqModelNumSites,MAX_NUM_MODELS,long);
		  /* grab the entire arg string */
		  if ( sscanf(argv[arg],"-d%s",argBuf) !=1 )
		    print_error(USAGE_START);

		  /* printf("got -d: %s\n",argBuf); */
		  /* now parse the entire arg string */
		  tokP = strtok(argBuf,",");
		  numNumSitesArgs = 0;
		  while ( tokP != NULL ) {
		    /* printf("token: %s\n",tokP); */
		    if ( sscanf(tokP,"%ld",
			       &freqModelNumSites[numNumSitesArgs++]) != 1 )
		      print_error(USAGE_START);
		    tokP = strtok(NULL,",");
		  }
#if 0
		  printf("done reading -d args; numNumSitesArgs = %d\n",
			 numNumSitesArgs);
		  for ( i = 0; i < numNumSitesArgs; i++ )
		    printf("  %ld",freqModelNumSites[i]);
		  printf("\n");
#endif 
		}
		break; /* case 'd' */


	      case 'N':                    /* BT 10/17/97 */
		if(sscanf(argv[arg],"-N%d",&nCount)!=1)
		  print_error(USAGE_START);
		nFlag = TRUE;
	        expect = DBL_MAX;
		break;
		
	      case 'R':			   /* BT 02/25/98 */
	         reverse=TRUE;
		 /* seqrev = (char *) tempnam( "/tmp", "seq" ); */
		 int fd = mkstemp(seqrev);
		 close(fd);
	         rev_compl(argv[1], seqrev); 
	         break;

	      case 'q':
		if(sscanf(argv[arg],"-q%ld",&MAX_IN_SEQS)!=1)
		  print_error(USAGE_START);
		break;

	      case 'v':                     /* BT 10/4/2000 */ 
                nNoOverlap = TRUE;
		break;

	      case 'B':                     /* BT 11/6/2001 */ 
                adjustN = TRUE;
		if( order )
		  {
		    fprintf( stderr, "****Warning -o and -B options are incompatable. -B ignored\n" );
		    adjustN = FALSE;
		  }
		if( repeats )
		  {
		    fprintf( stderr, "****Warning -r and -B options are incompatable. -B ignored\n" );
		    adjustN = FALSE;
		  }
		break;
			     
	      default: 
		print_error(USAGE_START);
	    }
	}

	A=MkAlpha(alphabet, Rmat);  /****/        /* BT 9/19/97 */

	if(seed == 18364592)  seed = (unsigned long) time(NULL)/2;
	sRandom(seed);
	if(reverse==TRUE)
	    number=GetFastaInfo(seqrev, MAX_IN_SEQS, &counts, &nsize, A);
	else    
	    number=GetFastaInfo(argv[1], MAX_IN_SEQS, &counts, &nsize, A);

	if( method == 'm' )    /* BT 8/3/2000 */
	  CalcPam1ScoreMatrix( A, counts );

	if(!mask){
	  /* printf( "dscan 1.0.11 04/29/2002\n" );  */
	  printf("dscan %s\n",revStr);
	   for(arg = 0; arg < argc; arg++) 
	     printf("%s ",argv[arg]); 
	   printf("\n\n");
	   for(total=0, i=1; i<=nAlpha(A); i++) total += counts[i];
	   /* fprintf(stdout,"%c: %d\n", AlphaChar(i,A),counts[i]); */  /* BT 10/14/97 */
	   for(i=1; i<=nAlpha(A); i++){
		fprintf(stdout,"%c: %ld (%f)\n",
			AlphaChar(i,A),counts[i],
			(double) counts[i]/(double) total);
	   }
	   fprintf(stdout,"\n");
	}

	if (reverse==TRUE){
          if((fptr = fopen(seqrev,"r")) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",seqrev);
                print_error("File does not exist!");
          }
        }
        else {
          if((fptr = fopen(argv[1],"r")) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",argv[1]);
                print_error("File does not exist!");
          }
        }

	F = MkScan(pseudo,argv[2],counts,A,hpsz,method, nFlag, nCount, nNoOverlap, adjustN,freqModelNumSites,numNumSitesArgs);

	if( F->N > 1 && F->adjustN )
	  {
	    fprintf( stderr, "****Warning -B option can only be used with single models, -B ignored.\n" );
	    F->adjustN = FALSE;
	  }

	OpenSeqFileScan(seqfile,F); 
	SetMaxEvalScan((double)expect, F);
	SetSinglePvalScan(singlePval,F);
/*** NEW ***/
	if(create_mtf) { OrderScanMtf(argv[1],F); }
/*** NEW ***/
	else if(mask) {
	   if(neg_mask) UseNegMaskScan(F);
	   if(order || repeats) {
		fprintf(stderr,"Masking not yet implemented for this option\n");
	   } else {
		MaskScan(fptr,number, nsize,F);
	   }
	} else {
	  if(order) {
		OrderScanScan(fptr,number, nsize,F); 
		PutScan(stdout, F);
	  } else if(repeats) {
		if(shuffle) ShuffleSegsScan(F);
		SetMinSegsRScan(min_segs,F);
		SetMaxSegsRScan(max_segs,F);
		printf("%ld sequences scanned\n",number);
		ScanRScan(fptr,number, nsize, F);
		PutRScan(stdout, F);
	  } else {
		ScanScan(fptr,number, nsize,F); 
		PutScan(stdout, F);
	  }
	  fprintf(stdout,"\ttime: %ld seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
	}
        fclose(fptr); NilAlpha(A); NilScan(F); free(counts);
	free(nsize);
	if( reverse )
	  {
	    remove( seqrev );
	    /* free( seqrev ); */
	  }
}

