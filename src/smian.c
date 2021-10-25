#include "scan.h"
#include "random.h"

#define	USAGE_START	"DScan 1.0.3 12/22/97\n\
   USAGE: dscan database snfile [options]\n\
   snfile = file with aligned segments generated using gibbs or asset\n\
   options:\n\
     -c         - create sequence file for output\n\
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
     -x         - mask out sequence regions matching motif(s)\n\n"

/**************************** Global Variables ******************************/
void	main(long argc,char *argv[])
{ 
	long	i,j,k,length,arg,number,hpsz=2000;
	long	min_segs=1,max_segs=10,*counts,time1,total;
	long	MAX_IN_SEQS=300000;
	/*	unsigned short	*nsize; */
	long	*nsize;          /* BT 9/26/97 */
	char	c,method='g';	/** default: Gribskov's method **/
	float	expect=0.01,pseudo=0.1,singlePval=0.0;
	a_type	A;
        FILE    *fptr;
	e_type	E;
	scn_typ F;
	Boolean seqfile=FALSE,order=FALSE,repeats=FALSE;
	Boolean	mask=FALSE,shuffle=FALSE,neg_mask=FALSE,create_mtf=FALSE;
	unsigned long seed=18364592;
	char   *alphabet = AMINO_ACIDS;
	char   *Rmat = PROT_BLOSUM62;
	short  nFlag = FALSE;     /* BT 10/17/97 */
	int    nCount = 0;

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
		if(sscanf(argv[arg],"-h%d",&hpsz)!=1)print_error(USAGE_START);
		break;
	      case 'E':
		if(sscanf(argv[arg],"-E%f", &singlePval)!= 1)
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
		if(sscanf(argv[arg],"-e%f", &expect)!= 1)
		  print_error(USAGE_START);
		break;
	      case 'M': 
		if(sscanf(argv[arg],"-M%d",&max_segs)!=1)
		  print_error(USAGE_START);
		break;
	      case 'm': 
		if(sscanf(argv[arg],"-m%d",&min_segs)!=1)
		  print_error(USAGE_START);
		break;
	      case 'P': method = 'r'; break;
	      case 'g': method = 'c'; break;          /* alternative Gribsov method */
	      case 'i': method = 'i'; break;          /* identity method */
	      case 'p':
		if(sscanf(argv[4],"-p%f",&pseudo)!=1)print_error(USAGE_START);
		break;
	      case 'o': order = TRUE; break;
	      case 'r': repeats = TRUE; break;
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
		method = 'i';               /* BT 10/14/97 */
		break;

	      case 'N':                    /* BT 10/17/97 */
		if(sscanf(argv[arg],"-N%d",&nCount)!=1)
		  print_error(USAGE_START);
		  nFlag = TRUE;
		  expect = FLT_MAX;
		break;
			     
		default: print_error(USAGE_START);
	    }
	}

	A=MkAlpha(alphabet, Rmat);  /****/        /* BT 9/19/97 */

	if(seed == 18364592)  seed = (unsigned long) time(NULL)/2;
	sRandom(seed);
	number = GetFastaInfo(argv[1], MAX_IN_SEQS, &counts, &nsize, A);
	if(!mask){
	   printf( "DScan 1.0.3 12/22/97\n" );                   /* BT 12/22/97 */
	   for(arg = 0; arg < argc; arg++) printf("%s ",argv[arg]); printf("\n");
	   for(total=0, i=1; i<=nAlpha(A); i++) total += counts[i];
	   /*	   fprintf(stdout,"%c: %d\n", AlphaChar(i,A),counts[i]); */  /* BT 10/14/97 */
	   for(i=1; i<=nAlpha(A); i++){
		fprintf(stdout,"%c: %d (%f)\n",
			AlphaChar(i,A),counts[i],
			(double) counts[i]/(double) total);
	   }
	   fprintf(stdout,"\n");
	}
        if((fptr = fopen(argv[1],"r")) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",argv[1]);
                print_error("File does not exist!");
        }
	F = MkScan(pseudo,argv[2],counts,A,hpsz,method, nFlag, nCount);
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
		printf("%d sequences scanned\n",number);
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
}

