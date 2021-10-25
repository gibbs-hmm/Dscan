/****
A Gibbs Sampler algorithm for finding multiple sites in multiple sequences 
****/
#include "gibbs.h"
unsigned long   gibbs_random_seed;

st_type	StartSitesGibbs(long argc, char *argv[])
/** get starting sites from argument strings **/
{
   st_type S;
   ss_type data;
   a_type  A=NULL;
   FILE	   *fptr;
   long    l;
   long     arg,t,i,n,ntyps;
   long     expect[MAX_NO_TYPESITES+2],sitelen[MAX_NO_TYPESITES+2];
   Boolean wilcox=FALSE,xnu=TRUE,shuffle=FALSE;
   Boolean msamp,interactive=FALSE;
   unsigned long seed=18364592;
   extern char GIBBS_USAGE[];

   if(argc < 3) print_error(GIBBS_USAGE);
   ntyps = ParseIntegers(argv[2], sitelen, GIBBS_USAGE);
   if(ntyps > MAX_NO_TYPESITES) print_error("too many element types");
   if(argc == 3 || argv[3][0] == '-'){
	arg = 3; msamp = FALSE;
   } else {	/** i.e., motif sampler **/
   	n = ParseIntegers(argv[3], expect, GIBBS_USAGE);
   	if(ntyps != n) print_error("inconsistent numbers of elements");
	arg = 4; msamp = TRUE;
   }
   for(  ; arg < argc; arg++){
       if(argv[arg][0] != '-') print_error(GIBBS_USAGE);
	switch(argv[arg][1]) {
	   case 'I': interactive=TRUE; break;
	   case 'n': A = MkAlpha("NACGT",NULL); xnu=FALSE; break;
	   case 'r': shuffle=TRUE; break;
	   case 's':
	  	  if(sscanf(argv[arg],"-s%ld", &l) == 1 ){
			seed = l;
		  } else  print_error(GIBBS_USAGE);
	   	break;
	   case 'w': wilcox = TRUE; 
		break;
	   case 'x': xnu = FALSE; break;
	   case 'k': print_error(GIBBS_USAGE); break;
	      	      
	   default: /** ignore **/ break;
	}
   }
   for(t=1; t<=ntyps; t++){
       if(sitelen[t] < 2 || sitelen[t] > 120)
                        print_error("element length out of range");
   }
   if(A==NULL) A=MkAlphabet(AMINO_ACIDS,PROT_BLOSUM62,SREL26_BLSM62);
   if(seed == 18364592) seed = (unsigned long) time(NULL)/2;
   gibbs_random_seed = seed;
   sRandom((unsigned long)seed);
   if(xnu) {
        if((data=MkXnuSeqSet(argv[1],A)) == NULL) 
		print_error("can not xnu sequences");
   } else data = SeqSet(argv[1],A);
   if(shuffle) ShuffleSeqSet(data);
   if(wilcox){
        for(t=1; t<=ntyps; t++) expect[t] *=2;
        fptr = tmpfile();
        PutSeqSetEs(fptr,data); ShuffleSeqSet(data); 
	PutSeqSetEs(fptr,data); NilSeqSet(data);
        rewind(fptr);
        data = SeqSet_fptr(fptr,A);
        fclose(fptr);
    }
    if(msamp){
    	S=MkSites(ntyps, sitelen, data);
        for(t=1; t<= ntyps; t++){
            if(!NRandomSites(t, expect[t], 500, S)){
                PutSites(stderr, t, S, NULL,NULL);
                print_error("too many sites; try a smaller number");
            }
	}
    } else {
      if(interactive){
	S = ReadRandomSites(ntyps, sitelen, data);
      } else {
    	S=MkSites(ntyps, sitelen, data);
	for(t = 1; t <= ntyps; t++){
            for(n = 1; n <= NSeqsSeqSet(data); n++) {
        	if(AddRandomSite(t,n,100000,S) == 0){
                      print_error("too many sites; try a smaller number");
		}
	    }
	}
      }
    }
    return S;
}

Boolean	OptionsGibbs(long argc, char *argv[], gs_type G)
{
   long     arg,i,t,w;
   float   f;
   double  d,T;
   Boolean b,msamp,flag=TRUE;
   Boolean wilcox=FALSE;
   extern char GIBBS_USAGE[];

   if(argc < 3) return FALSE;
   G->name = String(argv[1]);
   G->seed = gibbs_random_seed;
   if(argc == 3 || argv[3][0] == '-'){ 
	G->start = ArchiveSites(G->sites);
	arg = 3; msamp = FALSE; 
   } else {  /** use motif sampling **/
	NEW(G->expect,G->ntyps + 3, long);
   	i = ParseIntegers(argv[3], G->expect, GIBBS_USAGE);
   	if(G->ntyps != i) print_error("inconsistent numbers of elements");
	arg = 4; msamp = TRUE; 
   }
   for(; arg < argc; arg++){
       if(argv[arg][0] != '-') return FALSE;
	switch(argv[arg][1]) {
           case 'C':    /* Cutoff for print out */
                if(sscanf(argv[arg],"-C%lf", &d) != 1) return FALSE;
                if(d > 1.0 || d < 0.0) return FALSE;
		G->readcutoff = d;
                break;
	   case 'c':	/* number of cycles per shift */
	  	if(sscanf(argv[arg],"-c%d", &i) != 1) return FALSE;
		G->ncycles = MAX(long,1,i);
	   	break;
	   case 'd': G->fragment=FALSE; break;
	   case 'D': flag=FALSE; break;
	   case 'f': G->sfptr=open_file(NameSeqSet(G->data),".sn","w");
		break;
	   case 'I': /** ignore **/ break;
	   case 'L': if(sscanf(argv[arg],"-L%d", &i) != 1) return FALSE;
	        G->limit = i; G->stop=0;
		if(i < 1 || i >= 10000) print_error("stop out of range");
	        for(t=1;t<=G->ntyps;t++)G->stop+=NSeqsSeqSet(G->data)*i;
	   	break;
	   case 'M': G->mfptr=open_file(NameSeqSet(G->data),".mtf","w");
		break;
	   case 'm': if(sscanf(argv[arg],"-m%d", &i)!= 1) return FALSE;
		G->nconverge=MAX(long,1,i);
	   	break;
	   case 'n': /** do nothing **/ break;
	   case 'o': G->use_order = TRUE; break;
	   case 'p':	/* number of pseudo counts  for model */
	  	if(sscanf(argv[arg],"-p%lf", &d) != 1) return FALSE;
        	G->pseudo = MIN(double,10.0,d);
        	G->pseudo = MAX(double,0.0000001,G->pseudo);
	   	break;
	   case 'q':	/* number of pseudo counts  for model */
	  	if(sscanf(argv[arg],"-q%lf", &d) != 1) return FALSE;
		G->qseudo=d;
	   	break;
	   case 'R':
	  	if(sscanf(argv[arg],"-R%d",&i)!=1 || i<100) return FALSE;
		G->nread=i;
	   	break;
	   case 'r': /** do nothing **/ break;
	   case 's': /** do nothing **/ break;
	   case 'T': 
	   	if(sscanf(argv[arg],"-T%d",&i) != 1) return FALSE;
		if(i < 0 || i > 10) return FALSE;
		G->test[i] = TRUE; 
		break;
	   case 't': if(sscanf(argv[arg],"-t%d",&i) != 1) return FALSE;
		G->nruns=i;
	   	break;
	   case 'v': G->verbose=TRUE; break;
	   case 'W': 
	  	if(sscanf(argv[arg],"-W%lf", &d) != 1) return FALSE;
		if(d >= 1.0 || d <= 0.0) return FALSE;
		G->weight=d;
	   	break;
	   case 'w': 
   		wilcox=TRUE;
		if(NSeqsSeqSet(G->data) % 2 != 0)
			print_error("Wilcoxon: control set not present");
		G->wilcoxon=NSeqsSeqSet(G->data)/2; 
		break;
	   case 'x': /** do nothing **/ break;
	   default: return FALSE;
	}
    }
    if(wilcox) {
        for(t = 1; t <= G->ntyps; t++) G->expect[t] *= 2;
    }
    if(msamp){
	if(G->move) print_error(GIBBS_USAGE);
        NEW(G->best_p,G->ntyps+2,double);
        NEW(G->p,G->ntyps+2,double); NEW(G->map_p,G->ntyps+2,double);
        NEWPP(G->readprob,G->ntyps+2,double);
        NEW(G->sumprob,G->ntyps+2,double);
        for(t = 1; t <= G->ntyps; t++){
           NEWP(G->readprob[t],NSeqsSeqSet(G->data)+2,double);
           for(i = 1; i <= NSeqsSeqSet(G->data); i++) {
                NEW(G->readprob[t][i],SqLenSeqSet(i,G->data)+2,double);
           }
	}
	NEW(G->prior,G->ntyps+2,bp_type);
        if(flag) {
	   for(arg = 0; arg < argc; arg++) printf("%s ",argv[arg]); printf("\n");
	   PutSeqSetPIDs(stdout,G->data);
           fprintf(stdout,"\tseed: %d\n", G->seed);
	   fprintf(stdout, "         mean (+/-2sd):\n");
	}
	for(t = 1; t <= G->ntyps; t++){
           for(T=0.0, i = 1; i <= NSeqsSeqSet(G->data); i++) {
                w = SqLenSeqSet(i,G->data) - SiteLen(t,G->sites) + 1;
                if(w > 0) T += (double) w;
           }
	   if(flag){
             if(G->ntyps <=26) fprintf(stdout, "motif %c: ", (t+'A' - 1));
             else fprintf(stdout, "motif %d: ", t);
	   }
	   T -= SiteLen(t,G->sites)*G->expect[t]; 
           G->prior[t] = MkBPrior(G->expect[t], G->weight, T);
           if(flag) PutBPrior(stdout, G->prior[t]);
	}
    } else if(G->move){
	i= TotalSites(1,G->sites);
	for(t = 2; t <= G->ntyps; t++){
	  if(i != TotalSites(t,G->sites)){
	     print_error("The -a option requires the same # sites for all motifs");
	  }
	}
        if(flag) {
	   for(arg = 0; arg < argc; arg++) printf("%s ",argv[arg]); printf("\n");
	   PutSeqSetPIDs(stdout,G->data);
           fprintf(stdout,"\tseed: %d\n", G->seed);
        }
    }
    return TRUE;
}

