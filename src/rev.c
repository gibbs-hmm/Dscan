/* rev.c 
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/rev.c,v 2.0 2003/11/20 15:13:50 palumbo Exp $
 *
 */


#include "reverse.h"

        
void rev_compl(char *seqfile, char *seqtmp)
{	
	FILE 	    *fptr;
	FILE	    *fptr1;
	char	    id[128],rev_id[128],*seq,*rev;
	long	    length,i,leng;
	char 	    c='\n';
 	long	    MAX_IN_SEQS=300000;
 	int 	    title1,title2;
	int         nBlkSize;

 	if((fptr = fopen(seqfile,"r")) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",seqfile);
                print_error("File does not exist.");
        }
        if((fptr1 = fopen(seqtmp,"wr")) == NULL) {
                fprintf(stderr,"Could not create file \"%s\"\n",seqtmp);
                print_error("File creation error.");
        }
        
       
       if((c=fgetc(fptr))==EOF)
           print_error("File not in fasta format!");
           
       
       do{   
         if(c == '>')
              id[0]='>';
              
         for(i=1; (c=fgetc(fptr))!=EOF; i++){ 
	   
	   if(c=='\n') { 
	     if(i<73) {	      
	       id[i]= '\n'; 
	       title1=i;
	     }
	     break;
	   }
	   else if(i<72) {
	     id[i] = c;
	     title1=i;
	   }
	   else if(i==72) {
	     id[i] = '\n';
	     title1=i;
	   }
	 }

	for(i=0;i<title1;i++)
		rev_id[i]=id[i];
	rev_id[i]=' ';
	rev_id[i+1]='R';
	rev_id[i+2]='\n';
	title2=title1+2;
		
	for (i=0;i<=title1;i++)
		fputc(id[i],fptr1);

	/*	NEW(seq,MAX_IN_SEQS+1, char);
		NEW(rev,MAX_IN_SEQS+1, char); */

	nBlkSize = SEQ_ALLOC_SIZE+1;
	NEW(seq,nBlkSize, char);
	NEW(rev,nBlkSize, char); 	

        for(length=0,leng=0; (c=fgetc(fptr))!=EOF;){ 
                if(c == '>') {
		   ungetc(c,fptr);
		   break; 
		}
		else 
		   if(isalpha(c)) {
			if(islower(c))
			    c = toupper(c);
			switch(c){
			    case 'A':rev[leng++]='T'; break;
			    case 'T':rev[leng++]='A'; break;
			    case 'G':rev[leng++]='C'; break;
			    case 'C':rev[leng++]='G'; break;
			    case 'N':rev[leng++]='N'; break;
			}
			seq[length++] = c;
			if( length >= nBlkSize )
			  {
			    nBlkSize += SEQ_ALLOC_SIZE;
			    seq = (char *) realloc( seq, nBlkSize );
			    rev = (char *) realloc( rev, nBlkSize );
			  }
		   }
        }
	
	for(i=0;i<length;i++)
	     fputc(seq[i],fptr1);
	fputc('\n',fptr1);
	for (i=0;i<title2;i++)
		fputc(rev_id[i],fptr1);
	fputc('\n',fptr1);
	do{
	   leng--;}while (!rev[leng]) ; 
	      
	for(i=leng;i>=0;i--)
	     fputc(rev[i],fptr1);
	fputc('\n',fptr1);
	free(seq);
	free(rev);
    }while((c=fgetc(fptr))!=EOF);
    fclose(fptr);
    fclose(fptr1);
    /*    free(seq);
	  free(rev); */
   /*printf("\n\n A temporary file has been created in /tmp called %s\n\n", seqtmp);*/
}
