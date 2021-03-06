/* mheap.c 
 *
 *
 * CVS/RCS keyword:
 * $Header: /home/thompson/cvs/dscan/mheap.c,v 2.0 2003/11/20 15:13:24 palumbo Exp $
 *
 */

#include "mheap.h"

mh_type	Mheap(long hpsz, long d)
{
	mh_type H;
	long	i;

	NEW(H,1,mheap_type); 
	NEW(H->avail,hpsz+1,long);
	H->heap=dheap(hpsz,d);
	H->maxheap=dheap(hpsz,d);
	H->hpsz = hpsz;
	H->nfree = hpsz;
	for(i=1; i<=hpsz; i++) H->avail[i] = i;
	return H;
}

mh_type	NilMheap(mh_type H)
{
	if(H==NULL) return (mh_type) NULL;
	Nildheap(H->heap);
	Nildheap(H->maxheap);
	free(H->avail);
	free(H); 
	return (mh_type) NULL;
}

long	DelMinMheap(mh_type H)
{
	long i;
	
	if((i=delminHeap(H->heap))!= 0){
		rmHeap(i,H->maxheap); 
		H->nfree++;
		H->avail[H->nfree] = i;
		return i;
	} else return 0;
}

long	DelMaxMheap(mh_type H)
{
	long i;
	
	if((i=delminHeap(H->maxheap))!= 0){
		rmHeap(i,H->heap); 
		H->nfree++;
		H->avail[H->nfree] = i;
		return i;
	} else return 0;
}

long	RmMheap(long i, mh_type H)
{
	if(rmHeap(i,H->heap) == 0) return 0; 
	rmHeap(i,H->maxheap); 
	H->nfree++;
	H->avail[H->nfree] = i;
	return i;
}

long	InsertMheap(keytyp key, mh_type H)
/* NULL == not inserted; -1 == inserted but none deleted  */
{	
	long i;

	if(H->nfree > 0){
		i = H->avail[H->nfree];
		H->nfree--;
	} else if(minkeyHeap(H->maxheap) < -key) {
		i=delminHeap(H->maxheap);
		rmHeap(i,H->heap);
	} else return 0;
	insrtHeap(i,key,H->heap);
	insrtHeap(i,-key,H->maxheap);
	return i;
}

