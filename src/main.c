#include "gibbs.h"

char GIBBS_USAGE[] = GIBBS_USAGE0;

main(int argc, char *argv[])
{
	long	arg;
	ss_type	data;
	st_type	sites;
	gs_type	G;
	a_type	A=NULL;

	sites = StartSitesGibbs(argc, argv);
	data = SitesSeqSet(sites);
	A = SeqSetA(data);
	G=MkGibbs(argc, argv, sites);
	RunGibbs(stdout,G);
	NilGibbs(G);
	NilSeqSet(data);
	NilAlpha(A);
}

