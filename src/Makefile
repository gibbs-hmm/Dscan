# makefile for dscan
#

CC = cc

# gcc options
CFLAGS = -O3

INCDIRS = -I. 
LIBS = -lm 
STDINC = stdinc.h residues.h

SOBJS = smain.o afnio.o histogram.o probability.o sequence.o smatrix.o \
	smodel.o mheap.o sites.o alphabet.o karlin.o dheap.o scan.o \
	seqset.o random.o olist.o rev.o pam.o

SINC = afnio.h histogram.h probability.h sequence.h smatrix.h smodel.h \
	mheap.h sites.h residues.h alphabet.h karlin.h dheap.h scan.h \
	seqset.h random.h olist.h reverse.h pam.h


.c.o: 
	${CC} ${CFLAGS} -c $*.c


scan: ${SINC} ${SOBJS}
	${CC} ${CFLAGS} -o ../bin/dscan ${SOBJS} ${LIBS}

gibbs.o: ${INC} 
mtfsamp.o: ${INC} 
sitesamp.o: ${INC} 
columns.o: ${INC} 
options.o: ${INC} 
init.o: ${INC} 
map.o: ${INC} 
shift.o: ${INC} 
wilcoxon.o: ${INC} 
segibbs.o: ${INC} 
fmodel.o: fmodel.h $(INC)
betaprior.o: betaprior.h $(INC)

afnio.o: afnio.h $(STDINC)
alphabet.o: alphabet.h $(STDINC)
dheap.o: dheap.h $(STDINC)
gblast.o: gblast.h random.h alphabet.h dheap.h sequence.h mlist.h \
	afnio.h 
histogram.o: histogram.h $(STDINC)
karlin.o: alphabet.h karlin.h $(STDINC)
model.o: model.h alphabet.h $(STDINC)
mlist.o: mlist.h $(STDINC)
mheap.o: dheap.h mheap.h $(STDINC)
olist.o: olist.h $(STDINC)
order.o: order.h $(STDINC)
pairaln.o: pairaln.h alphabet.h sequence.h mheap.h $(STDINC)
probability.o: probability.h $(STDINC)
purge.o: $(PINC) $(STDINC)
random.o: random.h 
seqset.o: alphabet.h sequence.h seqset.h $(STDINC)
sequence.o: alphabet.h sequence.h dheap.h reverse.h $(STDINC)
sites.o: sites.h random.h dheap.h olist.h seqset.h $(STDINC)
smatrix.o: probability.h alphabet.h smatrix.h $(STDINC)
smodel.o: smatrix.h smodel.h alphabet.h $(STDINC)
swaln.o: swaln.h alphabet.h sequence.h $(STDINC)
rev.o: reverse.h $(STDINC)
smain.o: scan.h random.h reverse.h $(STDINC)
scan.o:  scan.h $(STDINC)
pam.o: pam.h $(STDINC) 

main.o: ${INC} 

clean: 
	@ \rm *.o 
	@ \rm ../bin/dscan
