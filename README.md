# Dscan
# dscan #

dscan searches a collection of sequences in FASTA format to find sites that match a given motif pattern.<br><br>
dscan also requires a
description of the motif model. This may be in the form of a collection of sites as output by the Gibbs or a matrix. The sites are used to create a position weight matrix (PWM). The PWM is used to sequentially search the FASTA sequences and return sites that have a statistically significant match. <br<
The entries in the matrix are
either the probabilities of finding each nucleotide or amino acid at each position or a
count of the nucleotides or amino acids at that position. The order of the columns is
ATCG for nucleotides and ACDEFGHIKLMNPQRSTVWY for amino acids<br>

See Neuwald, A., Liu, J., and Lawrence, C. 1995. Gibbs
motif sampling: Detection of bacterial outer
membrane protein repeats. Protein Sci. 4:1618-1632. https://onlinelibrary.wiley.com/doi/10.1002/pro.5560040820 for details of the method<br>

<pre><code>
*.*****.*.*..********
TTTTTGATCGTTTTCACAAAA
TATTTGCACGGCGTCACACTT
ACTGTGAGCATGGTCATATTT
TATGCAAAGGACGTCACATTA
AACGTGATCAATTTAACACCT
TATTTGAACCAGATCGCATTA
ACTTCGATACACATCACAATT
AATTTATTCCATGTCACACTT
AACGTGATCAACCCCTCAATT
TGAGTGAGCTAACTCACATTA
TCTGTAACAGAGATCACACAA
TTTGCAAGCAACATCACGAAA
AATTTGCACTGTGTCACAATT
TGCCTGACGGAGTTCACACTT
ATTGTGATTCGATTCACATTT
TGGGTTAACCACATCACAACA
TCTGTGCGGTATTTCACACCG
GTGCCGATCAACGTCTCATTT
TATGTGCGACCACTCACAAAT

dscan model consisting of a collection of sites.


*.*****.*.*..********
  7      11       0       1
  8       4       4       3
  1      13       3       2
  0       6       2      11
  0      15       4       0
  4       1       0      14
 14       1       4       0
  5       7       2       5
  3       1      12       3
  5       3       6       5
 13       2       1       3
  3       5       7       4
  5       5       3       6
  0      18       1       0
  1       0      18       0
 16       2       0       1
  0       1      18       0
 18       0       0       1
  7       6       6       0
  4      12       3       0
  7      11       0       1

The same data displayed as a count matrix.


*.*****.*.*..********
 0.37    0.58    0.00    0.05
 0.42    0.21    0.21    0.16
 0.05    0.68    0.16    0.11
 0.00    0.32    0.11    0.58
 0.00    0.79    0.21    0.00
 0.21    0.05    0.00    0.74
 0.74    0.05    0.21    0.00
 0.26    0.37    0.11    0.26
 0.16    0.05    0.63    0.16
 0.26    0.16    0.32    0.26
 0.68    0.11    0.05    0.16
 0.16    0.26    0.37    0.21
 0.26    0.26    0.16    0.32
 0.00    0.95    0.05    0.00
 0.05    0.00    0.95    0.00
 0.84    0.11    0.00    0.05
 0.00    0.05    0.95    0.00
 0.95    0.00    0.00    0.05
 0.36    0.32    0.32    0.00
 0.21    0.63    0.16    0.00
 0.37    0.58    0.00    0.05
The data displayed as a probability matrix. 
</code></pre> <br>
'*" indicates a site poistion. '.' indicates a don't care position.

## Installation ##

Clone this repo. The bin directory contains a binary compiled under Linux Ubuntu 20.04 LTS. To compile from source cd to the src directory run make.<br>

## Usage ##

<pre><code>
$ bin/dscan
dscan $Revision: 2.3 $
   USAGE: dscan database snfile [options]
   database = sequences to be searched, in fasta format
   snfile = file with aligned segments generated using Gibbs
   options:
     -c                  - create sequence file for output
     -C                  - create motif listing
     -e< float >         - maximum E-value detected
     -E< float >         - minimum -log10(P-value) required for each motif
     -h< int >           - size of heap for saving sequences
     -M< int >           - maximum # repeats per sequence (use with -r)
     -m< int >           - minimum # repeats per sequence (use with -r)
     -n                  - use DNA alphabet and Identity matrix
     -N< int >           - print top <int> values. Ignore -e.
     -o                  - require motifs to be in correct order
     -P                  - use product multinomial model instead of Gribskov method
     -g                  - use alternative Gribskov method
     -p< float >         - pseudo counts for product multinomial model
     -r                  - scan for repeats
     -S                  - shuffle input sequences
     -s< int >           - seed for random number generator
     -X                  - mask out sequence regions NOT matching motif(s)
     -x                  - mask out sequence regions matching motif(s)
     -R                  - reverse complement search of database file
                           (a default file 'seqXXXXXXX' is created in '/tmp')
     -a                  - use PAM1 DNA matrix
     -d< int >[,< int >] - num 'sites' for each *freq* model in snfile; comma
                           separated list (100)
     -q< in t>           - maximum number of input sequences
     -v                  - don't allow overlaping motifs when using multiple models
     -B                  - alternate Bonferroni adjustment for the size of the database (single model only)
     -h                - this message
</code></pre>

