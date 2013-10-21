

#ifndef BIOC_ALIGNMENTS
#define BIOC_ALIGNMENTS

struct ImageAlignment;
struct FastaAlign;

#include "fasta.h"
#include "tree.h"

/*
** Enums, defines, and whatnot
*/

#define GAPPED -1

enum AlignmentsStatus { ALIGNMENTS_SUCCESS, ALIGNMENTS_NO_FILE,
			ALIGNMENTS_BAD_FILE };

/*
** Structs, classes, and whatnot
*/

struct ImageAlignment {
  
  FastaSeq *firstSeq, *secSeq;
  int *firstImg, *secImg;

};

struct FastaAlign {

  int numSeqs, len;
  char **line;
  FastaSeq **seq;
  double **profile;

};

/*
** Function prototypes
*/



AlignmentsStatus readAvidAlignment( char *fileName, ImageAlignment &al );

void writeAvidAlignment( char *fileName, ImageAlignment &al );

AlignmentsStatus readBinaryAlignment( char *fileName, ImageAlignment &al );

void writeBinaryAlignment( char *fileName, ImageAlignment &al );

void writeAXTAlignment( char *fileName, ImageAlignment &al );

AlignmentsStatus readPhylipAlignment( char *fileName, FastaAlign &al,
				      MultiFastaSeq &seq );

void writePhylipAlignment( FastaAlign al, char **labels, char *fileName );

void writeClustalwAlignment( FastaAlign al, char **labels, char *fileName );

AlignmentsStatus readFastaAlignment( char *fileName, FastaAlign &al,
                                     MultiFastaSeq &seq );

AlignmentsStatus writeFastaAlignment( char *fileName, FastaAlign &al,
                                      MultiFastaSeq &seq );

AlignmentsStatus writeFastaAlignment( char *fileName, FastaAlign &al );

AlignmentsStatus printFastaAlignment( FastaAlign &al, MultiFastaSeq &seq );

FastaAlign extract_alignment( FastaAlign al, int startInd, int endInd );

FastaAlign extract_alignment( FastaAlign al, int numSeqs, int *indices, PhyloTree *tree = NULL );

FastaAlign extract_alignment( FastaAlign al, MultiFastaSeq &seqs, PhyloTree *sub_tree );

FastaAlign extract_alignment( FastaAlign &al, PhyloTree *sub_tree );

void unglue_alignment( FastaAlign al, PhyloTree *tree,
                       FastaAlign &al0, FastaAlign &al1, int *&image0, bool add_profiles );

void infer_profile( FastaAlign &al, PhyloTree *tree );

void delete_fasta_align( FastaAlign al );

double sp_score( FastaAlign al, double ****scores );

double sp_score( FastaAlign al, PhyloTree *tree );

#endif
