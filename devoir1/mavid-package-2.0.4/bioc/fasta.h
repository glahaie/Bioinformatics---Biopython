
#ifndef BIOC_FASTA
#define BIOC_FASTA

/*
** Enums, defines, and whatnot
*/

enum FastaStatus { FASTA_SUCCESS, FASTA_NO_SEQ_FILE, FASTA_BAD_SEQ_FILE,
		   FASTA_NO_MASK_FILE, FASTA_BAD_MASK_FILE };

enum FastaMaskStyle { FASTA_NO_MASKING, FASTA_N_MASKING, FASTA_SOFT_MASKING };

enum SeqType { DNA, PROTEIN };


/*
** Structs, classes, and whatnot
*/

struct FastaSeq {

  char *fastaLine;
  char *seq;
  char *mask;
  int len;
  SeqType type;

};

struct MultiFastaSeq {

  int numSeqs;
  FastaSeq *seqs;

};


/*
** Function prototypes
*/

FastaStatus readFastaSeq( char *fileName, FastaSeq &seq,
			  FastaMaskStyle maskStyle );

void writeFastaSeq( char *fileName, FastaSeq seq,
		    FastaMaskStyle maskStyle );

FastaStatus readMultiFastaSeq( char *fileName, MultiFastaSeq &seq,
			       FastaMaskStyle maskStyle );

void writeMultiFastaSeq( char *fileName, MultiFastaSeq seq,
			 FastaMaskStyle maskStyle );

void printMultiFastaSeq( MultiFastaSeq seq, FastaMaskStyle maskStyle );

void delete_fasta_seq( FastaSeq &seq );
void delete_multifasta_seq( MultiFastaSeq &seq );

#endif
