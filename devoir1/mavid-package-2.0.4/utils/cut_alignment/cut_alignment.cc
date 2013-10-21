#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "alignments.h"
#include "seq.h"

int main( int argc, char *argv[] )
{

  if( argc != 4 && argc != 5 ){
    printf( "Usage: cut_align <fasta alignment file> start end [sequence label]\n\n" );
    printf( "Returns the portion of the alignment between the given start and\n" );
    printf( "end positions. If a label is given then the coordinates are interpreted\n" );
    printf( "as being in that sequence, otherwise they are taken to be in the\n" );
    printf( "coordinates of the multiple alignment. (All coordinates begin at 1.)\n\n" );

    exit( 0 );
  }

  char *fileName, *outFileName;

  fileName = argv[1];
  outFileName = new char [strlen( fileName ) + 5];

  strcpy( outFileName, "cut." );
  strcat( outFileName, fileName );

  AlignmentsStatus stat;
  FastaAlign al;
  MultiFastaSeq seq;

  stat = readFastaAlignment( fileName, al, seq );
  if( stat != ALIGNMENTS_SUCCESS ){
    printf( "There was an error reading the file!\n\n" );
    exit( stat );
  }

  int start, end;

  start = atoi( argv[2] ) - 1;
  end = atoi( argv[3] ) - 1;

  // were we given a sequence label?
  if( argc == 5 ){
    int curSeq;

    //  associate label with sequence
    for( curSeq = 0; curSeq < seq.numSeqs; curSeq++ ){
      if( !strcmp( seq.seqs[curSeq].fastaLine, argv[4] )){
	break;
      }
    }
    if( curSeq == seq.numSeqs ){
      printf( "Error: Sequence label does not match any of the sequences.\n\n" );
      exit( 1 );
    }

    // get the starting and ending positions of the cut in alignment
    // coordinates
    int seqStart = start, seqLen = end - start;

    for( start = 0; seqStart > 0 && start < al.len; start++ ){
      if( al.line[curSeq][start] != DNA_GAP )
	seqStart--;
    }
    for( end = start + 1; seqLen > 0 && end < al.len; end++ ){
      if( al.line[curSeq][end] != DNA_GAP )
	seqLen--;
    }

    /*
    for( start = 0; al.line[curSeq][start] == DNA_GAP; start++ );
    for( int len = 0; len < seqStart; start++ )
      if( al.line[curSeq][start] != DNA_GAP )
	len++;
    end = start;
    for( int len = seqStart; len < seqEnd; end++ )
      if( al.line[curSeq][end] != DNA_GAP )
	len++;
    */

  }
  else {
    end++;
  }

  if( start < end ){
    for( int i = 0; i < al.numSeqs; i++ )
      al.line[i] += start;
    al.len = end - start;

    writeFastaAlignment( outFileName, al, seq );
    printf( "Output written to %s.\n\n", outFileName );
  }
  else {
    printf( "Error: The given range is empty.\n\n" );
  }

}
