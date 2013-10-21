#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "fasta.h"

int main( int argc, char *argv[] )
{

  if( argc != 2 ){
    printf( "Usage: checkfasta <file>\n\n" );
    printf( "If the given file is a valid fasta sequence file, checkfasta prints (and\n" );
    printf( "returns) the number of sequences in it. If not (or if the file is not found),\n" );
    printf( "then it prints (and returns) zero.\n" );

    exit( 0 );
  }

  // now let's do our thing.
  MultiFastaSeq seq;
  FastaStatus stat;

  stat = readMultiFastaSeq( argv[1], seq, FASTA_NO_MASKING );
  if( stat != FASTA_SUCCESS ){
    printf( "0\n" );
    return 0;
  }

  printf( "%u\n", seq.numSeqs );
  return seq.numSeqs;

}
