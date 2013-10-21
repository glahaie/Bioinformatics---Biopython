#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "tree.h"
#include "output.h"
#include "fasta.h"
#include "alignments.h"
#include <time.h>


int main( int argc, char *argv[] )
{

  // check to make sure there are the right number of arguments
  if( argc != 2 ){
    printf( "Usage: randtree seq_file\n" );
    printf( "Generates a random tree on the given sequences.\n\n" );
    exit( 0 );
  }

  srand(time(NULL));

  // read in the sequences
  FastaStatus seqStatus;
  MultiFastaSeq seqs;
  
  seqStatus = readMultiFastaSeq( argv[1], seqs, FASTA_NO_MASKING );
  if( seqStatus != FASTA_SUCCESS ){
    switch( seqStatus ){
    case FASTA_NO_SEQ_FILE:
      quit( MAVID_NO_SEQ_FILE );
      break;
    case FASTA_BAD_SEQ_FILE:
      quit( MAVID_BAD_SEQ_FILE );
      break;
    default:
      break;
    }
  }

  if( seqs.numSeqs < 2 )
    quit( MAVID_BAD_SEQ_FILE );
  
  PhyloTree *res = random_tree( seqs );

  printPhyloTree( res );

  return 0;

}


