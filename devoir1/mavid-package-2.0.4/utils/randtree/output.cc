#include <stdio.h>
#include <stdlib.h>
#include "output.h"

void quit( int status )
{

  switch( status ){
  case MAVID_BAD_ARGS:
    printf( "Error: Bad arguments!\n\n" );
    break;
  case MAVID_BAD_TREE:
    printf( "Error: Either the tree file could not be found or it was not in Newick format.\n\n" );
    break;
  case MAVID_NON_BINARY:
    printf( "Error: The phylogenetic tree must be a binary tree.\n\n" );
    break;
  case MAVID_NO_SEQ_FILE:
    printf( "Error: The sequence file could not be found.\n\n" );
    break;
  case MAVID_BAD_SEQ_FILE:
    printf( "Error: The sequence file is not in multi-FASTA format.\n\n" );
    break;
  case MAVID_NO_MASK_FILE:
    printf( "Error: The masking file could not be found.\n\n" );
    break;
  case MAVID_BAD_MASK_FILE:
    printf( "Error: The masking file is not in multi-FASTA format.\n\n" );
    break;
  case MAVID_TREE_SEQ_DIFF:
    printf( "Error: The tree file does not correspond to the sequence file!\n\n" );
    break;
  }
  exit( status );

}
