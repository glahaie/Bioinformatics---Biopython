#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "alignments.h"
#include "seq.h"


// THIS UTILITY ISN'T WRITTEN YET.

void usage( )
{

  printf( "Usage: translate_coords <fasta alignment file> position_1 [position_2 position_3...] [label]\n\n" );

}

int main( int argc, char *argv[] )
{

  if( argc < 3 ){
    usage();
    exit( 0 );
  }

  char *align_filename, *label;
  int *pos, num_pos;

  align_filename = argv[1];
  // is the last argument a number or a label?
  if( isdigit(argv[argc - 1][0]) ){
    label = NULL;
    num_pos = argc - 2;
  }
  else {
    label = argv[argc - 1];
    num_pos = argc - 3;
  }

  // read in the positions
  if( num_pos == 0 ){
    printf( "Error: At least one position must be specified.\n\n" );
    exit( 0 );
  }

  pos = new int [num_pos];

  if( label )
    printf( "%s\t", label );

  for( int i = 0; i < num_pos; i++ ){
    char *end_ptr;

    pos[i] = strtol( argv[i + 2], &end_ptr, 10 );
    if( *end_ptr != '\0' || pos[i] < 0 ){
      printf( "Error: \"%s\" is not a valid position.\n\n", argv[i + 2] );
      exit( 0 );
    }
    printf( "%d\t", pos[i] );
  }
  printf( "\n" );

  exit( 0 );

}
