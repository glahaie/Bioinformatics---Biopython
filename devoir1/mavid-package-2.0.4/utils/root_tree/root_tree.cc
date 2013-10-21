#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "tree.h"


int main( int argc, char *argv[] )
{

  if( argc != 2 ){
    printf( "Usage: root_tree <tree file>\n\n" );
    printf( "Roots the given tree at the midpoint between the two most distant leaves.\n\n" );
    exit( 0 );
  }

  PhyloTree *tree, *a, *b, *rooted;
  float dist;

  tree = readPhyloTree( argv[1] );

  if( !tree ){
    printf( "Error reading tree!\n\n" );
    exit( 0 );
  }

  farthestPair( tree, a, b, dist );
  rooted = midpointRoot( tree, a, b, dist );
  printPhyloTree( rooted );

}
