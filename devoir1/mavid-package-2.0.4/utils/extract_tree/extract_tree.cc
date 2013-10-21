#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "tree.h"


int main( int argc, char *argv[] )
{

  if( argc < 4 ){
    printf( "Usage: extract_tree <tree file> label1 label2 [label3...]\n\n" );
    printf( "Extracts the subtree of the given tree with the given labels.\n\n" );
    exit( 0 );
  }

  PhyloTree *tree, *to_ret;
  
  tree = readPhyloTree( argv[1] );

  if( !tree ){
    printf( "Error reading tree!\n\n" );
    exit( 0 );
  }

  to_ret = extract_tree( tree, argv + 2, argc - 2 );

  // now make sure we found all the requested leaves.
  int num_leaves = countPhyloTreeNodes( to_ret );

  if( num_leaves == 0 ){
    printf( "Error: None of the labels were found.\n\n" );
    exit(1);
  }
  else if( num_leaves != argc - 2 ){
    char **labels = new char* [num_leaves];

    getLabels( labels, to_ret );
    printf( "Error: Some labels not found.\n" );
    printf( "The leaves which were found are: %s", labels[0] );
    for( int i = 1; i < num_leaves; i++ ){
      printf( ", %s", labels[i] );
    }
    printf( "\n\n" );
    exit(1);
  }

  printPhyloTree( to_ret );

  return 0;

}
