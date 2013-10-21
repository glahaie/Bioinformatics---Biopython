#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "tree.h"


int main( int argc, char *argv[] )
{

  PhyloTree *tree, *sub_tree;
  int num_leaves;

  if( argc < 2 || argc == 3 ){
    printf( "Usage: tree_dists <tree file> [label1 label2 [label3 ... ]]\n\n" );
    printf( "Computes pairwise distances in a tree. If no labels are given then\n" );
    printf( "all distances are computed. If only two labels are given then their\n" );
    printf( "pairwise distance is printed. If more than two labels are given then\n" );
    printf( "the pairwise distance matrix for those leaves is printed.\n\n" );
    exit( 0 );
  }

  tree = readPhyloTree( argv[1] );

  if( !tree ){
    printf( "Error reading tree!\n\n" );
    exit( 0 );
  }


  // get the distance matrix
  char **ordered_labels;
  double **dist;

  // get the appropriate subtree
  if( argc == 2 ){
    num_leaves = countPhyloTreeNodes( tree );
    dist = get_dists( tree, ordered_labels, num_leaves );
  }
  else {
    num_leaves = argc - 2;
    ordered_labels = argv + 2;
    sub_tree = extract_tree( tree, ordered_labels, num_leaves );
    if( countPhyloTreeNodes( sub_tree ) != num_leaves ){
      printf( "Error: Not all labels were found.\n\n" );
      exit(0);
    }
    dist = get_ordered_dists( sub_tree, ordered_labels, num_leaves );
  }

  // if the user requested the distance between two leaves, just print that.
  // if they didn't specifically request two leaves, print out the result in
  // matrix form even if there's only two leaves.
  if( num_leaves == 2 && argc == 4 ){
    printf( "%f\n", dist[0][1] );
  }
  else {
    // print the labels
    for( int i = 0; i < num_leaves; i++ )
      printf( "%s\n", ordered_labels[i] );
    // print the distance matrix
    for( int i = 0; i < num_leaves; i++ ){
      for( int j = 0; j < num_leaves; j++ )
	printf( "%f\t", dist[i][j] );
      printf( "\n" );
    }
  }

  return 0;

}
