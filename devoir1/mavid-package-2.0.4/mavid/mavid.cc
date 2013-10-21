#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "tree.h"
#include "output.h"
#include "malign.h"
#include "fasta.h"
#include "alignments.h"
#include "matrices.h"
#include "align.h"
#include "constraints.h"
#include "refine.h"


void parse_args( int &argc, char **&argv, char *&constraints_file, char *&map_file, bool &use_refinement )
{

  if( argc == 1 )
    return;

  while( argc > 1 && argv[1][0] == '-' ){
    if( !strcmp( argv[1], "-m" ) ){
      map_file = argv[2];
      argc -= 2;
      argv += 2;
    }
    else if( !strcmp( argv[1], "-c" ) ){
      constraints_file = argv[2];
      argc -= 2;
      argv += 2;
    }
    else if( !strcmp( argv[1], "-r" ) ){
      use_refinement = true;
      argc--;
      argv++;
    }
    else
      return;

  }

}


int main( int argc, char *argv[] )
{

  PhyloTree *tree;
  int numSeqs;
  char *constraints_file = NULL, *map_file = NULL;
  Map *map = NULL;
  ConstraintsTree *cons = NULL;
  bool use_refinement = false;


  splash();

#ifdef PROFILES
  printf( "USING PROFILES\n\n" );
#endif

  // check to make sure there are the right number of arguments
  parse_args( argc, argv, constraints_file, map_file, use_refinement );

  if( argc != 3 )
    quit( MAVID_BAD_ARGS );

  // read in the tree or exit if the file isn't found
  tree = readPhyloTree( argv[1] );

  if( !tree )
    quit( MAVID_BAD_TREE );

  // make sure tree is binary
  if( !isBinaryPhyloTree( tree ) )
    quit( MAVID_NON_BINARY );

  // read in the sequences
  FastaStatus seqStatus;
  MultiFastaSeq seq;
  
  seqStatus = readMultiFastaSeq( argv[2], seq, FASTA_N_MASKING );
  if( seqStatus != FASTA_SUCCESS ){
    switch( seqStatus ){
    case FASTA_NO_SEQ_FILE:
      quit( MAVID_NO_SEQ_FILE );
      break;
    case FASTA_BAD_SEQ_FILE:
      quit( MAVID_BAD_SEQ_FILE );
      break;
    case FASTA_NO_MASK_FILE:
      quit( MAVID_NO_MASK_FILE );
      break;
    case FASTA_BAD_MASK_FILE:
      quit( MAVID_BAD_MASK_FILE );
      break;
    case FASTA_SUCCESS:
      break;
    }
  }

  // count the sequences in the file
  numSeqs = seq.numSeqs;

  if( numSeqs < 2 )
    quit( MAVID_BAD_SEQ_FILE );

  // associate tree labels with seqs
  if( !getPhyloTreeIndices( tree, seq ) )
    quit( MAVID_TREE_SEQ_DIFF );

  // make sure the tree has the same number of nodes as we have sequences
  if( numSeqs != countPhyloTreeNodes( tree ) )
    quit( MAVID_TREE_SEQ_DIFF );



  if( map_file ){
    map = read_map_file( map_file );
    if( !map ){
      printf( "Error: Map \"%s\" could not be read!\n\n", map_file );
      exit(0);
    }
  }

  if( constraints_file ){
    FILE *in = fopen( constraints_file, "r" );

    if( in ){
      cons = make_constraintsTree( in, tree );
      fclose( in );
    }

    if( !cons ){
      printf( "Error: Constraints file \"%s\" could not be read!\n\n", constraints_file );
      exit(0);
    }
  }
  else {
    cons = phylo_to_constraints( tree );
  }

  // do the alignment
  FastaAlign al;
  double **profile;

  al = multiAlign( seq, tree, cons, map );
  if( use_refinement ){
    FastaAlign refined_al;

    printf( "\nInitial alignment finished, now refining the alignment...\n" );
    refined_al = refine_leaves( al, seq, tree, profile, cons, map );
    al = refined_al;
  }



  // output the results
  // first let's get the labels for the sequences
  char **labels = new char* [al.numSeqs];

  getLabels( labels, tree );

  writePhylipAlignment( al, labels, "mavid.phy" );
  writeFastaAlignment( "mavid.mfa", al );
	  
  // done!
  printf( "MAVID worked!\n\n" );

  delete_fasta_align( al );
  delete_multifasta_seq( seq );
  for( int i = 0; i < numSeqs; i++ ){
    delete[] labels[i];
  }
  delete[] labels;

  deletePhyloTree( tree );
  delete_constraints_tree( cons );

  if( map ){
    for( int i = 0; i < map->num; i++ )
      delete[] map->segment[i];
    delete[] map;
  }

  return 0;

}
