#include <stdio.h>
#include "common.h"
#include "output.h"
#include "malign.h"
#include "fasta.h"

#define VERSION "(version 2.0, build 4)"

void splash()
{

  cout << endl
       << "*******************************************************" << endl
       << "*                                                     *" << endl
       << "*                  Welcome to MAVID.                  *" << endl
       << "*                  " << VERSION;

  for( int i = 35 - strlen(VERSION); i > 0; i-- )
    cout << " ";

  cout << "*" << endl
       << "*                                                     *" << endl
       << "*******************************************************" << endl
       << endl << endl;
}


void usage( )
{

  cout << "Usage:     mavid [args] tree_file seq_file     " << endl << endl
       << "    tree_file: phylogenetic tree file in Newick format." << endl
       << "    seq_file: file containing the sequences to be aligned in "
       << "multiFasta format." << endl << endl
       << "Arguments to mavid:" << endl << endl

       << "-h : Print this message." << endl 
       << "-r : Perform iterative refinement at the leaves." << endl 
       << "-c <file>: Read constraints from the given file." << endl
       << "-m <file>: Read map from the given file." << endl << endl


       << endl << endl;

}

void quit( int status )
{

  switch( status ){
  case MAVID_BAD_ARGS:
    usage( );
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


bool outputClustalw( FastaAlign al, PhyloTree *tree, char *fileName )
{

  FILE *out;

  out = fopen( fileName, "w" );
  if( !out )
    return false;

  // first let's get the labels for the sequences
  char **labels = new char* [al.numSeqs];

  getLabels( labels, tree );

  // get the max label length and then the line length
  unsigned int maxll = 0, linelen;
  for( int i = 0; i < al.numSeqs; i++ )
    if( strlen( labels[i] ) > maxll )
      maxll = strlen( labels[i] );

  // make room for the labels
  linelen = 80 - maxll - 1;
  // then round to a multiple of 10
  linelen -= linelen % 10;
  // now adjust maxll to use up the rest of the space
  maxll = 80 - linelen;

  for( int i = 0, block = 0; i < al.len; i += 60, block++ ){
    for( int seq = 0; seq < al.numSeqs; seq++ ){
      fprintf( out, "%s", labels[seq] );

      for( int k = 0; k < 60; k++ )
	fprintf( out, "%c", ALPHA[(int)al.line[seq][i + k]] );

      fprintf( out, "\n" );
    }
    for( int k = 0; k < 60; k++ ){
      bool allSame = true;
      for( int seq = 0; seq < al.numSeqs; seq++ ){
	if( seq > 0 && al.line[seq - 1][i + k] != al.line[seq][i + k] ){
	  allSame = false;
	  break;
	}
      }
      if( allSame )
	fprintf( out, "*" );
      else
	fprintf( out, " " );
    }

    fprintf( out, "\n" );
  }

  return true;

}
