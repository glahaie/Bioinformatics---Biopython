#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "alignments.h"

int main( int argc, char *argv[] )
{

  if( argc < 2 ){
    printf( "Usage: fasta2phylip <fasta alignment file> [output file name]\n\n" );
    exit( 0 );
  }

  char *fileName, *outFileName;

  // first let's prepare the output file name.
  fileName = argv[1];
  if( argc > 2 )
    outFileName = argv[2];
  else {
    outFileName = new char [strlen( fileName ) + 6];
    strcpy( outFileName, fileName);
    if( strcmp( fileName + strlen( fileName ) - 4, ".mfa" ) == 0 )
      outFileName[strlen(fileName) - 4] = '\0';
    strcat( outFileName, ".phy" );
  }

  // now let's do our thing.
  AlignmentsStatus stat;
  FastaAlign al;
  MultiFastaSeq seq;

  stat = readFastaAlignment( fileName, al, seq );
  if( stat != ALIGNMENTS_SUCCESS ){
    printf( "There was an error reading the file!\n\n" );
    exit( stat );
  }

  // make labels
  char **labels = new char* [al.numSeqs];

  for( int i = 0; i < al.numSeqs; i++ ){
    labels[i] = new char [11];
    strncpy( labels[i], seq.seqs[i].fastaLine, 10 );
    labels[i][10] = '\0';
    for( int j = 0; j < 10; j++ ){
      if( labels[i][j] == ' ' ){
	labels[i][j] = '\0';
	break;
      }
    }
  }

  // write the output
  writePhylipAlignment( al, labels, outFileName );

}
