#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "alignments.h"

int main( int argc, char *argv[] )
{

  if( argc < 2 ){
    printf( "Usage: phylip2fasta <phylip alignment file> [output file name]\n\n" );
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
    if( strcmp( fileName + strlen( fileName ) - 4, ".phy" ) == 0 )
      outFileName[strlen(fileName) - 4] = '\0';
    strcat( outFileName, ".mfa" );
  }

  // now let's do our thing.
  AlignmentsStatus stat;
  FastaAlign al;
  MultiFastaSeq seq;

  stat = readPhylipAlignment( fileName, al, seq );
  if( stat != ALIGNMENTS_SUCCESS ){
    printf( "There was an error reading the file!\n\n" );
    exit( stat );
  }

  // write the output
  writeFastaAlignment( outFileName, al, seq );

}
