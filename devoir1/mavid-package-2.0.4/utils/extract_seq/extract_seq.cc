#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "tree.h"


bool parse_args( int argc, char *argv[], char *&seq_file,
		 int &num_labels, char **&labels, int *&starts, int *&ends,
		 bool &inclusive )
{

  if( argc < 3 )
    return false;

  if( strcmp( argv[1], "-i" ) == 0 ){
    inclusive = true;
    argc--;
    argv++;
  }
  else
    inclusive = false;

  seq_file = argv[1];
  num_labels = argc - 2;
  labels = new char* [num_labels];
  starts = new int [num_labels];
  ends = new int [num_labels];

  for( int i = 0; i < num_labels; i++ ){
    char *string = argv[i + 2];
    int k;

    labels[i] = string;
    
    for( k = 0; string[k] != '\0' && string[k] != '['; k++ );
    if( string[k] == '\0' ){
      starts[i] = 1;
      ends[i] = -1;
    }
    else {
      if( sscanf( string + k, "[%d,%d]", &starts[i], &ends[i] ) != 2 )
	return false;
      string[k] = '\0';
    }
  }

  return true;

}


int *associate_labels( MultiFastaSeq seq, int num_labels, char **labels )
{

  int *indices = new int [num_labels];

  for( int i = 0; i < num_labels; i++ ){
    indices[i] = -1;
    for( int j = 0; j < seq.numSeqs; j++ ){
      if( !strcmp( labels[i], seq.seqs[j].fastaLine ) ){
	indices[i] = j;
	break;
      }
    }
  }

  return indices;

}


int main( int argc, char *argv[] )
{

  bool inclusive;
  char *seq_file, **labels;
  int num_labels, *indices, *starts, *ends;

  if( !parse_args( argc, argv, seq_file, num_labels, labels, starts, ends, inclusive ) ){
    printf( "Usage: extract_seq [-i] <sequence file> region1 [region2 ... ]\n\n" );
    printf( "Extracts the specified regions of the sequences in the given file. A region\n" );
    printf( "is of the form \"label\", \"label[start,end]\", or \"label[start,-1]\" (in the latter\n" );
    printf( "case, the end position is taken to be the end of the sequence). If the -i\n" );
    printf( "option is specified then any sequence which does not have a region listed will\n" );
    printf( "be included in its entirety.\n\n" );
    exit( 0 );
  }

  FastaStatus status;
  MultiFastaSeq input, output;
  bool *included, missing_labels;

  status = readMultiFastaSeq( seq_file, input, FASTA_SOFT_MASKING );
  if( status != FASTA_SUCCESS ){
    printf( "Error: There was an error reading the sequence file.\n\n" );
    exit(0);
  }
  indices = associate_labels( input, num_labels, labels );

  missing_labels = false;
  for( int i = 0; i < num_labels; i++ )
    if( indices[i] == -1 ){
      printf( "Error: Label \"%s\" was not found.\n", labels[i] );
      missing_labels = true;
    }
  if( missing_labels )
    exit(0);

  // the included array keeps track of which sequences have been included in the output.
  included = new bool [input.numSeqs];
  for( int i = 0; i < input.numSeqs; i++ )
    included[i] = false;

  // prepare the output
  output.numSeqs = num_labels;
  output.seqs = new FastaSeq [num_labels + input.numSeqs];
  for( int i = 0; i < num_labels; i++ ){
    included[indices[i]] = true;
    output.seqs[i] = input.seqs[indices[i]];
    // if the end of the sequence was requested, substitute that.
    if( ends[i] == -1 )
      ends[i] = input.seqs[indices[i]].len;

    // make sure the region actually makes sense.
    if( starts[i] > ends[i] ){
      printf( "Error: For label \"%s\", the given region is empty.\n", labels[i] );
      exit(0);
    }
    if( ends[i] > input.seqs[indices[i]].len ){
      printf( "Error: The requested ending point for label \"%s\" is after the end of the sequence.\n\n", labels[i] );
      exit(0);
    }
    if( starts[i] < 1 ){
      printf( "Error: The requested start point for label \"%s\" is before the beginning of the sequence.\n\n", labels[i] );
      exit(0);
    }

    if( starts[i] != 1 || ends[i] != input.seqs[indices[i]].len ){
      output.seqs[i].fastaLine = new char [strlen(input.seqs[indices[i]].fastaLine) + 30];
      sprintf( output.seqs[i].fastaLine, "%s from %u to %u", input.seqs[indices[i]].fastaLine, starts[i], ends[i] );
    }
    output.seqs[i].seq += starts[i] - 1;
    output.seqs[i].mask += starts[i] - 1;
    output.seqs[i].len = ends[i] - starts[i] + 1;
  }

  // if in inclusive mode, add in any sequences which were left out.
  if( inclusive ){
    for( int i = 0; i < input.numSeqs; i++ )
      if( !included[i] ){
	output.seqs[output.numSeqs] = input.seqs[i];
	output.numSeqs++;
      }
  }

  printMultiFastaSeq( output, FASTA_SOFT_MASKING );

  return 0;

}
