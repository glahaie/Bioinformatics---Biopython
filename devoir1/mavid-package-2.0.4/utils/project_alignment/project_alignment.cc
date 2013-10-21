#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "alignments.h"
#include "seq.h"




int main( int argc, char *argv[] )
{

  if( argc < 4 ){
    printf( "Usage: project_alignment [-so] <fasta alignment file> <label 1> <label 2> [label 3]...\n\n" );
    exit( 0 );
  }

  AlignmentsStatus stat;
  FastaAlign al;
  MultiFastaSeq seq;
  bool to_standard_output = false;

  if( strcmp( argv[1], "-so" ) == 0 ){
    to_standard_output = true;
    argc--;
    argv++;
  }

  stat = readFastaAlignment( argv[1], al, seq );
  if( stat != ALIGNMENTS_SUCCESS ){
    printf( "There was an error reading the file!\n\n" );
    exit( stat );
  }

  // associate the sequence labels with sequences
  int num_to_extract = argc - 2;
  int *indices = new int [num_to_extract];

  for( int i = 0; i < num_to_extract; i++ ){
    int len = strlen( argv[i + 2] );
    bool was_exact, found_twice;

    indices[i] = -1;
    found_twice = was_exact = false;
    for( int j = 0; j < seq.numSeqs; j++ ){
      if( !strncmp( seq.seqs[j].fastaLine, argv[i + 2], len )){
	if( seq.seqs[j].fastaLine[len] == '\0' ){
	  if( was_exact ){
	    found_twice = true;
	    break;
	  }
	  else {
	    indices[i] = j;
	    was_exact = true;
	    found_twice = false;
	  }
	}
	else if( indices[i] != -1 && !was_exact ){
	  found_twice = true;
	}
	else {
	  indices[i] = j;
	}
      }
    }

    if( found_twice ){
      printf( "Error: The sequence label \"%s\" did not match uniquely.\n",
	      argv[i + 2] );
      exit(0);
    }
    if( indices[i] == -1 ){
      printf( "Error: The label \"%s\" does not match any sequence.\n",
	      argv[i + 2] );
      exit(0);
    }
  }

  // write the alignment
  FastaAlign extracted;
  MultiFastaSeq temp;

  temp.numSeqs = num_to_extract;
  temp.seqs = new FastaSeq [num_to_extract];

  for( int i = 0; i < num_to_extract; i++ )
    temp.seqs[i] = seq.seqs[indices[i]];

  extracted = extract_alignment( al, num_to_extract, indices );
  if( to_standard_output ){
    printFastaAlignment( extracted, temp );
  }
  else {
    int name_len = 0;
    char *out_file_name;

    for( int i = 0; i < num_to_extract; i++ )
      name_len += strlen( argv[i + 2] ) + 1;
    name_len += strlen( argv[1] ) + 1;
    out_file_name = new char [name_len];

    out_file_name[0] = '\0';
    for( int i = 0; i < num_to_extract; i++ ){
      strcat( out_file_name, argv[i + 2] );
      out_file_name[strlen( out_file_name )] = '_';
    }
    strcat( out_file_name, argv[1] );
    printf( "%s\n", out_file_name );

    writeFastaAlignment( out_file_name, extracted, temp );
  }

}
