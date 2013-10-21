
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "fasta.h"
#include "seq.h"
#include "file.h"


FastaStatus readFastaSeq( char *fileName, FastaSeq &seq,
			  FastaMaskStyle maskStyle )
{

  MultiFastaSeq temp;
  FastaStatus stat;

  stat = readMultiFastaSeq( fileName, temp, maskStyle );

  if( stat )
    return stat;

  if( temp.numSeqs > 1 )
    return FASTA_BAD_SEQ_FILE;

  seq = temp.seqs[0];

  delete temp.seqs;

  return FASTA_SUCCESS;

}


FastaStatus readMultiFastaSeq( char *fileName, MultiFastaSeq &seq,
			       FastaMaskStyle maskStyle )
{

  char *file = readFile( fileName );

  if( !file )
    return FASTA_NO_SEQ_FILE;

  // strip off any leading whitespace
  char *origfile = file;

  for( ; isspace(file[0]); file++ );

  // should start with a '>'
  if( file[0] != '>' ){
    delete[] origfile;
    return FASTA_BAD_SEQ_FILE;
  }

  // count the number of seqs in the file
  seq.numSeqs = 0;

  for( int i = 0; file[i]; i++ )
    if( file[i] == '>' ){
      seq.numSeqs++;
      for( ; file[i] && file[i] != '\n'; i++ );
    }

  // start reading in the seqs
  seq.seqs = new FastaSeq [seq.numSeqs];

  for( int curSeq = 0, i = 0; curSeq < seq.numSeqs; curSeq++ ){
    // advance past the '>'
    i++;
    // get to the beginning of the fasta line
    while( file[i] && isspace(file[i]) && file[i] != '\n')
      i++;

    // copy the fasta line
    seq.seqs[curSeq].fastaLine = new char [1024];
    for( int len = 0; file[i] && file[i] != '\n'; len++, i++ ){
      seq.seqs[curSeq].fastaLine[len] = file[i];
      seq.seqs[curSeq].fastaLine[len + 1] = 0;
    }

    // don't allow an empty fasta line
    if( seq.seqs[curSeq].fastaLine[0] == '\0' ){
      delete[] origfile;
      return FASTA_BAD_SEQ_FILE;
    }

    // get to the beginning of the sequence
    while( file[i] && isspace(file[i]) )
      i++;

    // find the sequence length
    int seqStart = i;
    seq.seqs[curSeq].len = 0;
    for( ; file[i] && file[i] != '>'; i++ )
      if( !isspace(file[i]) )
	seq.seqs[curSeq].len++;

    // don't allow an empty sequence
    if( seq.seqs[curSeq].len == 0 ){
      delete[] origfile;
      return FASTA_BAD_SEQ_FILE;
    }

    // copy the seq
    i = seqStart;
    seq.seqs[curSeq].seq = new char [seq.seqs[curSeq].len];
    for( int p = 0; p < seq.seqs[curSeq].len; i++ )
      if( !isspace(file[i]) )
	seq.seqs[curSeq].seq[p++] = file[i];

    // get to the beginning of the next seq
    while( file[i] && file[i] != '>' )
      i++;

  }

  // now to read in the masks
  // first let's open the right file
  if( !(maskStyle == FASTA_NO_MASKING || maskStyle == FASTA_SOFT_MASKING) ){
    char maskFileName[1024];

    delete[] origfile;
    
    strcpy( maskFileName, fileName );
    strcat( maskFileName, ".masked" );
    file = readFile( maskFileName );
    if( !file )
      return FASTA_NO_MASK_FILE;

    origfile = file;
    // strip off any leading whitespace
    for( ; isspace(file[0]); file++ );
    // should start with a '>'
    if( file[0] != '>' ){
      delete[] origfile;
      return FASTA_BAD_MASK_FILE;
    }

  }

  // now let's copy the mask
  for( int curSeq = 0, i = 0; curSeq < seq.numSeqs; curSeq++ ){

    // allocate space for the mask
    seq.seqs[curSeq].mask = new char [seq.seqs[curSeq].len];
    // skip past the fasta line
    for( ; file[i] && file[i] != '\n'; i++ );
    // copy the mask
    for( int p = 0; p < seq.seqs[curSeq].len; i++ ){
      if( !isspace(file[i]) ){
	if( !file[i] || file[i] == '>' )
	  return FASTA_BAD_MASK_FILE;
	seq.seqs[curSeq].mask[p++] = file[i];
      }
    }

    // get to the beginning of the next seq
    while( file[i] && file[i] != '>' ){
      if( !isspace( file[i] ) )
	return FASTA_BAD_MASK_FILE;
      i++;
    }
  }

  // now let's figure out what kind(i.e. DNA or PROTEIN) each sequence is
  for( int curSeq = 0; curSeq < seq.numSeqs; curSeq++ ){
    int DNACount = 0, ProtCount = 0;
    for( int i = 0; i < seq.seqs[curSeq].len; i++ ){
      if( isDNA( seq.seqs[curSeq].seq[i] ) )
	DNACount++;
      if( isProt( seq.seqs[curSeq].seq[i] ) )
	ProtCount++;
    }
    seq.seqs[curSeq].type = (DNACount >= ProtCount) ? DNA : PROTEIN;
  }

  // now let's really make the masks work and mask out bad characters
  for( int curSeq = 0; curSeq < seq.numSeqs; curSeq++ ){
    for( int i = 0; i < seq.seqs[curSeq].len; i++ ){    
      if( seq.seqs[curSeq].type == DNA ){
	if( !isDNA( seq.seqs[curSeq].seq[i] ) )
	  seq.seqs[curSeq].seq[i] = 'N';
	if( !isDNA( seq.seqs[curSeq].mask[i] ) )
	  seq.seqs[curSeq].mask[i] = 'N';
	if( maskStyle == FASTA_SOFT_MASKING && 
	    islower( seq.seqs[curSeq].mask[i] ) )
	  seq.seqs[curSeq].mask[i] = 'N';
      }
      else if( seq.seqs[curSeq].type == PROTEIN ){
	if( !isProt( seq.seqs[curSeq].seq[i] ) )
	  seq.seqs[curSeq].seq[i] = 'X';
	if( !isProt( seq.seqs[curSeq].mask[i] ) )
	  seq.seqs[curSeq].mask[i] = 'X';
	if( maskStyle == FASTA_SOFT_MASKING && 
	    islower( seq.seqs[curSeq].mask[i] ) )
	  seq.seqs[curSeq].mask[i] = 'X';
      }
    }
  }

  // convert seqs to internal format
  for( int curSeq = 0; curSeq < seq.numSeqs; curSeq++ ){
    for( int i = 0; i < seq.seqs[curSeq].len; i++ ){
      if( seq.seqs[curSeq].type == DNA ){
	seq.seqs[curSeq].seq[i] = alpha2DNA( seq.seqs[curSeq].seq[i] );
	seq.seqs[curSeq].mask[i] = alpha2DNA( seq.seqs[curSeq].mask[i] );
      }
      else if( seq.seqs[curSeq].type == PROTEIN ){
	seq.seqs[curSeq].seq[i] = alpha2Prot( seq.seqs[curSeq].seq[i] );
	seq.seqs[curSeq].mask[i] = alpha2Prot( seq.seqs[curSeq].mask[i] );
      }
    }
  }

  delete[] origfile;

  return FASTA_SUCCESS;

}


void writeFastaSeq( char *fileName, FastaSeq seq,
		    FastaMaskStyle maskStyle )
{

  MultiFastaSeq temp;

  temp.numSeqs = 1;
  temp.seqs = &seq;

  writeMultiFastaSeq( fileName, temp, maskStyle );

}


void writeMultiFastaSeq( char *fileName, MultiFastaSeq seq,
			 FastaMaskStyle maskStyle )
{

  FILE *out = fopen( fileName, "w" );

  // first output the seqs
  for( int curSeq = 0; curSeq < seq.numSeqs; curSeq++ ){
    bool empty_line = true;

    fprintf( out, ">%s\n", seq.seqs[curSeq].fastaLine );
    for( int i = 0, lC = 0; i < seq.seqs[curSeq].len; i++, lC++ ){
      char c = seq.seqs[curSeq].seq[i];

      if( seq.seqs[curSeq].type == DNA )
	c = DNA2alpha( c );
      else
	c = Prot2alpha( c );
      if( maskStyle == FASTA_SOFT_MASKING && seq.seqs[curSeq].mask[i] != seq.seqs[curSeq].seq[i] )
	c = tolower( c );
      fprintf( out, "%c", c );
      empty_line = false;
      if( lC == 50 ){
	lC = 0;
	fprintf( out, "\n" );
	empty_line = true;
      }
    }
    if( !empty_line )
      fprintf( out, "\n" );
  }
  fclose( out );

  // now output the masks
  if( maskStyle == FASTA_N_MASKING ){
    char *maskedName = new char [strlen( fileName ) + 8 ];

    strcpy( maskedName, fileName );
    strcat( maskedName, ".masked" );
    out = fopen( maskedName, "w" );

    for( int curSeq = 0; curSeq < seq.numSeqs; curSeq++ ){
      fprintf( out, ">%s\n", seq.seqs[curSeq].fastaLine );
      for( int i = 0, lC = 0; i < seq.seqs[curSeq].len; i++, lC++ ){
	char c = seq.seqs[curSeq].mask[i];

	if( seq.seqs[curSeq].type == DNA )
	  c = DNA2alpha( c );
	else
	  c = Prot2alpha( c );
	if( maskStyle == FASTA_SOFT_MASKING )
	  c = tolower( c );
	fprintf( out, "%c", c );
	if( lC == 50 ){
	  lC = 0;
	  fprintf( out, "\n" );
	}
      }
    }
    fclose( out );
  }

}


void printMultiFastaSeq( MultiFastaSeq seq, FastaMaskStyle maskStyle )
{

  // first output the seqs
  for( int curSeq = 0; curSeq < seq.numSeqs; curSeq++ ){
    bool empty_line = true;

    printf( ">%s\n", seq.seqs[curSeq].fastaLine );
    for( int i = 0, lC = 0; i < seq.seqs[curSeq].len; i++, lC++ ){
      char c = seq.seqs[curSeq].seq[i];

      if( seq.seqs[curSeq].type == DNA )
	c = DNA2alpha( c );
      else
	c = Prot2alpha( c );
      if( maskStyle == FASTA_SOFT_MASKING && seq.seqs[curSeq].mask[i] != seq.seqs[curSeq].seq[i] )
	c = tolower( c );
      printf( "%c", c );
      empty_line = false;
      if( lC == 50 ){
	lC = 0;
	printf( "\n" );
	empty_line = true;
      }
    }
    if( !empty_line )
      printf( "\n" );
  }

}


void delete_fasta_seq( FastaSeq &seq )
{

  delete[] seq.fastaLine;
  delete[] seq.seq;
  delete[] seq.mask;

}


void delete_multifasta_seq( MultiFastaSeq &seq )
{

  for( int i = 0; i < seq.numSeqs; i++ )
    delete_fasta_seq( seq.seqs[i] );
  delete[] seq.seqs;

}
