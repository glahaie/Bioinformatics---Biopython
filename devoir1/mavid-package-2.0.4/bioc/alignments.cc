

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "alignments.h"
#include "malign.h"
#include "seq.h"
#include "fasta.h"
#include "file.h"
#include "matrices.h"


AlignmentsStatus readAvidAlignment( char *fileName, ImageAlignment &al )
{

  char *file = readFile( fileName );

  if( !file )
    return ALIGNMENTS_NO_FILE;

  al.firstSeq = new FastaSeq;
  al.secSeq = new FastaSeq;
  al.firstSeq->len = 0;
  al.secSeq->len = 0;

  // first let's get the length of each sequence
  for( int i = 0; file[i]; i++ ){
    if( isdigit( file[i] ) ){
      if( i == 0 || file[i - 1] == '\n' )
	al.firstSeq->len++;
      if( !file[i + 1] || file[i + 1] == '\n' )
	al.secSeq->len++;
    }
  }
  // allocate everything
  al.firstSeq->seq = new char [al.firstSeq->len];
  al.firstSeq->mask = new char [al.firstSeq->len];
  al.firstImg = new int [al.firstSeq->len];
  al.secSeq->seq = new char [al.secSeq->len];
  al.secSeq->mask = new char [al.secSeq->len];
  al.secImg = new int [al.secSeq->len];
  // just use empty fasta lines
  al.firstSeq->fastaLine = new char [1024];
  al.firstSeq->fastaLine[0] = 0;
  al.secSeq->fastaLine = new char [1024];
  al.secSeq->fastaLine[0] = 0;

  // now read in the alignment
  for( int i = 0, p1 = 0, p2 = 0; file[i]; ){
    // so we'll do this so at this point we're always at the start of a line
    char firstBase;
    // first base isn't gapped
    if( isdigit( file[i] ) ){
      // skip to the character for the first base
      while( isdigit( file[i] ) )
	i++;
      while( isspace( file[i] ) )
	i++;
      // if we got the the end of the line(or file) without finding a character
      // the file is not in the right format
      if( !isalpha( file[i] ) )
	return ALIGNMENTS_BAD_FILE;
      firstBase = file[i];
      i++;
      // now go to the second base or the end of the line/file
      while( file[i] && file[i] != '\n' && !isalpha(file[i]) )
	i++;
      // second base wasn't gapped
      if( isalpha( file[i] ) ){
	al.firstSeq->seq[p1] = firstBase;
	al.secSeq->seq[p2] = file[i];
	al.firstImg[p1] = p2;
	al.secImg[p2] = p1;
	p1++; p2++;
      }	
      // if it was
      else {
	al.firstSeq->seq[p1] = firstBase;
	al.firstImg[p1++] = GAPPED;
      }
    }
    // if the first base was gapped
    else {
      // now go to the second base or the end of the line/file
      while( file[i] && file[i] != '\n' && !isalpha(file[i]) )
	i++;
      // if we didn't find a character, then file is bad
      if( !isalpha( file[i] ) )
	return ALIGNMENTS_BAD_FILE;
      al.secSeq->seq[p2] = file[i];
      al.secImg[p2++] = GAPPED;
    }
    while( file[i] && file[i] != '\n' )
      i++;
    while( file[i] && isspace( file[i] ) )
      i++;
  }

  // now let's figure out what kind of sequence each one was(of course they
  // have to be the same kind).
  // TEMP - right now just make them both DNA.
  al.firstSeq->type = DNA;
  al.secSeq->type = DNA;

  // mask out any bad chars and convert from plain text.
  if( al.firstSeq->type == DNA ){
    for( int i = 0; i < al.firstSeq->len; i++ ){
      if( !isDNA( al.firstSeq->seq[i] ) )
	al.firstSeq->seq[i] = 'N';
      al.firstSeq->seq[i] = alpha2DNA( al.firstSeq->seq[i] );
    }
    for( int i = 0; i < al.secSeq->len; i++ ){
      if( !isDNA( al.secSeq->seq[i] ) )
	al.secSeq->seq[i] = 'N';
      al.secSeq->seq[i] = alpha2DNA( al.secSeq->seq[i] );
    }
  }

  strncpy( al.firstSeq->mask, al.firstSeq->seq, al.firstSeq->len );
  strncpy( al.secSeq->mask, al.secSeq->seq, al.secSeq->len );

  return ALIGNMENTS_SUCCESS;

}


void writeAvidAlignment( char *fileName, ImageAlignment &al )
{

  FILE *out = fopen( fileName, "w" );

  for( int p1 = 0, p2 = 0; p1 < al.firstSeq->len || p2 < al.secSeq->len; ){

    if( p1 > al.firstSeq->len || al.secImg[p2] < 0 ){
      fprintf( out, "     |   %c  %u\n", DNA2alpha(al.secSeq->seq[p2]),
	       p2 + 1 );
      p2++;
    }
    else if( p2 > al.secSeq->len || al.firstImg[p1] < 0 ){
      fprintf( out, "%u  %c   |\n", p1 + 1, DNA2alpha(al.firstSeq->seq[p1]) );
      p1++;
    }
    else {
      fprintf( out, "%u  %c", p1 + 1, DNA2alpha(al.firstSeq->seq[p1]) );
      if( al.firstSeq->seq[p1] == al.secSeq->seq[p2] &&
	  al.firstSeq->seq[p1] != DNA_N )
	fprintf( out, " - " );
      else
	fprintf( out, "   " );
      fprintf( out, "%c  %u\n", DNA2alpha(al.secSeq->seq[p2]), p2 + 1 );
      p1++; p2++;
    }
  }

  fclose( out );

}

#define BOTBITS 15


// TEMP - only works for DNA
AlignmentsStatus readBinaryAlignment( char *fileName, ImageAlignment &al )
{

  FILE *in = fopen( fileName, "r" );

  if( !in )
    return ALIGNMENTS_NO_FILE;

  al.firstSeq = new FastaSeq;
  al.secSeq = new FastaSeq;
  al.firstSeq->len = 0;
  al.secSeq->len = 0;

  // first let's get the length of the sequences and verify the file is valid.
  for( char c = fgetc( in ); c != EOF; c = fgetc( in ) ){
    char bot = (c & BOTBITS), top = (c >> 4);

    if( 0 > bot || bot > 5 || 0 > top || top > 5 || (top == 0 && bot == 0) ){
      delete al.firstSeq;
      delete al.secSeq;
      return ALIGNMENTS_BAD_FILE;
    }
    if( top > 0 )
      al.firstSeq->len++;
    if( bot > 0 )
      al.secSeq->len++;
  }

  // allocate everything
  al.firstSeq->seq = new char [al.firstSeq->len];
  al.firstSeq->mask = new char [al.firstSeq->len];
  al.firstImg = new int [al.firstSeq->len];
  al.secSeq->seq = new char [al.secSeq->len];
  al.secSeq->mask = new char [al.secSeq->len];
  al.secImg = new int [al.secSeq->len];
  // just use empty fasta lines
  al.firstSeq->fastaLine = new char [1024];
  al.firstSeq->fastaLine[0] = 0;
  al.secSeq->fastaLine = new char [1024];
  al.secSeq->fastaLine[0] = 0;

  // now let's go back and read in the seqs and the alignment.
  int p1 = 0, p2 = 0;
  // the 'Q' value should never be accessed.
  char bin2DNA[6] = { 'Q', DNA_A, DNA_C, DNA_T, DNA_G, DNA_N };

  fseek( in, 0, SEEK_SET );
  for( char c = fgetc( in ); c != EOF; c = fgetc( in ) ){
    char bot = (c & BOTBITS), top = (c >> 4);
    if( top > 0 ){
      al.firstSeq->seq[p1] = bin2DNA[(int)top];
      al.firstSeq->mask[p1] = bin2DNA[(int)top];
      if( bot > 0 ){
	al.secSeq->seq[p2] = bin2DNA[(int)bot];
	al.secSeq->mask[p2] = bin2DNA[(int)bot];
	al.firstImg[p1] = p2;
	al.secImg[p2] = p1;
	p1++; p2++;
      }
      else {
	al.firstImg[p1] = GAPPED;
	p1++;
      }
    }
    else {
      al.secSeq->seq[p2] = bin2DNA[(int)bot];
      al.secSeq->mask[p2] = bin2DNA[(int)bot];
      al.secImg[p2] = GAPPED;
      p2++;
    }
  }

  return ALIGNMENTS_SUCCESS;

}


// TEMP - only works for DNA
void writeBinaryAlignment( char *fileName, ImageAlignment &al )
{

  char DNA2Bin[5] = {1, 2, 4, 3, 5};
  FILE *out = fopen( fileName, "w" );
  int baseLen = al.firstSeq->len, secLen = al.secSeq->len;

  for( int p1 = 0, p2 = 0; p1 < baseLen || p2 < secLen; ){

    if( p1 > baseLen || al.secImg[p2] < 0 ){
      fputc( DNA2Bin[(int)al.secSeq->seq[p2]], out );
      p2++;
    }
    else if( p2 > secLen || al.firstImg[p1] < 0 ){
      fputc( (DNA2Bin[(int)al.firstSeq->seq[p1]] << 4), out );
      p1++;
    }
    else {
      fputc( ((DNA2Bin[(int)al.firstSeq->seq[p1]] << 4) +
	      DNA2Bin[(int)al.secSeq->seq[p2]]), out );
      p1++;
      p2++;
    }
  }

}


void writeAXTAlignment( char *fileName, ImageAlignment &al )
{

  FILE *out = fopen( fileName, "w" );

  fprintf( out, "1 base 1 %u sec 1 %u +\n", al.firstSeq->len, al.secSeq->len );

  for( int p1 = 0, p2 = 0; p1 < al.firstSeq->len || p2 < al.secSeq->len; ){
    if( p1 > al.firstSeq->len || (p2 < al.secSeq->len && al.secImg[p2] < 0) ){
      fprintf( out, "-" );
      p2++;
    }
    else {
      fprintf( out, "%c", DNA2alpha(al.firstSeq->seq[p1]) );
      p1++; p2++;
    }
  }
  fprintf( out, "\n" );

  for( int p1 = 0, p2 = 0; p1 < al.firstSeq->len || p2 < al.secSeq->len; ){
    if( p2 > al.secSeq->len || (p1 < al.firstSeq->len && al.firstImg[p1] < 0) ){
      fprintf( out, "-" );
      p1++;
    }
    else {
      fprintf( out, "%c", DNA2alpha(al.secSeq->seq[p2]) );
      p1++; p2++;
    }
  }
  fprintf( out, "\n\n" );

}


AlignmentsStatus readPhylipAlignment( char *fileName, FastaAlign &al,
                                      MultiFastaSeq &seq )
{

  FILE *in = fopen( fileName, "r" );

  if( !in )
    return ALIGNMENTS_NO_FILE;

  fscanf( in, "%u %u", &al.numSeqs, &al.len );

  for( int i = 0; i < al.numSeqs; i++ )
    al.line[i] = new char [al.len];

  

  return ALIGNMENTS_SUCCESS;

}


void writePhylipAlignment( FastaAlign al, char **labels, char *fileName )
{

  FILE *out;

  out = fopen( fileName, "w" );
  if( !out )
    return;

  // now start printing stuff out
  fprintf( out, "%6d %6d\n", al.numSeqs, al.len );

  for( int i = 0, block = 0; i < al.len; i += 50, block++ ){
    for( int seq = 0; seq < al.numSeqs; seq++ ){
      if( block == 0 )
        fprintf( out, "%-10s", labels[seq] );
      else
        fprintf( out, "          " );
      for( int k = 0; k < 50 && i + k < al.len; k++ ){
        if( (k % 10) == 0 )
          fprintf( out, " " );
        fprintf( out, "%c", DNA2alpha( al.line[seq][i + k] ) );
      }
      fprintf( out, "\n" );
    }
    fprintf( out, "\n" );
  }

}


void writeClustalwAlignment( FastaAlign al, char **labels, char *fileName )
{

  FILE *out;

  out = fopen( fileName, "w" );
  if( !out )
    return;

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
        fprintf( out, "%c", DNA2alpha( al.line[seq][i + k] ) );

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

}


AlignmentsStatus writeFastaAlignment( char *fileName, FastaAlign &al,
				      MultiFastaSeq &seq )
{

  FILE *out;

  out = fopen( fileName, "w" );

  for( int curSeq = 0; curSeq < al.numSeqs; curSeq++ ){
    fprintf( out, ">%s\n", seq.seqs[curSeq].fastaLine );
    for( int i = 0; i < al.len; i++ ){
      fprintf( out, "%c", DNA2alpha( al.line[curSeq][i] ) );
      if( (i + 1) % 60 == 0 )
	fprintf( out, "\n" );
    }
    fprintf( out, "\n" );
  }

  return ALIGNMENTS_SUCCESS;

}


AlignmentsStatus writeFastaAlignment( char *fileName, FastaAlign &al )
{

  FILE *out;

  out = fopen( fileName, "w" );

  for( int curSeq = 0; curSeq < al.numSeqs; curSeq++ ){
    fprintf( out, ">%s\n", al.seq[curSeq]->fastaLine );
    for( int i = 0; i < al.len; i++ ){
      fprintf( out, "%c", DNA2alpha( al.line[curSeq][i] ) );
      if( (i + 1) % 60 == 0 )
	fprintf( out, "\n" );
    }
    fprintf( out, "\n" );
  }

  return ALIGNMENTS_SUCCESS;

}



AlignmentsStatus printFastaAlignment( FastaAlign &al, MultiFastaSeq &seq )
{

  for( int curSeq = 0; curSeq < al.numSeqs; curSeq++ ){
    printf( ">%s\n", seq.seqs[curSeq].fastaLine );
    for( int i = 0; i < al.len; i++ ){
      printf( "%c", DNA2alpha( al.line[curSeq][i] ) );
      if( (i + 1) % 60 == 0 )
	printf( "\n" );
    }
    printf( "\n" );
  }

  return ALIGNMENTS_SUCCESS;

}



AlignmentsStatus readFastaAlignment( char *fileName, FastaAlign &al,
				     MultiFastaSeq &seq )
{

  char *file = readFile( fileName );

  if( !file ){
    al.numSeqs = 0;
    al.len = 0;
    al.line = NULL;
    al.seq = NULL;
    seq.numSeqs = 0;
    seq.seqs = NULL;
    return ALIGNMENTS_NO_FILE;
  }

  // strip off any leading whitespace
  for( ; isspace(file[0]); file++ );

  // should start with a '>'
  if( file[0] != '>' )
    return ALIGNMENTS_BAD_FILE;

  // count the number of seqs in the alignment
  al.numSeqs = 0;

  for( int i = 0; file[i]; i++ ){
    if( file[i] == '>' ){
      al.numSeqs++;
      for( ; file[i] && file[i] != '\n'; i++ );
    }
  }
  seq.numSeqs = al.numSeqs;

  // start reading in the alignment
  seq.seqs = new FastaSeq [seq.numSeqs];
  al.line = new char* [al.numSeqs];
  al.len = -1;

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

    // get to the beginning of the sequence
    while( file[i] && isspace(file[i]) )
      i++;

    // find the length of the alignment line and the sequence.
    int seqStart = i, len = 0;

    seq.seqs[curSeq].len = 0;
    for( ; file[i] && file[i] != '>'; i++ ){
      if( !isspace(file[i]) ){
        len++;
	if( file[i] != '-' )
	  seq.seqs[curSeq].len++;
      }
    }

    // all lines must have the same length
    if( al.len != -1 && al.len != len )
      return ALIGNMENTS_BAD_FILE;
    al.len = len;

    // copy the alignment line and the sequence
    i = seqStart;
    al.line[curSeq] = new char [len];
    seq.seqs[curSeq].seq = new char [seq.seqs[curSeq].len];
    seq.seqs[curSeq].mask = new char [seq.seqs[curSeq].len];
    seq.seqs[curSeq].type = DNA;

    seq.seqs[curSeq].len = 0;
    for( int p = 0; p < len; i++ ){
      if( !isspace(file[i]) ){
        al.line[curSeq][p++] = file[i];
	if( file[i] != '-' ){
	  seq.seqs[curSeq].seq[seq.seqs[curSeq].len++] = file[i];
	}
      }
    }

    // get to the beginning of the next seq
    while( file[i] && file[i] != '>' )
      i++;

  }

  // now let's go back through and convert the sequences and alignment to
  // our internal format
  for( int curSeq = 0; curSeq < al.numSeqs; curSeq++ ){
    for( int i = 0; i < al.len; i++ )
      al.line[curSeq][i] = alpha2DNA( al.line[curSeq][i] );
    for( int i = 0; i < seq.seqs[curSeq].len; i++ ){
      seq.seqs[curSeq].seq[i] = alpha2DNA( seq.seqs[curSeq].seq[i] );
      seq.seqs[curSeq].mask[i] = seq.seqs[curSeq].seq[i];
    }
  }

  // and add pointers to the sequences to the alignment
  al.seq = new FastaSeq* [al.numSeqs];

  for( int i = 0; i < al.numSeqs; i++ ){
    al.seq[i] = &seq.seqs[i];
  }

  return ALIGNMENTS_SUCCESS;
  
}


// Extracts the sequences in [startInd,endInd) from the alignment.
FastaAlign extract_alignment( FastaAlign al, int startInd, int endInd )
{

  int numSeqs = endInd - startInd;
  int *indices = new int [numSeqs];
  FastaAlign toRet;

  for( int i = startInd; i < endInd; i++ ){
    indices[i - startInd] = i;
  }

  toRet = extract_alignment( al, numSeqs, indices );

  delete[] indices;

  return toRet;

}


// Extracts the sequences indicated by indices from the alignment. If tree is not NULL then it infers the profile
// for the extracted alignment as well.
FastaAlign extract_alignment( FastaAlign al, int numSeqs, int *indices, PhyloTree *tree )
{

  FastaAlign toRet;

  if( numSeqs == 1 ){
    MultiFastaSeq temp;
    int index = indices[0];

    temp.numSeqs = 1;
    temp.seqs = al.seq[index];

    return make_trivial_alignment( temp, 0 );
  }

  toRet.numSeqs = numSeqs;
  toRet.len = 0;

  for( int i = 0; i < al.len; i++ ){
    bool allGaps = true;

    for( int j = 0; j < numSeqs; j++ ){
      allGaps &= (al.line[indices[j]][i] == DNA_GAP);
    }
    if( !allGaps ){
      toRet.len++;
    }
  }
  
  toRet.line = new char* [numSeqs];
  toRet.seq = new FastaSeq* [numSeqs];
  toRet.profile = NULL;

  for( int j = 0; j < numSeqs; j++ ){
    toRet.line[j] = new char [toRet.len + 1];
    toRet.seq[j] = al.seq[indices[j]];
  }

  toRet.len = 0;
  for( int i = 0; i < al.len; i++ ){
    bool allGaps = true;

    for( int j = 0; j < numSeqs; j++ ){
      toRet.line[j][toRet.len] = al.line[indices[j]][i];
      allGaps &= (al.line[indices[j]][i] == DNA_GAP);
    }
    if( !allGaps ){
      toRet.len++;
    }
  }

  if( tree )
    infer_profile( toRet, tree );

  return toRet;

}


FastaAlign extract_alignment( FastaAlign al, MultiFastaSeq &seqs, PhyloTree *sub_tree )
{

  int num_leaves = countPhyloTreeNodes( sub_tree );
  PhyloTree **leaves = new PhyloTree* [num_leaves];
  int *indices = new int [num_leaves];

  num_leaves = 0;
  get_leaves( sub_tree, leaves, num_leaves );

  // for each leaf, figure out which row of the alignment it corresponds to.
  for( int i = 0; i < num_leaves; i++ ){
    for( int j = 0; j < al.numSeqs; j++ ){
      if( al.seq[j] == &seqs.seqs[leaves[i]->index] ){
        indices[i] = j;
        break;
      }
    }
  }

  return extract_alignment( al, num_leaves, indices, sub_tree );

}


FastaAlign extract_alignment( FastaAlign &al, PhyloTree *sub_tree )
{

  int num_leaves = countPhyloTreeNodes( sub_tree );
  PhyloTree **leaves = new PhyloTree* [num_leaves];
  int *indices = new int [num_leaves];

  num_leaves = 0;
  get_leaves( sub_tree, leaves, num_leaves );

  // for each leaf, figure out which row of the alignment it corresponds to.
  for( int i = 0; i < num_leaves; i++ ){
    for( int j = 0; j < al.numSeqs; j++ ){
      if( al.seq[j] == leaves[i]->seq ){
        indices[i] = j;
        break;
      }
    }
  }

  return extract_alignment( al, num_leaves, indices, sub_tree );

}


void unglue_alignment( FastaAlign al, PhyloTree *tree,
		       FastaAlign &al0, FastaAlign &al1, int *&image0, bool add_profiles )
{

  // for each row of the alignment, figure out which child it belongs to
  int num_left_leaves = countPhyloTreeNodes( tree->child[0] );
  PhyloTree **left_leaves = new PhyloTree* [num_left_leaves];

  num_left_leaves = 0;
  get_leaves( tree->child[0], left_leaves, num_left_leaves );

  bool *is_left = new bool [al.numSeqs];

  for( int i = 0; i < al.numSeqs; i++ ){
    is_left[i] = false;
    for( int j = 0; j < num_left_leaves; j++ ){
      if( al.seq[i] == left_leaves[j]->seq ){
	is_left[i] = true;
        break;
      }
    }
  }

  // calculate image0. al.len is larger than we need but this shouldn't be too big of a problem.
  image0 = new int [al.len];

  for( int p = 0, cur_pos = 0, cur_image = 0; p < al.len; p++ ){
    bool cur_pos_gapped = true, image_gapped = true;

    for( int i = 0; i < al.numSeqs; i++ ){
      if( al.line[i][p] != DNA_GAP )
	if( is_left[i] )
	  cur_pos_gapped = false;
	else
	  image_gapped = false;
    }

    if( image_gapped ){
      image0[cur_pos] = -1;
      cur_pos++;
    }
    else if( cur_pos_gapped ){
      cur_image++;
    }
    else {
      image0[cur_pos] = cur_image;
      cur_pos++; cur_image++;
    }
  }

  // extract the two subalignments
  int *indices = new int [al.numSeqs];

  for( int i = 0, count = 0; i < al.numSeqs; i++ ){
    if( is_left[i] )
      indices[count++] = i;
  }
  al0 = extract_alignment( al, num_left_leaves, indices, tree->child[0] );
  
  for( int i = 0, count = 0; i < al.numSeqs; i++ ){
    if( !is_left[i] )
      indices[count++] = i;
  }
  al1 = extract_alignment( al, al.numSeqs - num_left_leaves, indices, tree->child[1] );

}


void infer_profile( FastaAlign &al, PhyloTree *tree )
{

  FastaAlign al0, al1, temp_al;
  int *image0;

  if( al.profile )
    return;

  // the unglue_alignment function will add profiles to al0 and al1
  unglue_alignment( al, tree, al0, al1, image0, true );
  // we're going to use the glue_alignments function but we need to make some fake input for it
  int *shift0, *shift1;

  shift0 = new int [al0.len];
  shift1 = new int [al1.len];

  for( int i = 0; i < al0.len; i++ )
    shift0[i] = i;
  for( int i = 0; i < al1.len; i++ )
    shift1[i] = i;

  temp_al = glue_alignments( al0, image0, shift0, al0.len, tree->child[0]->parentLen,
			     al1, NULL, shift1, al1.len, tree->child[1]->parentLen );
  al.profile = temp_al.profile;
  temp_al.profile = NULL;

  // clean up
  delete_fasta_align( temp_al );
  delete[] shift0;
  delete[] shift1;

}


void delete_fasta_align( FastaAlign al )
{

  for( int i = 0; i < al.numSeqs; i++ )
    delete[] al.line[i];
  delete[] al.line;
  delete[] al.seq;

  if( al.profile ){
    // remember the crazy way profiles were allocated.
    delete[] al.profile[0];
    delete[] al.profile;
  }

}


// yes, a quadruple pointer. frickin sweet.
double sp_score( FastaAlign al, double ****scores )
{

  double to_ret;

  to_ret = 0;

  // for each column of the alignment
  for( int p = 0; p < al.len; p++ ){
    char test;

    // for each pair of sequences without gaps in this column
    for( int i = 0; i < al.numSeqs - 1; i++ ){
      if( al.line[i][p] == DNA_GAP )
	continue;
      for( int j = i + 1; j < al.numSeqs; j++ ){
	if( al.line[j][p] == DNA_GAP )
	  continue;
	test = al.line[i][p];
	test = al.line[j][p];	
	to_ret += scores[i][j][al.line[i][p]][al.line[j][p]];
      }
    }

  }

  return to_ret;

}


// given a tree, generate the score matrices and then call the above routine.
double sp_score( FastaAlign al, PhyloTree *tree )
{

  double ****all_scores, **dists;
  double score;

  // first get the pairwise distances
  dists = get_ordered_dists( tree, al );

  // now get the pairwise score matrices
  all_scores = new double*** [al.numSeqs];
  for( int i = 0; i < al.numSeqs; i++ ){
    all_scores[i] = new double** [al.numSeqs];
    all_scores[i][i] = NULL;
    for( int j = 0; j < al.numSeqs; j++ ){
      if( i == j )
	continue;
      score_matrix( all_scores[i][j], dists[i][j] );
    }
  }

  score = sp_score( al, all_scores );

  for( int i = 0; i < al.numSeqs; i++ ){
    for( int j = i + 1; j < al.numSeqs; j++ ){
      for( int k = 0; k < 6; k++ )
	delete[] all_scores[i][j][k];
      delete[] all_scores[i][j];
    }
    delete[] all_scores[i];
    delete[] dists[i];
  }
  delete[] all_scores;
  delete[] dists;

  return score;

}
