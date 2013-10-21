#include <stdio.h>
#include <time.h>
#include "align.h"
#include "malign.h"
#include "tree.h"
#include "seq.h"
#include "constraints.h"
#include "matrices.h"

//#define MAX_PROB_GAPPED 1
#define MAX_PROB_GAPPED 0.75

FastaAlign multiAlign( MultiFastaSeq &seqs, PhyloTree *tree,
		       ConstraintsTree *cons, Map *map )
{

  FastaAlign al0, al1, to_ret;

  if( tree->numChild == 0 ){
    to_ret = make_trivial_alignment( seqs, tree->index );
    return to_ret;
  }

  al0 = multiAlign( seqs, tree->child[0], cons->child[0], map );
  al1 = multiAlign( seqs, tree->child[1], cons->child[1], map );

  //  int start_time, stop_time;

  printf( "Aligning %u versus %u\n", al0.numSeqs, al1.numSeqs );
  //  start_time = time(NULL);
  to_ret = mergeAlign( al0, al1, seqs, tree,
		      cons, map );
  //  stop_time = time(NULL);
  //  printf( "\t(Time taken: %u seconds)\n", stop_time - start_time );

  delete_fasta_align( al0 );
  delete_fasta_align( al1 );

  return to_ret;

}


FastaAlign mergeAlign( FastaAlign al0, FastaAlign al1, MultiFastaSeq &seqs, PhyloTree *tree,
		       ConstraintsTree *cons, Map *map )
{

  bp *an0, *an1;
  char *mask0, *mask1;
  int len0, len1;
  int *shift0, *shift1;
  double **aprof0 = NULL, **aprof1 = NULL;

  WeakList *a_weak = NULL;
  bool nonoverlapping, left_first;
  WeakList *weak = get_weak( map, tree, nonoverlapping, left_first );

  if( nonoverlapping ){
    return make_nonoverlapping_alignment( al0, al1, left_first );
  }


  reconstruct( al0, tree->child[0], seqs, an0, mask0, len0, shift0,
	       aprof0 );
  reconstruct( al1, tree->child[1], seqs, an1, mask1, len1, shift1,
	       aprof1 );

  // Convert the constraints to be between the ancestral sequences
  ConstraintsList *a_cons = NULL;

  ancestralize_constraints( cons, weak, al0, shift0, len0, al1, shift1, len1,
			    seqs, a_cons, a_weak );

  // generate the alignment parameters
  double **scores, gap_open, gap_extend;

  //  score_matrix( scores, 100*(tree->lens[0] + tree->lens[1]) );
  //scores = hoxd_scores( );
  scores = simple_scores();
  gap_open = scores[DNA_A][DNA_N]*5.0;
  gap_extend = gap_open/20.0;

  // align the ancestral sequences
  int *img0, *img1;
  //  int start_time, stop_time;

  //  printf( "Aligning alignments..." );
  //  start_time = time(NULL);
  globalAlign( an0, an1, mask0, mask1, len0, len1, false,
	       scores, gap_open, gap_extend, aprof0, aprof1, &img0, &img1, a_cons, a_weak );
  //  stop_time = time(NULL);
  //  printf( "(time taken: %u seconds)\n", stop_time - start_time );

  // now glue the alignments together
  FastaAlign toRet;

  toRet = glue_alignments( al0, img0, shift0, len0, tree->child[0]->parentLen,
			   al1, img1, shift1, len1, tree->child[1]->parentLen );

  delete[] an0;
  delete[] mask0;
  delete[] shift0;
  delete[] an1;
  delete[] mask1;
  delete[] shift1;
  delete[] img0;
  delete[] img1;
  delete[] aprof0;
  delete[] aprof1;

  return toRet;

}


void getIndices( PhyloTree *tree, int *indices, int &baseNum )
{

  if( tree->numChild == 0 ){
    indices[baseNum] = tree->index;
    baseNum++;
    return;
  }

  getIndices( tree->child[0], indices, baseNum );
  getIndices( tree->child[1], indices, baseNum );

}


void reconstruct( FastaAlign al, PhyloTree *tree, MultiFastaSeq &seqs,
		  bp *&seq, char *&mask, int &len, int *&shift,
		  double **&a_profile )
{

  seq = new bp [al.len];
  mask = new char [al.len];
  shift = new int [al.len];

  // the one sequence case is simple.
  if( al.numSeqs == 1 ){
    len = al.len;
    for( int i = 0; i < len; i++ ){
      seq[i] = al.line[0][i];
      mask[i] = (seqs.seqs[tree->index].mask[i] == BASE_N) ? MASKED : UNMASKED;
      shift[i] = i;
    }
    a_profile = new double* [len];
    for( int i = 0; i < len; i++ )
      a_profile[i] = al.profile[i];
    return;
  }

  // for more than one sequence, use the most likely character in the ancestral profile
  int *indices = new int [al.numSeqs];
  int junk = 0;

  getIndices( tree, indices, junk );


  // first pass, just figure out which character to select from each column
  int *pos = new int [al.numSeqs];

  for( int j = 0; j < al.numSeqs; j++ ){
    pos[j] = 0;
  }

  for( int i = 0; i < al.len; i++ ){
    double max_val = 0;
    char max_char;

    for( int k = 0; k < 4; k++ ){
      if( al.profile[i][k] > max_val ){
	max_val = al.profile[i][k];
	max_char = (char)k;
      }
    }

    mask[i] = UNMASKED;
    for( int j = 0; j < al.numSeqs; j++ ){
      if( al.line[j][i] != DNA_GAP ){
	if( pos[j] < seqs.seqs[indices[j]].len && seqs.seqs[indices[j]].mask[pos[j]] == BASE_N ){
	  mask[i] = MASKED;
	}
      }
    }

    if( al.profile[i][4] > MAX_PROB_GAPPED )
      seq[i] = DNA_GAP;
    else
      seq[i] = max_char;

  }

  // second pass, remove gaps and create shift matrix
  len = 0;
  for( int i = 0; i < al.len; i++ ){
    seq[len] = seq[i];
    mask[len] = mask[i];
    shift[len] = i;
    if( seq[i] != BASE_GAP )
      len++;
  }

  a_profile = new double* [len];
  for( int i = 0; i < len; i++ )
    a_profile[i] = al.profile[shift[i]];

  delete[] indices;
  delete[] pos;

#ifdef DEBUG
  // sanity check the results
  assert( len > 0 );
  for( int i = 0; i < len; i++ ){
    assert( 0 <= seq[i] && seq[i] <= 4 );
    assert( mask[i] == MASKED || mask[i] == UNMASKED );
    assert( shift[i] <= al.len );
    assert( i == 0 || shift[i - 1] < shift[i] );
  }
#endif

}



FastaAlign glue_alignments( FastaAlign al0, int *img0, int *shift0, int len0, float dist0,
			    FastaAlign al1, int *img1, int *shift1, int len1, float dist1 )
{

  FastaAlign toRet;

  toRet.numSeqs = al0.numSeqs + al1.numSeqs;

  // find the length of the final alignment. first count the matches in the
  // alignment of the ancestral sequences.
  int matches = 0;

  for( int i = 0; i < len0; i++ ){
    if( img0[i] >= 0 )
      matches++;
  }

  // length of the whole alignment equals length of the ancestral alignment
  // plus the number of columns left out of the ancestors.
  toRet.len = (len0 + len1 - matches) + (al0.len - len0) + (al1.len - len1);

  toRet.line = new char* [toRet.numSeqs];
  for( int i = 0; i < toRet.numSeqs; i++ )
    toRet.line[i] = new char [toRet.len];

  // add the sequences to toRet
  toRet.seq = new FastaSeq* [toRet.numSeqs];

  for( int i = 0; i < al0.numSeqs; i++ )
    toRet.seq[i] = al0.seq[i];
  for( int i = 0; i < al1.numSeqs; i++ )
    toRet.seq[al0.numSeqs + i] = al1.seq[i];

  // build a profile for the glued alignment
  double *gapped_profile = new double [5];
  double **matrix0, **matrix1;

  // prepare memory to store the profiles
  toRet.profile = new double* [toRet.len];
  toRet.profile[0] = new double [toRet.len*5];
  for( int i = 1; i < toRet.len; i++ )
    toRet.profile[i] = toRet.profile[i - 1] + 5;

  // prepare stuff to do profile inference
  gapped_profile[0] = gapped_profile[1] = gapped_profile[2] = gapped_profile[3] = 0.0;
  gapped_profile[4] = 1.0;
  transition_matrix_gapped( matrix0, 100.0*dist0 );
  transition_matrix_gapped( matrix1, 100.0*dist1 );

  // now convert the alignment of the ancestral sequences to an alignment
  // of the alignments. we only need the image of al0 in al1.
  int *tempImg0 = new int [al0.len];

  for( int i = 0; i < al0.len; i++ )
    tempImg0[i] = -1;

  for( int i = 0; i < len0; i++ ){
    if( img0[i] >= 0 )
      tempImg0[shift0[i]] = shift1[img0[i]];
  }

  toRet.len = 0;
  int lastImg = -1;
  for( int i = 0; i < al0.len; i++ ){
    if( tempImg0[i] == -1 ){
      for( int j = 0; j < al0.numSeqs; j++ )
	toRet.line[j][toRet.len] = al0.line[j][i];
      for( int j = al0.numSeqs; j < toRet.numSeqs; j++ )
	toRet.line[j][toRet.len] = BASE_GAP;
      infer_profile( al0.profile[i], matrix0, gapped_profile, matrix1, toRet.profile[toRet.len] );
      toRet.len++;
    }
    else {
      for( int pos = lastImg + 1; pos < tempImg0[i]; pos++ ){
	for( int j = 0; j < al0.numSeqs; j++ )
	  toRet.line[j][toRet.len] = BASE_GAP;
	for( int j = al0.numSeqs; j < toRet.numSeqs; j++ )
	  toRet.line[j][toRet.len] = al1.line[j - al0.numSeqs][pos];
	infer_profile( gapped_profile, matrix0, al1.profile[pos], matrix1, toRet.profile[toRet.len] );
	toRet.len++;
      }
      lastImg = tempImg0[i];
      for( int j = 0; j < al0.numSeqs; j++ )
	toRet.line[j][toRet.len] = al0.line[j][i];
      for( int j = al0.numSeqs; j < toRet.numSeqs; j++ )
	toRet.line[j][toRet.len] = al1.line[j - al0.numSeqs][tempImg0[i]];
      infer_profile( al0.profile[i], matrix0, al1.profile[tempImg0[i]], matrix1, toRet.profile[toRet.len] );
      toRet.len++;
    }      
  }
  for( int i = lastImg + 1; i < al1.len; i++ ){
    for( int j = 0; j < al0.numSeqs; j++ )
      toRet.line[j][toRet.len] = BASE_GAP;
    for( int j = al0.numSeqs; j < toRet.numSeqs; j++ )
      toRet.line[j][toRet.len] = al1.line[j - al0.numSeqs][i];
    infer_profile( gapped_profile, matrix0, al1.profile[i], matrix1, toRet.profile[toRet.len] );
    toRet.len++;
  }

  for( int i = 0; i < 5; i++ ){
    delete[] matrix0[i];
    delete[] matrix1[i];
  }
  delete[] matrix0;
  delete[] matrix1;
  delete[] gapped_profile;
  delete[] tempImg0;

#ifdef DEBUG
  // extract the two old alignments from the new one and check that they match.
  FastaAlign old;

  // check first alignment
  old = extract_alignment( toRet, 0, al0.numSeqs );

  assert( old.numSeqs == al0.numSeqs && old.len == al0.len );
  for( int i = 0; i < al0.len; i++ )
    for( int j = 0; j < al0.numSeqs; j++ )
      assert( old.line[j][i] == al0.line[j][i] );

  for( int i = 0; i < old.numSeqs; i++ )
    delete[] old.line[i];
  delete[] old.line;
  delete[] old.seq;

  // check second alignment
  old = extract_alignment( toRet, al0.numSeqs, toRet.numSeqs );

  assert( old.numSeqs == al1.numSeqs && old.len == al1.len );
  for( int i = 0; i < al1.len; i++ )
    for( int j = 0; j < al1.numSeqs; j++ )
      assert( old.line[j][i] == al1.line[j][i] );

  for( int i = 0; i < old.numSeqs; i++ )
    delete[] old.line[i];
  delete[] old.line;
  delete[] old.seq;
#endif

  return toRet;

}



FastaAlign make_trivial_alignment( MultiFastaSeq &seqs, int index )
{

  FastaAlign to_ret;

  to_ret.numSeqs = 1;
  to_ret.line = new char* [1];
  to_ret.len = seqs.seqs[index].len;
  to_ret.line[0] = new char [to_ret.len];
  for( int i = 0; i < to_ret.len; i++ )
    to_ret.line[0][i] = seqs.seqs[index].seq[i];
  to_ret.seq = new FastaSeq* [1];
  to_ret.seq[0] = &seqs.seqs[index];
  to_ret.profile = new double* [to_ret.len];
  to_ret.profile[0] = new double [to_ret.len*5];
  for( int i = 0; i < to_ret.len; i++ ){
    if( i != 0 )
      to_ret.profile[i] = to_ret.profile[i - 1] + 5;
    for( int j = 0; j < 5; j++ )
      to_ret.profile[i][j] = 0.0;
    if( seqs.seqs[index].seq[i] == BASE_N ){
      to_ret.profile[i][BASE_A] = to_ret.profile[i][BASE_C] = to_ret.profile[i][BASE_G] = to_ret.profile[i][BASE_T] = 
	0.25;
    }
    else {
      to_ret.profile[i][seqs.seqs[index].seq[i]] = 1.0;
    }
  }

  return to_ret;

}



FastaAlign make_nonoverlapping_alignment( FastaAlign al0, FastaAlign al1, bool left_first )
{

  FastaAlign to_ret;

  to_ret.numSeqs = al0.numSeqs + al1.numSeqs;
  to_ret.len = al0.len + al1.len;
  to_ret.line = new char* [to_ret.numSeqs];
  for( int i = 0; i < to_ret.numSeqs; i++ )
    to_ret.line[i] = new char [to_ret.len];

  int left_offset, right_offset;

  if( left_first ){
    left_offset = 0;
    right_offset = al0.len;
  }
  else {
    left_offset = al1.len;
    right_offset = 0;
  }

  for( int p = left_offset; p < left_offset + al0.len; p++ ){
    for( int i = 0; i < al0.numSeqs; i++ )
      to_ret.line[i][p] = al0.line[i][p - left_offset];
    for( int i = 0; i < al1.numSeqs; i++ )
      to_ret.line[al0.numSeqs + i][p] = DNA_GAP;
  }
  for( int p = right_offset; p < right_offset + al1.len; p++ ){
    for( int i = 0; i < al0.numSeqs; i++ )
      to_ret.line[i][p] = DNA_GAP;
    for( int i = 0; i < al1.numSeqs; i++ )
      to_ret.line[al0.numSeqs + i][p] = al1.line[i][p - right_offset];
  }

  return to_ret;

}



void infer_profile( double *prof0, double **matrix0, double *prof1, double **matrix1, double *dest )
{

  /*  for( int i = 0; i < 5; i++ ){
    dest[i] = (prof0[i] + prof1[i])/2.0;
    }*/


  for( int i = 0; i < 5; i++ ){
    dest[i] = 0;
  }

  for( int j = 0; j < 5; j++ ){
    for( int k = 0; k < 5; k++ ){
      double norm_factor = 0;

      for( int i = 0; i < 5; i++ )
	norm_factor += matrix0[i][j]*matrix1[i][k];

      for( int i = 0; i < 5; i++ )
	dest[i] += matrix0[i][j]*matrix1[i][k]*prof0[j]*prof1[k]/norm_factor;
    }
  }


}
