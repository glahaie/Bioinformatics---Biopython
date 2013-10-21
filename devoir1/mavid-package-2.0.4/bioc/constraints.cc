#include <string.h>
#include <stdio.h>
#include "tree.h"
#include "constraints.h"
#include "match.h"


void delete_constraints_tree( ConstraintsTree *tree )
{

  for( int i = 0; i < tree->numChild; i++ )
    delete_constraints_tree( tree->child[i] );
  if( tree->numChild > 0 ){
    delete[] tree->child;
  }
  delete tree;

}


void get_leaves( ConstraintsTree *tree, ConstraintsTree **leaves, int &num )
{

  if( !tree )
    return;

  if( tree->numChild == 0 ){
    leaves[num] = tree;
    num++;
  }

  for( int i = 0; i < tree->numChild; i++ )
    get_leaves( tree->child[i], leaves, num );

}


Constraint read_constraint( FILE *in, ConstraintsTree **leaves, int num_seqs )
{

  //    char first_label[1024], second_label[1024];
  char *first_label = new char [1024], *second_label = new char [1024];
  Constraint res;

  fscanf( in, "%s %u %u %s %u %u\n", first_label, &res.first_start,
	  &res.first_end, second_label, &res.second_start,
          &res.second_end );

  for( int i = 0; i < num_seqs; i++ ){
    if( strcmp( first_label, leaves[i]->label ) == 0 )
      res.first_index = leaves[i]->index;
    if( strcmp( second_label, leaves[i]->label ) == 0 )
      res.second_index = leaves[i]->index;
  }

  return res;

}


ConstraintsTree *phylo_to_constraints( PhyloTree *tree, ConstraintsTree *parent )
{

  if( !tree )
    return NULL;

  ConstraintsTree *res = new ConstraintsTree;

  if( tree->label ){
    res->label = new char [strlen(tree->label) + 1];
    strcpy( res->label, tree->label );
  }
  res->index = tree->index;
  res->parent = parent;
  res->list = NULL;
  res->numChild = tree->numChild;
  res->child = new ConstraintsTree* [res->numChild];
  for( int i = 0; i < tree->numChild; i++ ){
    res->child[i] = phylo_to_constraints( tree->child[i], res );
  }

  return res;

}


ConstraintsTree *get_lca( ConstraintsTree *leaf1, int depth1,
			  ConstraintsTree *leaf2, int depth2 )
{

  while( leaf1 != leaf2 ){
    if( depth1 > depth2 ){
      leaf1 = leaf1->parent;
      depth1--;
    }
    else {
      leaf2 = leaf2->parent;
      depth2--;
    }
  }

  return leaf1;

}


ConstraintsTree ***make_lca_matrix( ConstraintsTree **leaves, int num_seqs )
{

  int *depths = new int [num_seqs];
  ConstraintsTree ***res = new ConstraintsTree** [num_seqs];

  for( int i = 0; i < num_seqs; i++ )
    res[i] = new ConstraintsTree* [num_seqs];

  for( int i = 0; i < num_seqs; i++ ){
    depths[i] = 0;
    for( ConstraintsTree *temp = leaves[i]; temp->parent != NULL;
	 temp = temp->parent )
      depths[i]++;
  }

  for( int i = 0; i < num_seqs; i++ ){
    res[leaves[i]->index][leaves[i]->index] = leaves[i];
    for( int j = i + 1; j < num_seqs; j++ )
      res[leaves[i]->index][leaves[j]->index] = 
	res[leaves[j]->index][leaves[i]->index] = 
	get_lca( leaves[i], depths[i], leaves[j], depths[j] );
  }

  return res;

}


void add_constraint( ConstraintsList *&list, Constraint cons )
{

  ConstraintsList *to_add = new ConstraintsList;

  to_add->cons = cons;
  if( list )
    to_add->next = list;
  else
    to_add->next = NULL;
  
  list = to_add;

}

ConstraintsTree *make_constraintsTree( FILE *in, PhyloTree *tree )
{

  // make empty tree
  ConstraintsTree *res;

  res = phylo_to_constraints( tree );

  // make index pair to tree pointer translator (i.e. lca matrix)
  int num_seqs = countPhyloTreeNodes( tree ), junk = 0;
  ConstraintsTree **leaves = new ConstraintsTree* [num_seqs];
  ConstraintsTree ***lca_matrix;

  get_leaves( res, leaves, junk );
  lca_matrix = make_lca_matrix( leaves, num_seqs );

  // parse the contraints file
  int thig = fgetc( in );

  if( !feof(in) )
    ungetc( thig, in );
  while( !feof(in) ){
    Constraint cons;
    ConstraintsTree *lca;

    cons = read_constraint( in, leaves, num_seqs );
    lca = lca_matrix[cons.first_index][cons.second_index];
    add_constraint( lca->list, cons );
  }

  return res;

}


int **make_ancestral_map( FastaAlign al, int *shift, int len,
			  ConstraintsTree **leaves, int *lens )
{

  // translation from multi-alignment coord to ancestral coord
  int *unshift = new int [al.len];

  for( int i = 0, last = -1; i < len; i++ ){
    for( int p = last + 1; p <= shift[i]; p++ )
      unshift[p] = i;
    last = shift[i];
  }
  for( int p = shift[len - 1]; p < al.len; p++ )
    unshift[p] = len - 1;

#ifdef DEBUG
  for( int i = 0; i < len; i++ )
    assert( unshift[shift[i]] == i );
#endif

  // shifts will translate seq -> ancestral
  int **shifts = new int* [al.numSeqs], *cur_pos = new int [al.numSeqs];  

  for( int i = 0; i < al.numSeqs; i++ ){
    shifts[i] = new int [lens[leaves[i]->index]];
    cur_pos[i] = 0;
  }

  for( int p = 0; p < al.len; p++ ){
    for( int i = 0; i < al.numSeqs; i++ ){
      if( al.line[i][p] != DNA_GAP ){
	shifts[i][cur_pos[i]++] = unshift[p];
      }
    }
  }

  delete[] unshift;

  return shifts;

}



int cons_compare( const void *a, const void *b )
{

  ConstraintsList *c1 = *((ConstraintsList **)a),
    *c2 = *((ConstraintsList **)b);

  if( c1->cons.first_start > c2->cons.first_start )
    return 1;
  else if( c1->cons.first_start < c2->cons.first_start )
    return -1;
  else {
    if( c1->cons.second_start > c2->cons.second_start )
      return 1;
    else
      return -1;
  }
}


void choose_consistent( ConstraintsList *&list )
{

  if( !list )
    return;

  int num_anchors = 0;

  for( ConstraintsList *temp = list; temp != NULL; temp = temp->next )
    num_anchors++;

  ConstraintsList **lists = new ConstraintsList* [num_anchors];

  num_anchors = 0;
  for( ConstraintsList *temp = list; temp != NULL; temp = temp->next ){
    lists[num_anchors] = temp;
    num_anchors++;
  }

  qsort( lists, num_anchors, sizeof( ConstraintsList* ), cons_compare );
  list = lists[0];

  ConstraintsList *tail = list;

  for( int i = 1; i < num_anchors; i++ ){
    if( lists[i]->cons.first_start >= tail->cons.first_end &&
	lists[i]->cons.second_start >= tail->cons.second_end ){
      tail->next = lists[i];
      tail = lists[i];
    }
    else
      delete lists[i];
  }
  tail->next = NULL;

}


void ancestralize_constraints( ConstraintsTree *cons, WeakList *weak,
			       FastaAlign al0, int *shift0, int len0,
			       FastaAlign al1, int *shift1, int len1,
			       MultiFastaSeq &seqs,
			       ConstraintsList *&a_cons, WeakList *&a_weak )
{

  a_cons = NULL;
  a_weak = NULL;

  if( !cons->list && !weak )
    return;

  // get leaves from the constraint tree
  ConstraintsTree **leaves =
    new ConstraintsTree* [al0.numSeqs + al1.numSeqs];
  int num_leaves = 0, max_idx = 0;

  get_leaves( cons, leaves, num_leaves );
  for( int i = 0; i < num_leaves; i++ )
    if( leaves[i]->index > max_idx )
      max_idx = leaves[i]->index;

  // to reduce code complexity, make an array of sequence lengths
  int *lens = new int [seqs.numSeqs];

  for( int i = 0; i < seqs.numSeqs; i++ ){
    lens[i] = seqs.seqs[i].len;
  }

  // get the maps from the sequences to their ancestral sequences.
  int **shifts0, **shifts1;

  shifts0 = make_ancestral_map( al0, shift0, len0, leaves, lens );
  shifts1 = make_ancestral_map( al1, shift1, len1, leaves + al0.numSeqs,
				lens );

  int **shifts = new int* [max_idx + 1], *which = new int [max_idx + 1];

  for( int i = 0; i < al0.numSeqs; i++ ){
    shifts[leaves[i]->index] = shifts0[i];
    which[leaves[i]->index] = 0;
  }
  for( int i = 0; i < al1.numSeqs; i++ ){
    shifts[leaves[i + al0.numSeqs]->index] = shifts1[i];
    which[leaves[i + al0.numSeqs]->index] = 1;
  }

  // now ancestralize the constraints
  Constraint to_add;

  for( ConstraintsList *temp = cons->list; temp != NULL; temp = temp->next ){
    to_add.first_start =
      shifts[temp->cons.first_index][temp->cons.first_start];
    to_add.first_end = shifts[temp->cons.first_index][temp->cons.first_end - 1] + 1;
    to_add.second_start =
      shifts[temp->cons.second_index][temp->cons.second_start];
    to_add.second_end = shifts[temp->cons.second_index][temp->cons.second_end - 1] + 1;
#ifdef DEBUG
    assert( which[temp->cons.first_index] != which[temp->cons.second_index] );
#endif
    to_add.first_index = 0;
    to_add.second_index = 1;
    if( which[temp->cons.first_index] == 1 ){
      swap( to_add.first_start, to_add.second_start );
      swap( to_add.first_end, to_add.second_end );
      swap( to_add.first_index, to_add.second_index );
    }

    // DELETE ME
    //    printf( "%u %u %u %u %u %u -> %u %u an0 %u %u an1\n", temp->cons.first_start, temp->cons.first_end, temp->cons.first_index, temp->cons.second_start, temp->cons.second_end, temp->cons.second_index, to_add.first_start, to_add.first_end, to_add.second_start, to_add.second_end );

    add_constraint( a_cons, to_add );

  }

  unionize( a_cons );

  for( ConstraintsList *temp = a_cons; temp != NULL; temp = temp->next )
    printf( "%u %u %u %u\n", temp->cons.first_start, temp->cons.first_end,
	    temp->cons.second_start, temp->cons.second_end );

  // now for the other one
  WeakList *tail = NULL;

  for( WeakList *temp = weak; temp != NULL; temp = temp->next ){
    WeakList *to_add = new WeakList;

    to_add->weak.i = shifts[temp->weak.i_index][temp->weak.i];
    to_add->weak.j= shifts[temp->weak.j_index][temp->weak.j];
    to_add->weak.i_index = 0;
    to_add->weak.j_index = 0;
    to_add->weak.type = temp->weak.type;
    to_add->next = NULL;

    if( !a_weak )
      a_weak = tail = to_add;
    else {
      tail->next = to_add;
      tail = to_add;
    }
  }

  for( int i = 0; i < num_leaves; i++ ){
    delete[] shifts[leaves[i]->index];
  }
  delete[] shifts0;
  delete[] shifts1;
  delete[] shifts;
  delete[] which;

}


Map *read_map_file( char *file_name )
{

  FILE *in = fopen( file_name, "r" );

  if( !in ){
    return NULL;
  }

  // count the number of segments and the width
  Map *to_ret = new Map;

  to_ret->width = 1;
  for( char c = fgetc( in ); c != '\n'; c = fgetc( in ) ){
    if( c == '\t' )
      to_ret->width++;
  }

  to_ret->num = 1; // we already read in the first line above   
  while( !feof(in) ){
    if( fgetc( in ) == '\n' )
      to_ret->num++;
  }

  // rewind to the beginning, allocate space and read in the segments
  fseek( in, 0, SEEK_SET );

  to_ret->segment = new Range* [to_ret->num];
  for( int i = 0; i < to_ret->num; i++ ){
    to_ret->segment[i] = new Range [to_ret->width];

    fscanf( in, "%d %d", &to_ret->segment[i][0].start,
	    &to_ret->segment[i][0].end );
    for( int j = 1; j < to_ret->width; j++ ){
      fscanf( in, "\t%d %d", &to_ret->segment[i][j].start,
	      &to_ret->segment[i][j].end );
    }
    fscanf( in, "\n" );
  }

  fclose( in );

  return to_ret;

}


WeakList *get_weak( Map *map, PhyloTree *tree, bool &nonoverlapping,
		    bool &left_first )
{

  if( !map ){
    nonoverlapping = false;
    return NULL;
  }

  // get the children of this tree
  int left_num, right_num, junk = 0;
  PhyloTree **left, **right;

  left_num = countPhyloTreeNodes( tree->child[0] );
  right_num = countPhyloTreeNodes( tree->child[1] );
  left = new PhyloTree* [left_num];
  right = new PhyloTree* [right_num];
  get_leaves( tree->child[0], left, junk );
  junk = 0;
  get_leaves( tree->child[1], right, junk );

  // run through each line in the map
  Range **range = new Range* [map->num];
  bool *is_left = new bool [map->num];
  int *index = new int [map->num];
  int *map_pos = new int [map->num];
  bool found_anchor = false;
  int num = 0;

  for( int i = 0; i < map->num; i++ ){
    bool found_left, found_right;

    found_left = found_right = false;
    // run through all the left children and see if they're in this segment
    for( int k = 0; k < left_num; k++ ){
      if( map->segment[i][left[k]->index].start != -1 ){
	range[num] = &map->segment[i][left[k]->index];
	is_left[num] = true;
	index[num] = left[k]->index;
	map_pos[num] = i;

	found_left = true;
	break;
      }
    }
	
    // run through all the right children and see if they're in this segment
    for( int k = 0; k < right_num; k++ ){
      if( map->segment[i][right[k]->index].start != -1 ){
	range[num] = &map->segment[i][right[k]->index];
	is_left[num] = false;
	index[num] = right[k]->index;
	map_pos[num] = i;

	found_right = true;
	break;
      }
    }

    if( found_left && found_right ){
      found_anchor = true;
      continue;
    }
    if( !found_left && !found_right )
      continue;

    num++;
  }

  // check if they're non-overlapping
  // !!!!!!!!!! we really should check whether they have an anchor that's
  // not from the map !!!!!!!!!!!

  nonoverlapping = false;
  if( !found_anchor ){
    //    bool found_transition = false, overlapping = false;
    int num_transitions = 0, last_transition = -1;

    for( int i = 1; i < num; i++ ){
      if( is_left[i - 1] != is_left[i] ){
	num_transitions++;
	last_transition = i;
      }
    }

    if( num_transitions == 1 && (map_pos[last_transition] > map_pos[last_transition - 1] + MAX_MISSING) ){
      left_first = is_left[last_transition - 1];
      nonoverlapping = true;
      return NULL;
    }
  }

  // turn the ranges into constraints
  WeakList *to_ret = NULL;

  for( int i = num - 1; i > 0; i-- ){
    WeakList *to_add = new WeakList;

    if( is_left[i] == is_left[i - 1] )
      continue;
    if( is_left[i - 1] ){
      to_add->weak.i = range[i - 1]->end;
      to_add->weak.j = range[i]->start;
      to_add->weak.type = LEQ;
      to_add->weak.i_index = index[i - 1];
      to_add->weak.j_index = index[i];
    }
    else {
      to_add->weak.i = range[i]->start;
      to_add->weak.j = range[i - 1]->end;
      to_add->weak.type = GEQ;
      to_add->weak.i_index = index[i];
      to_add->weak.j_index = index[i - 1];
    }
    to_add->next = to_ret;
    to_ret = to_add;
  }

  return to_ret;

}



void filter_by_weak( matchlist *&head, WeakList *weaks )
{

  if( !weaks )
    return;

  // the N constraints divide each sequence into N + 1 regions
  int num_weaks = 0;

  for( WeakList *temp = weaks; temp != NULL; temp = temp->next )
    num_weaks++;

  int *first_pos = new int [num_weaks + 2],
    *second_pos = new int [num_weaks + 2];
  WeakType *type = new WeakType [num_weaks];

  first_pos[0] = second_pos[0] = -1;
  first_pos++; second_pos++;
  num_weaks = 0;
  for( WeakList *temp = weaks; temp != NULL; temp = temp->next ){
    first_pos[num_weaks] = temp->weak.i;
    second_pos[num_weaks] = temp->weak.j;
    type[num_weaks] = temp->weak.type;
    num_weaks++;
  }
  first_pos[num_weaks] = second_pos[num_weaks] = 1000000000;
  

  // run through all the matches
  for( matchlist *cur = head, *last = NULL; cur != NULL; ){    
    int bin;
    bool is_rejected;

    // bin it in the first sequence
    for( bin = 0; first_pos[bin] < cur->basePos; bin++ );

    // if it violates the constraints throw it out
    is_rejected = false;

    if( cur->secPos < second_pos[bin - 1] ){
      if( cur->secPos < second_pos[bin - 2] )
	is_rejected = true;
      else
	is_rejected = (type[bin - 1] != LEQ);
    }
    else if( cur->secPos > second_pos[bin] ){
      if( cur->secPos > second_pos[bin + 1] )
	is_rejected = true;
      else
	is_rejected = (type[bin] != GEQ);
    }

    if( is_rejected ){
      if( last ){
	// in the middle of the list. splice this guy out.
	last->next = cur->next;
	delete cur;
	cur = last->next;
      }
      else {
	// we must be at the head of the list
	head = cur->next;
	delete cur;
	cur = head;
      }
    }
    else {
      last = cur;
      cur = cur->next;
    }
  }
}



int *consposvals;
int *conslenvals;
int byConsPosAndLen(const void *i, const void *j)
{

  if( consposvals[*((int *)i)] > consposvals[*((int *)j)] )
    return 1;
  if( consposvals[*((int *)i)] < consposvals[*((int *)j)] )
    return -1;
  if( conslenvals[*((int *)i)] > conslenvals[*((int *)j)] )
    return 1;

  return -1;

}

int consfindUniqueOrder( int *locs, int *lens, int num,
			 int *outLocs, int *outLens, int *order )
{

  // first let's get the ordering when sorted by ending position.
  int *sortOrder = new int [num];

  for( int i = 0; i < num; i++ ){
    sortOrder[i] = i;
    locs[i] += lens[i];
  }
  consposvals = locs;
  conslenvals = lens;
  qsort( sortOrder, num, sizeof(int), byConsPosAndLen );
  for( int i = 0; i < num; i++ )
    locs[i] -= lens[i];

  // now let's compute the output arrays
  int outNum = 0;
  outLocs[0] = locs[sortOrder[0]];
  outLens[0] = lens[sortOrder[0]];
  order[sortOrder[0]] = 0;

  // for the duration of the next loop, outNum will be the *index* of the
  // last loc recorded. after the loop, it should be incremented.
  for( int i = 1; i < num; i++ ){
    // if this loc/len pair is different than the last one, add it in.
    if( locs[sortOrder[i]] != outLocs[outNum] ||
	lens[sortOrder[i]] != outLens[outNum] ){
      outNum++;
      outLocs[outNum] = locs[sortOrder[i]];
      outLens[outNum] = lens[sortOrder[i]];
    }
    order[sortOrder[i]] = outNum;
  }

  // verify acceptability
#ifdef DEBUG
  for( int i = 0; i < num; i++ ){
    assert( locs[i] == outLocs[order[i]] );
    assert( lens[i] == outLens[order[i]] );
    if( 0 < i && i < outNum + 1 )
      assert( outLocs[i - 1] + outLens[i - 1] <=
	      outLocs[i] + outLens[i] );
  }
#endif

  delete[] sortOrder;
  return outNum + 1;

}


#define MATCH 0
#define GAP_BASE 1
#define GAP_SECOND 2

#define MAVID_MAX_NUM_MATCH 8000

int select_consistent( matchlist *maskingMatches,
		       int* &outBaseLocs, int* &outSecLocs, int* &outLens )
{

  // If we don't have any matches, return 0.
  if( maskingMatches == NULL )
    return 0;

  // Grab as many as possible up to MAVID_MAX_NUM_MATCH
  matchlist **best = new matchlist* [MAVID_MAX_NUM_MATCH];
  //  matchlist *best [MAVID_MAX_NUM_MATCH];
  int num = 0;
  matchlist *temp;

  temp = maskingMatches;
  while( num < MAVID_MAX_NUM_MATCH && temp != NULL ){
    best[num] = temp;
    num++;
    temp = temp->next;
  }

  // Now let's get the unique positions and the correct orderings
  // First for the base
  int *locs = new int [MAVID_MAX_NUM_MATCH], *lens = new int [MAVID_MAX_NUM_MATCH];
  int *baseLocs = new int [MAVID_MAX_NUM_MATCH], *baseLens = new int [MAVID_MAX_NUM_MATCH],
    *baseOrder = new int [MAVID_MAX_NUM_MATCH];
  /*int locs[MAVID_MAX_NUM_MATCH], lens[MAVID_MAX_NUM_MATCH];
  int baseLocs[MAVID_MAX_NUM_MATCH], baseLens[MAVID_MAX_NUM_MATCH],
  baseOrder[MAVID_MAX_NUM_MATCH];*/

  for( int i = 0; i < num; i++ ){
    locs[i] = best[i]->basePos;
    lens[i] = best[i]->len;
  }
  int baseNum = consfindUniqueOrder( locs, lens, num, baseLocs, baseLens, 
				     baseOrder );

  // And then for the second
  int *secLocs = new int [MAVID_MAX_NUM_MATCH], *secLens = new int [MAVID_MAX_NUM_MATCH],
    *secOrder = new int [MAVID_MAX_NUM_MATCH];
  /*int secLocs[MAVID_MAX_NUM_MATCH], secLens[MAVID_MAX_NUM_MATCH],
    secOrder[MAVID_MAX_NUM_MATCH];*/

  for( int i = 0; i < num; i++ ){
    locs[i] = best[i]->secPos;
  }
  int secNum = consfindUniqueOrder( locs, lens, num, secLocs, secLens, 
				secOrder );

  // first we need to make the score matrix, and the LNO(last non-overlapping)
  // arrays.

  int *baseLNO = new int [MAVID_MAX_NUM_MATCH];
  int *secLNO = new int [MAVID_MAX_NUM_MATCH];
  /*  int baseLNO[MAVID_MAX_NUM_MATCH];
      int secLNO[MAVID_MAX_NUM_MATCH];*/

  // as in the loops below, this is increased by one so that it represents
  // a position in the scoreMatrix(where things are increased by one due to
  // the boundary conditions).
  baseLNO[0] = secLNO[0] = 0;

  for( int i = 1; i < baseNum; i++ ){
    int j;
    for( j = i - 1; j >= 0 && (baseLocs[j] + baseLens[j] > baseLocs[i]); j-- )
      {}
    baseLNO[i] = j + 1;
  }
  for( int i = 1; i < secNum; i++ ){
    int j;
    for( j = i - 1; j >= 0 && (secLocs[j] + secLens[j] > secLocs[i]); j-- )
      {}
    secLNO[i] = j + 1;
  }

  float **scores = new float* [MAVID_MAX_NUM_MATCH];
  //  float scores[MAVID_MAX_NUM_MATCH][MAVID_MAX_NUM_MATCH];

  for( int i = 0; i < baseNum; i++ ){
    scores[i] = new float [MAVID_MAX_NUM_MATCH];
    for( int j = 0; j < secNum; j++ )
      scores[i][j] = 0;
  }
  for( int i = 0; i < num; i++ ){
    scores[baseOrder[i]][secOrder[i]] = 1;
  }

  // now perform an adapted Smith-Waterman procedure to find the maximum set
  // of consistent matches.
  float **scoreMatrix = new float* [MAVID_MAX_NUM_MATCH + 1];
  char **ptrMatrix = new char* [MAVID_MAX_NUM_MATCH];
  /*  float scoreMatrix[MAVID_MAX_NUM_MATCH + 1][MAVID_MAX_NUM_MATCH + 1];
      char ptrMatrix[MAVID_MAX_NUM_MATCH][MAVID_MAX_NUM_MATCH];*/

  for( int i = 0; i < baseNum + 1; i++ ){
    scoreMatrix[i] = new float [MAVID_MAX_NUM_MATCH + 1];
    if( i < baseNum )
      ptrMatrix[i] = new char [MAVID_MAX_NUM_MATCH];
    scoreMatrix[i][0] = 0;
  }
  /*  for( int i = 0; i < baseNum + 1; i++ )
      scoreMatrix[i][0] = 0;*/

  for( int i = 0; i < secNum + 1; i++ )
    scoreMatrix[0][i] = 0;

  for( int i = 1; i < baseNum + 1; i++ ){
    for( int j = 1; j < secNum + 1; j++ ){
      scoreMatrix[i][j] = scoreMatrix[baseLNO[i - 1]][secLNO[j - 1]] + 
	scores[i - 1][j - 1];
      ptrMatrix[i - 1][j - 1] = MATCH;
      // note that these ">="s are critical! otherwise it will match things
      // which are not matches!
      if( scoreMatrix[i - 1][j] >= scoreMatrix[i][j] ){
	scoreMatrix[i][j] = scoreMatrix[i - 1][j];
	ptrMatrix[i - 1][j - 1] = GAP_BASE;
      }
      if( scoreMatrix[i][j - 1] >= scoreMatrix[i][j] ){
	scoreMatrix[i][j] = scoreMatrix[i][j - 1];
	ptrMatrix[i - 1][j - 1] = GAP_SECOND;
      }
    }
  }

  outBaseLocs = new int [num];
  outSecLocs = new int [num];
  outLens = new int [num];
  int numMatches = 0;

  for( int i = baseNum - 1, j = secNum - 1; i >= 0 && j >= 0; ){
    if( ptrMatrix[i][j] == MATCH ){
      outBaseLocs[numMatches] = baseLocs[i];
      outSecLocs[numMatches] = secLocs[j];
      outLens[numMatches] = baseLens[i];
      numMatches++;
      // subtract one to account for the adding of one above.
      i = baseLNO[i] - 1;
      j = secLNO[j] - 1;
    }
    else if( ptrMatrix[i][j] == GAP_BASE ){
      i--;
    }
    else
      j--;
  }

  // reverse the arrays so that they are in ascending order.
  int mid = (numMatches >> 1);
  for( int i = 0; i < mid; i++ ){
    swap( outBaseLocs[i], outBaseLocs[numMatches - i - 1] );
    swap( outSecLocs[i], outSecLocs[numMatches - i - 1] );
    swap( outLens[i], outLens[numMatches - i - 1] );
  }

  delete[] best;
  delete[] locs;
  delete[] lens;
  delete[] baseLocs;
  delete[] baseLens;
  delete[] baseOrder;
  delete[] secLocs;
  delete[] secLens;
  delete[] secOrder;
  delete[] baseLNO;
  delete[] secLNO;

  for( int i = 0; i < baseNum + 1; i++ )
    delete[] scoreMatrix[i];
  for( int i = 0; i < baseNum; i++ ){
    delete[] ptrMatrix[i];
    delete[] scores[i];
  }
  delete[] scoreMatrix;
  delete[] ptrMatrix;
  delete[] scores;

  return numMatches;

}




void unionize( ConstraintsList *&list )
{

  if( !list )
    return;

  int num_anchors = 0;

  for( ConstraintsList *temp = list; temp != NULL; temp = temp->next )
    num_anchors++;

  ConstraintsList **lists = new ConstraintsList* [num_anchors];

  num_anchors = 0;
  for( ConstraintsList *temp = list; temp != NULL; temp = temp->next ){
    lists[num_anchors] = temp;
    num_anchors++;
  }

  qsort( lists, num_anchors, sizeof( ConstraintsList* ), cons_compare );
  list = lists[0];

  ConstraintsList *tail = list;

  for( int i = 0; i < num_anchors; i++ ){
    // if this guy has been deleted, skip it.
    if( !lists[i] )
      continue;
    // else, first extend it
    for( int j = i + 1; j < num_anchors; j++ ){
      if( !lists[j] )
	continue;
      // if they don't overlap in the first sequence, give up. none of the
      // later ones will either.
      if( lists[j]->cons.first_start >= lists[i]->cons.first_end )
	break;
      // if they overlap in the second sequence, extend i and delete j.
      if( lists[j]->cons.second_start >= lists[i]->cons.second_start &&
	  lists[j]->cons.second_start < lists[i]->cons.second_end ){
	lists[i]->cons.first_end = MAX( lists[i]->cons.first_end,
					lists[j]->cons.first_end );
	lists[i]->cons.second_end = MAX( lists[i]->cons.second_end,
					 lists[j]->cons.second_end );
	delete lists[j];
	lists[j] = NULL;
      }
    }

    // then add it
    tail->next = lists[i];
    tail = lists[i];
    tail->next = NULL;
  }

  delete[] lists;

}
