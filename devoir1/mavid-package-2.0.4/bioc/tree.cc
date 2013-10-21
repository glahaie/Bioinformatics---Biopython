#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include "tree.h"
#include "seq.h"
#include "matrices.h"

#define MAX_LEN 30

void deletePhyloTree( PhyloTree *tree )
{

  for( int i = 0; i < tree->numChild; i++ )
    deletePhyloTree( tree->child[i] );
  if( tree->numChild > 0 ){
    delete[] tree->child;
    delete[] tree->lens;
  }
  delete tree;

}


void deleteNode( PhyloTree *tree )
{
  
  if( tree->numChild > 0 ){
    delete[] tree->child;
    delete[] tree->lens;
  }
  delete tree;

}


PhyloTree *copyPhyloTree( PhyloTree *tree )
{

  PhyloTree *toRet = new PhyloTree;

  if( tree->label ){
    toRet->label = new char [strlen(tree->label) + 1];
    strcpy( toRet->label, tree->label );
  }
  else {
    toRet->label = NULL;
  }

  toRet->index = tree->index;
  toRet->seq = tree->seq;
  toRet->parent = tree->parent;
  toRet->parentLen = tree->parentLen;
  toRet->numChild = tree->numChild;

  if( toRet->numChild > 0 ){
    toRet->child = new PhyloTree* [toRet->numChild];
    toRet->lens = new float [toRet->numChild];
    for( int i = 0; i < toRet->numChild; i++ ){
      toRet->child[i] = copyPhyloTree( tree->child[i] );
      toRet->child[i]->parent = toRet;
      toRet->child[i]->parentLen = tree->lens[i];
      toRet->lens[i] = tree->lens[i];
    }
  }

  return toRet;
  
}

PhyloTree *readPhyloTree( char *fileName )
{

  int strlen;
  char *string;
  FILE *in;
  PhyloTree *toRet;

  in = fopen( fileName, "r" );
  if( !in )
    return NULL;

  for( strlen = 0; !feof( in ); ){
    char c = getc( in );

    if( !isspace( c ) && c != EOF )
      strlen++;
  }
  string = new char [strlen + 1];

  rewind( in );
  for( strlen = 0; !feof( in ); ){
    char c = getc( in );

    if( !isspace( c ) && c != EOF )
      string[strlen++] = c;
  }
  string[strlen] = '\0';

  if( string[0] != '(' ){
    delete[] string;
    fclose( in );
    return NULL;
  }
  
  toRet = stringToPhyloTree( string + 1 );

  delete[] string;
  fclose( in );

  return toRet;

}

PhyloTree *stringToPhyloTree( char *string )
{

  PhyloTree *toRet = new PhyloTree;
  toRet->numChild = 1;
  toRet->label = NULL;
  toRet->index = -1;
  toRet->seq = NULL;
  toRet->parent = NULL;
  toRet->parentLen = 0;

  for( int depth = 0, i = 0; depth >= 0; i++ ){
    if( string[i] == '\0' )
      return NULL;
    if( string[i] == '(' )
      depth++;
    if( string[i] == ')' )
      depth--;
    if( string[i] == ',' && depth == 0 )
      toRet->numChild++;
  }

  toRet->child = new PhyloTree* [toRet->numChild];
  toRet->lens = new float [toRet->numChild];
  toRet->parent = NULL;

  toRet->numChild = 0;
  for( int i = 0; string[i] != ')'; ){
    if( string[i] == '(' ){
      toRet->child[toRet->numChild] = stringToPhyloTree( string + i + 1 );
      if( !toRet->child[toRet->numChild] ){
	for( int j = 0; j < toRet->numChild; j++ )
	  delete toRet->child[j];
	delete toRet;
	return NULL;
      }
      toRet->child[toRet->numChild]->parent = toRet;
      if( toRet->child[toRet->numChild] == NULL )
	return NULL;
      i++;
      for( int depth = 1; depth > 0; i++ ){
	if( string[i] == '(' )
	  depth++;
	if( string[i] == ')' )
	  depth--;
      }
      if( string[i] != ':' )
	return NULL;
      sscanf( string + i + 1, "%f", &toRet->lens[toRet->numChild] );
      toRet->child[toRet->numChild]->parentLen = toRet->lens[toRet->numChild];
      toRet->numChild++;
      while( string[i] != ',' && string[i] != ')' )
	i++;
    }
    else {
      int temp = i;

      while( string[i] != ':' && string[i] != ',' && string[i] != ')' )
	i++;
      if( string[i] != ':' || temp == i )
	return NULL;
      toRet->child[toRet->numChild] = new PhyloTree;
      toRet->child[toRet->numChild]->numChild = 0;
      toRet->child[toRet->numChild]->parent = toRet;
      toRet->child[toRet->numChild]->lens = NULL;
      toRet->child[toRet->numChild]->child = NULL;
      toRet->child[toRet->numChild]->index = -1;
      toRet->child[toRet->numChild]->seq = NULL;
      toRet->child[toRet->numChild]->label = new char [i - temp + 1];
      i = temp;
      while( string[i] != ':' ){
	toRet->child[toRet->numChild]->label[i - temp] = string[i];
	i++;
      }
	
      toRet->child[toRet->numChild]->label[i - temp] = '\0';
      sscanf( string + i + 1, "%f", &toRet->lens[toRet->numChild] );
      toRet->child[toRet->numChild]->parentLen = toRet->lens[toRet->numChild];
      toRet->numChild++;
      while( string[i] != ',' && string[i] != ')' )
	i++;
    }
    if( string[i] == ',' )
      i++;
  }

  return toRet;

}


void printPhyloTree( PhyloTree *tree )
{

  if( tree->numChild == 0 ){
    for( int i = 0; ; i++ ){
      if( isspace(tree->label[i]) || tree->label[i] == ',' || tree->label[i] == ':' || tree->label[i] == '\0' )
	return;
      printf( "%c", tree->label[i] );
    }
    return;
  }

  printf( "(" );
  printPhyloTree( tree->child[0] );
  printf( ":%f", tree->lens[0] );
  for( int i = 1; i < tree->numChild; i++ ){
    printf( "," );
    printPhyloTree( tree->child[i] );
    printf( ":%f", tree->lens[i] );
  }
  printf( ")" );
  if( tree->parent == NULL )
    printf( ";\n" );

}

bool isBinaryPhyloTree( PhyloTree *tree )
{

  if( !tree )
    return false;

  if( tree->numChild == 0 )
    return true;

  if( tree->numChild != 2 )
    return false;

  return isBinaryPhyloTree( tree->child[0] ) &&
    isBinaryPhyloTree( tree->child[1] );

}


int countPhyloTreeNodes( PhyloTree *tree )
{

  if( !tree )
    return 0;

  if( tree->numChild == 0 )
    return 1;

  int num = 0;

  for( int i = 0; i < tree->numChild; i++ )
    num += countPhyloTreeNodes( tree->child[i] );

  return num;

}


bool getPhyloTreeIndices( PhyloTree *tree, char **fastaLines, int numSeqs )
{

  if( tree->numChild == 0 ){
    unsigned int len = strlen( tree->label );

    tree->index = -1;
    for( int i = 0; i < numSeqs; i++ ){
      if( strlen( fastaLines[i] ) >= len &&
	  strncmp( tree->label, fastaLines[i], len ) == 0 &&
	  (len == MAX_LEN || fastaLines[i][len] == '\0' || fastaLines[i][len] == ' ') ){
	tree->index = i;
      }
    }
    if( tree->index == -1 )
      return false;
  }
    
  for( int i = 0; i < tree->numChild; i++ )
    if( !getPhyloTreeIndices( tree->child[i], fastaLines, numSeqs ) )
      return false;

  return true;

}


bool getPhyloTreeIndices( PhyloTree *tree, MultiFastaSeq &seq )
{

  if( tree->numChild == 0 ){
    unsigned int len = strlen( tree->label );

    tree->index = -1;
    for( int i = 0; i < seq.numSeqs; i++ ){
      if( strlen( seq.seqs[i].fastaLine ) >= len &&
	  strncmp( tree->label, seq.seqs[i].fastaLine, len ) == 0 &&
	  (len == MAX_LEN || seq.seqs[i].fastaLine[len] == '\0' || seq.seqs[i].fastaLine[len] == ' ') ){
	tree->index = i;
	tree->seq = &seq.seqs[i];
      }
    }
    if( tree->index == -1 )
      return false;
  }
    
  for( int i = 0; i < tree->numChild; i++ )
    if( !getPhyloTreeIndices( tree->child[i], seq ) )
      return false;

  return true;

}


int getLabels( char **labels, PhyloTree *tree )
{

  if( tree->numChild == 0 ){
    labels[0] = tree->label;
    return 1;
  }

  int numLabels = 0;

  for( int i = 0; i < tree->numChild; i++ )
    numLabels += getLabels( labels + numLabels, tree->child[i] );

  return numLabels;

}


void permuteMultiFasta( MultiFastaSeq &src, MultiFastaSeq &dest,
			PhyloTree *tree, int &curIdx )
{

  if( tree->numChild == 0 ){
    dest.seqs[curIdx] = src.seqs[tree->index];
    curIdx++;
  }
  else {
    permuteMultiFasta( src, dest, tree->child[0], curIdx );
    permuteMultiFasta( src, dest, tree->child[1], curIdx );
  }

}


float **leftProfs, **rightProfs;
int lastDepth = -1;

int mlRoot( PhyloTree *tree, char *leaves, float *profile, int curDepth,
	    int maxDepth )
{

  if( curDepth == 0 && maxDepth > lastDepth ){
    for( int i = 0; i < lastDepth; i++ ){
      delete[] leftProfs;
      delete[] rightProfs;
    }
    delete[] leftProfs;
    delete[] rightProfs;

    leftProfs = new float* [maxDepth];
    rightProfs = new float* [maxDepth];
    for( int i = 0; i < maxDepth; i++ ){
      leftProfs[i] = new float [4];
      rightProfs[i] = new float [4];
    }
    lastDepth = maxDepth;
  }

  if( tree->numChild == 0 ){
    if( leaves[0] < DNA_N ){
      for( int i = 0; i < 4; i++ )
	profile[i] = 0.0;
      profile[leaves[0]] = 1.0;
    }
    else {
      profile[0] = 0.25;
      profile[1] = 0.25;
      profile[2] = 0.25;
      profile[3] = 0.25;
    }

    return 1;
  }
  else {
    int numLeaves;

    numLeaves = mlRoot( tree->child[0], leaves, leftProfs[curDepth],
			curDepth + 1, maxDepth );
    numLeaves += mlRoot( tree->child[0], leaves + numLeaves,
			 leftProfs[curDepth], curDepth + 1, maxDepth );

    return numLeaves;
  }

}


// finds the pair of leaves in tree which are furthest apart.
void farthestPair( PhyloTree *tree, PhyloTree *&first, PhyloTree *&second,
                   float &dist )
{

  PhyloTree *deepest;
  float depth;

  farthestPairRecurse( tree, first, second, dist, deepest, depth );

}

// finds the pair of leaves in tree which are furthest apart as well as the
// deepest node in tree.
void farthestPairRecurse( PhyloTree *tree, PhyloTree *&first,
                          PhyloTree *&second, float &dist, PhyloTree *&deepest,
                          float &depth )
{

  if( tree->numChild == 0 ){
    first = NULL;
    second = NULL;
    dist = 0;
    deepest = tree;
    depth = 0;
    return;
  }

  PhyloTree *nextDeepest = NULL;
  float nextDepth;

  dist = -FLT_MAX;
  depth = nextDepth = -FLT_MAX;
  for( int i = 0; i < tree->numChild; i++ ){
    PhyloTree *tempFirst, *tempSecond;
    float tempDist;
    PhyloTree *tempDeepest;
    float tempDepth;

    farthestPairRecurse( tree->child[i], tempFirst, tempSecond, tempDist,
			 tempDeepest, tempDepth );

    tempDepth += tree->lens[i];
    if( tempDepth > depth ){
      if( depth > 0 ){
	nextDepth = depth;
	nextDeepest = deepest;
      }
      depth = tempDepth;
      deepest = tempDeepest;
    }
    else if( tempDepth > nextDepth ){
      nextDepth = tempDepth;
      nextDeepest = tempDeepest;
    }

    if( tempDist > dist ){
      first = tempFirst;
      second = tempSecond;
      dist = tempDist;
    }

  }

  if( depth + nextDepth >= dist ){
    first = deepest;
    second = nextDeepest;
    dist = depth + nextDepth;
  }

}


PhyloTree *midpointRoot( PhyloTree *tree )
{

  PhyloTree *a, *b;
  float dist;

  farthestPair( tree, a, b, dist );
  return midpointRoot( tree, a, b, dist );

}


PhyloTree *midpointRoot( PhyloTree *tree, PhyloTree *first, PhyloTree *second,
			 float dist )
{

  // first let's find the depth of the two leaves
  int firstDepth, secondDepth;

  firstDepth = 0;
  for( PhyloTree *temp = first; temp != NULL; temp = temp->parent )
    firstDepth++;
  secondDepth = 0;
  for( PhyloTree *temp = second; temp != NULL; temp = temp->parent )
    secondDepth++;

  PhyloTree *currFirst, *currSec;
  float firstLen, secondLen;

  for( firstLen = secondLen = 0, currFirst = first, currSec = second;
       (firstLen + currFirst->parentLen < dist/2.0) &&
	 (secondLen + currSec->parentLen < dist/2.0); ){
    if( secondDepth >= firstDepth ){
      secondLen += currSec->parentLen;
      currSec = currSec->parent;
      secondDepth--;
    }
    else {
      firstLen += currFirst->parentLen;
      currFirst = currFirst->parent;
      firstDepth--;
    }
  }

  PhyloTree *toRet;

  toRet = new PhyloTree;
  toRet->numChild = 2;
  toRet->child = new PhyloTree* [2];
  toRet->lens = new float [2];
  toRet->parent = NULL;
  toRet->parentLen = 0;

  if( firstLen + currFirst->parentLen > dist/2.0 ){
    toRet->child[0] = copyPhyloTree( currFirst );
    toRet->lens[0] = dist/2.0 - firstLen;
    toRet->child[0]->parent = toRet;
    toRet->child[0]->parentLen = toRet->lens[0];
    toRet->child[1] = root_at_node( currFirst->parent, currFirst );
    toRet->lens[1] = currFirst->parentLen - toRet->lens[0];
    toRet->child[1]->parent = toRet;
    toRet->child[1]->parentLen = toRet->lens[1];
  }
  else {
    toRet->child[0] = copyPhyloTree( currSec );
    toRet->lens[0] = dist/2.0 - secondLen;
    toRet->child[0]->parent = toRet;
    toRet->child[0]->parentLen = toRet->lens[0];
    toRet->child[1] = root_at_node( currSec->parent, currSec );
    toRet->lens[1] = currSec->parentLen - toRet->lens[0];
    toRet->child[1]->parent = toRet;
    toRet->child[1]->parentLen = toRet->lens[1];
  }

  simplifyTree( toRet );

  return toRet;

}



// Returns a tree rooted at distance dist along the edge to tree's parent.
PhyloTree *root_on_edge( PhyloTree *tree, float dist )
{

  PhyloTree *toRet;

  toRet = new PhyloTree;
  toRet->numChild = 2;
  toRet->child = new PhyloTree* [2];
  toRet->lens = new float [2];
  toRet->parent = NULL;
  toRet->parentLen = 0;
  toRet->index = tree->index;
  toRet->seq = tree->seq;
  if( tree->label ){
    toRet->label = new char [strlen(tree->label) + 1];
    strcpy( toRet->label, tree->label );
  }
  else {
    toRet->label = NULL;
  }

  toRet->child[0] = copyPhyloTree( tree );
  toRet->lens[0] = dist;
  toRet->child[0]->parent = toRet;
  toRet->child[0]->parentLen = dist;
  toRet->child[1] = root_at_node( tree->parent, tree );
  toRet->lens[1] = tree->parentLen - dist;
  toRet->child[1]->parent = toRet;
  toRet->child[1]->parentLen = toRet->lens[1];

  simplifyTree( toRet );

  return toRet;

}


// Root at a given node in a tree, optionally removing one of its children
// in the process
PhyloTree *root_at_node( PhyloTree *tree, PhyloTree *toRemove )
{

  PhyloTree *toRet = new PhyloTree;

  toRet->numChild = tree->numChild;
  toRet->index = tree->index;
  toRet->seq = tree->seq;
  if( tree->label ){
    toRet->label = new char [strlen(tree->label) + 1];
    strcpy( toRet->label, tree->label );
  }
  else {
    toRet->label = NULL;
  }

  if( tree->parent != NULL ){
    toRet->numChild++;
  }
  if( toRemove != NULL ){
    toRet->numChild--;
  }

  toRet->child = new PhyloTree* [toRet->numChild];
  toRet->lens = new float [toRet->numChild];
  for( int i = 0, curChild = 0; i < tree->numChild; i++ ){
    if( tree->child[i] != toRemove ){
      toRet->child[curChild] = copyPhyloTree( tree->child[i] );
      toRet->child[curChild]->parent = toRet;
      toRet->child[curChild]->parentLen = tree->lens[i];
      toRet->lens[curChild] = tree->lens[i];
      curChild++;
    }
  }

  if( tree->parent != NULL ){
    toRet->child[toRet->numChild - 1] = root_at_node( tree->parent, tree );
    toRet->child[toRet->numChild - 1]->parent = toRet;
    toRet->child[toRet->numChild - 1]->parentLen = tree->parentLen;
    toRet->lens[toRet->numChild - 1] = tree->parentLen;
  }

  simplifyTree( toRet );

  return toRet;

}


void simplifyTree( PhyloTree *tree )
{

  for( int i = 0; i < tree->numChild; i++ ){
    if( tree->child[i]->numChild == 1 ){
      PhyloTree *temp = tree->child[i]->child[0];

      tree->lens[i] += tree->child[i]->lens[0];
      deleteNode( tree->child[i] );
      tree->child[i] = temp;
    }
    simplifyTree( tree->child[i] );
  }

}


PhyloTree *extract_tree( PhyloTree *tree, char **labels, int count )
{

  // is this a leaf?
  if( tree->numChild == 0 ){
    // is this one of the leaves we're extracting? if so, return it.
    for( int i = 0; i < count; i++ ){
      if( !strcmp( labels[i], tree->label ) )
	return copyPhyloTree( tree );
    }
    // else, return a NULL pointer.
    return NULL;
  }

  // internal node
  PhyloTree **rets = new PhyloTree* [tree->numChild];
  PhyloTree *to_ret = new PhyloTree;

  to_ret->numChild = 0;
  for( int i = 0; i < tree->numChild; i++ ){
    rets[to_ret->numChild] = extract_tree( tree->child[i], labels, count );
    if( rets[to_ret->numChild] )
      to_ret->numChild++;
  }

  if( to_ret->numChild == 0 ){
    delete to_ret;
    return NULL;
  }
  if( to_ret->numChild == 1 ){
    delete to_ret;
    to_ret = rets[0];
    to_ret->parentLen += tree->parentLen;
    to_ret->parent = NULL;
    delete rets;
    return to_ret;
  }

  to_ret->child = new PhyloTree* [to_ret->numChild];
  to_ret->lens = new float [to_ret->numChild];
  for( int i = 0; i < to_ret->numChild; i++ ){
    to_ret->child[i] = rets[i];
    to_ret->child[i]->parent = tree;
    to_ret->lens[i] = rets[i]->parentLen;
  }
  to_ret->parent = NULL;
  to_ret->parentLen = tree->parentLen;
  return to_ret;

}



// Generates a random tree labeled with the fasta lines from seqs.
PhyloTree *random_tree( MultiFastaSeq &seqs )
{

  char **fasta_ptrs = new char* [seqs.numSeqs];

  srand(time(NULL));

  for( int i = 0; i < seqs.numSeqs; i++ ){
    fasta_ptrs[i] = seqs.seqs[i].fastaLine;
  }
  for( int i = 0; i < seqs.numSeqs; i++ ){
    int j = (int)((float)seqs.numSeqs*(float)rand()/((float)RAND_MAX + 1.0));

    swap( fasta_ptrs[i], fasta_ptrs[j] );
  }

  return random_tree( fasta_ptrs, seqs.numSeqs );

}


// Generates a random tree with the leaves *in the order given by fasta_ptrs*.
PhyloTree *random_tree( char **fasta_ptrs, int numSeqs )
{

  PhyloTree *to_ret = new PhyloTree;

  if( numSeqs == 1 ){
    to_ret->numChild = 0;
    to_ret->child = NULL;
    to_ret->label = fasta_ptrs[0];
    to_ret->parent = NULL;
  }
  else {
    int num_left = 1 + (int)((float)(numSeqs - 1)*(float)rand()/((float)RAND_MAX + 1.0));

    to_ret->numChild = 2;
    to_ret->parent = NULL;
    to_ret->child = new PhyloTree* [2];
    to_ret->lens = new float [2];
    to_ret->child[0] = random_tree( fasta_ptrs, num_left );
    to_ret->lens[0] = 1.0;
    to_ret->child[0]->parent = to_ret;
    to_ret->child[1] = random_tree( fasta_ptrs + num_left, numSeqs - num_left );
    to_ret->lens[1] = 1.0;
    to_ret->child[1]->parent = to_ret;
  }

  return to_ret;

}



PhyloTree **get_leaves( PhyloTree *tree )
{

  PhyloTree **leaves;
  int num;

  num = countPhyloTreeNodes( tree );
  leaves = new PhyloTree* [num];
  num = 0;
  get_leaves( tree, leaves, num );

  return leaves;

}


void get_leaves( PhyloTree *tree, PhyloTree **leaves, int &num )
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




double get_dist( PhyloTree *first, PhyloTree *sec )
{

  int firstDepth, secDepth;

  firstDepth = 0;
  for( PhyloTree *temp = first; temp != NULL; temp = temp->parent )
    firstDepth++;

  secDepth = 0;
  for( PhyloTree *temp = sec; temp != NULL; temp = temp->parent )
    secDepth++;

  double currDist = 0;
  PhyloTree *currFirst = first, *currSec = sec;

  while( currFirst != currSec ){
    if( firstDepth >= secDepth ){
      currDist += currFirst->parentLen;
      currFirst = currFirst->parent;
      firstDepth--;
    }
    else {
      currDist += currSec->parentLen;
      currSec = currSec->parent;
      secDepth--;
    }
  }

  return currDist;

}


// Returns an array of distances and fills the labels array with the leaf
// labels in the same order as the array.
double **get_dists( PhyloTree *tree, char **&labels, int num_leaves = 0 )
{

  PhyloTree **leaves;
  double **dists;

  // get the leaves of the tree
  leaves = get_leaves( tree );
  num_leaves = countPhyloTreeNodes( tree );

  // get the labels
  labels = new char* [num_leaves];
  for( int i = 0; i < num_leaves; i++ )
    labels[i] = leaves[i]->label;

  // now actually compute the distances.
  dists = new double* [num_leaves];
  for( int i = 0; i < num_leaves; i++ )
    dists[i] = new double [num_leaves];

  for( int i = 0; i < num_leaves - 1; i++ ){
    dists[i][i] = 0.0;
    for( int j = i + 1; j < num_leaves; j++ ){
      dists[i][j] = get_dist( leaves[i], leaves[j] );
      dists[j][i] = dists[i][j];
    }
  }

  delete[] leaves;
  return dists;

}


// Returns an array of distances with the leaves ordered according to the 
// fastaLines array.
double **get_ordered_dists( PhyloTree *tree, char **fastaLines, int num_seqs )
{

  PhyloTree **leaves, **reordered_leaves;
  double **dists;

  // get the leaves of the tree
  leaves = get_leaves( tree );
  if( num_seqs == 0 )
    num_seqs = countPhyloTreeNodes( tree );

  // reorder them so their order agrees with the order of the
  // fastaLines
  reordered_leaves = new PhyloTree* [num_seqs];
  for( int i = 0; i < num_seqs; i++ ){
    for( int j = 0; j < num_seqs; j++ ){
      if( strcmp( leaves[j]->label, fastaLines[i] ) == 0 ){
	reordered_leaves[i] = leaves[j];
	break;
      }
    }
  }

  // now actually compute the distances.
  dists = new double* [num_seqs];
  for( int i = 0; i < num_seqs; i++ )
    dists[i] = new double [num_seqs];

  for( int i = 0; i < num_seqs - 1; i++ ){
    dists[i][i] = 0.0;
    for( int j = i + 1; j < num_seqs; j++ ){
      dists[i][j] = get_dist( leaves[i], leaves[j] );
      dists[j][i] = dists[i][j];
    }
  }

  return dists;

}


double **get_ordered_dists( PhyloTree *tree, MultiFastaSeq &seqs )
{

  char **fastaLines = new char* [seqs.numSeqs];
  double **dists;

  for( int i = 0; i < seqs.numSeqs; i++ )
    fastaLines[i] = seqs.seqs[i].fastaLine;

  dists = get_ordered_dists( tree, fastaLines, seqs.numSeqs );

  delete[] fastaLines;
  return dists;

}


double **get_ordered_dists( PhyloTree *tree, FastaAlign al )
{

  char **fastaLines = new char* [al.numSeqs];
  double **dists;

  for( int i = 0; i < al.numSeqs; i++ )
    fastaLines[i] = al.seq[i]->fastaLine;

  dists = get_ordered_dists( tree, fastaLines, al.numSeqs );

  delete[] fastaLines;
  return dists;

}


void add_cond_matrices( PhyloTree *tree )
{

  if( tree->numChild == 0 )
    return;

  tree->cond_matrices = new double** [tree->numChild];
  for( int i = 0; i < tree->numChild; i++ ){
    transition_matrix_gapped( tree->cond_matrices[i], 100.0*tree->lens[i] );
    add_cond_matrices( tree->child[i] );
  }

}
