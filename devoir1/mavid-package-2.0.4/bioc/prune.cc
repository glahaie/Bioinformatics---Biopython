#include "seq.h"
#include "alignments.h"
#include "tree.h"
#include "prune.h"


// calculates the likelhoods of the columns of al using the given
// tree. assumes that the rows of al are in the same order as the
// BFS of tree.
double *column_likelihoods( FastaAlign al, PhyloTree *tree )
{

  bp *col = new bp [al.numSeqs];
  double *to_ret;

  add_cond_matrices( tree );

  to_ret = new double [al.len];
  for( int p = 0; p < al.len; p++ ){
    for( int i = 0; i < al.numSeqs; i++ )
      col[i] = al.line[i][p];
    to_ret[p] = column_likelihood( col, tree );
  }

  return to_ret;

}


// uses the Felsenstein pruning algorithm to calculate the likelihood of a single
// column.
double column_likelihood( bp *col, PhyloTree *tree )
{

  double *message;

  prune( col, tree, message );

  // SHOULD MODIFY TO USE ARBITRARY/CORRECT EQUI DISTS
  return 1/5.0*( message[0] + message[1] + message[2] + message[3] + message[4] );

}


// recursively implements Felsenstein's pruning algorithm. fills message
// with the message to return to the current node's parent and returns
// the number of children of the current node.
int prune( bp *col, PhyloTree *tree, double *&message )
{

  message = new double [5];

  if( tree->numChild == 0 ){
    for( int i = 0; i < 5; i++ )
      message[i] = 0;
    if( col[0] == DNA_GAP )
      message[col[0] - 1] = 1.0;
    else
      message[col[0]] = 1.0;
    return 1;
  }

  int num_children0, num_children1;
  double *message0, *message1;

  num_children0 = prune( col, tree->child[0], message0 );
  num_children1 = prune( col + num_children0, tree->child[1], message1 );

  for( int i = 0; i < 5; i++ ){
    double factor0, factor1;

    factor0 = factor1 = 0;
    for( int k = 0; k < 5; k++ ){
      factor0 += tree->cond_matrices[0][i][k]*message0[k];
      factor1 += tree->cond_matrices[1][i][k]*message1[k];
    }
    message[i] = factor0*factor1;
  }

  return num_children0 + num_children1;

}
