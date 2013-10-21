#include "alignments.h"
#include "malign.h"
#include "tree.h"
#include "constraints.h"



FastaAlign refine_alignment( FastaAlign al, PhyloTree *edge, MultiFastaSeq &seqs, PhyloTree *tree,
			     double **&profile, ConstraintsTree *cons, Map *map )
{

  PhyloTree *rerooted_tree;

  // produce the new guide tree
  rerooted_tree = root_on_edge( edge, edge->parentLen/2 );

  // extract the old alignments
  FastaAlign al0, al1, to_ret;

  al0 = extract_alignment( al, seqs, rerooted_tree->child[0] );
  al1 = extract_alignment( al, seqs, rerooted_tree->child[1] );

  // align them to each other
  to_ret = mergeAlign( al0, al1, seqs, rerooted_tree,
		       cons, map );

  return to_ret;

}


FastaAlign refine_leaves( FastaAlign al, MultiFastaSeq &seqs, PhyloTree *tree,
			  double **&profile, ConstraintsTree *cons, Map *map )
{

  FastaAlign to_ret;
  int num_leaves = countPhyloTreeNodes( tree );
  PhyloTree **leaves = new PhyloTree* [num_leaves];

  num_leaves = 0;
  get_leaves( tree, leaves, num_leaves );

  to_ret = al;
  for( int i = 0; i < num_leaves; i++ ){
    FastaAlign temp;

    temp = refine_alignment( to_ret, leaves[i], seqs, tree, profile, cons, map );
    if( i > 0 ){
      for( int j = 0; j < to_ret.numSeqs; j++ )
	delete[] to_ret.line[j];
      delete[] to_ret.line;
      delete[] to_ret.seq;
    }
    to_ret = temp;
  }

  return to_ret;

}
