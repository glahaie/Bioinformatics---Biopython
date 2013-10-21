
#ifndef BIOC_TREE_H
#define BIOC_TREE_H

struct PhyloTree;

#include "common.h"
#include "alignments.h"
#include "fasta.h"


struct PhyloTree {
  // node data
  char *label;
  FastaSeq *seq;
  int index;
  PhyloTree *parent;
  float parentLen;

  // children data
  int numChild;
  PhyloTree **child;
  float *lens;
  double ***cond_matrices;
};

void deletePhyloTree( PhyloTree *tree );

void deleteNode( PhyloTree *tree );

PhyloTree *copyPhyloTree( PhyloTree *tree );

PhyloTree *readPhyloTree( char *fileName );

PhyloTree *stringToPhyloTree( char *string );

void printPhyloTree( PhyloTree *tree );

bool isBinaryPhyloTree( PhyloTree *tree );

int countPhyloTreeNodes( PhyloTree *tree );

bool getPhyloTreeIndices( PhyloTree *tree, char **fastaLines, int numSeqs );

bool getPhyloTreeIndices( PhyloTree *tree, MultiFastaSeq &seq );

int getLabels( char **labels, PhyloTree *tree );

int mlRoot( PhyloTree *tree, char *leaves, float *profile );

void permuteMultiFasta( MultiFastaSeq &src, MultiFastaSeq &dest,
                        PhyloTree *tree, int &curIdx );

void farthestPair( PhyloTree *tree, PhyloTree *&first, PhyloTree *&second,
                   float &dist );

void farthestPairRecurse( PhyloTree *tree, PhyloTree *&first,
                          PhyloTree *&second, float &dist, PhyloTree *&deepest,
                          float &depth );

PhyloTree *midpointRoot( PhyloTree *tree );

PhyloTree *midpointRoot( PhyloTree *tree, PhyloTree *first, PhyloTree *second,
                         float dist );

PhyloTree *root_on_edge( PhyloTree *tree, float dist );

PhyloTree *root_at_node( PhyloTree *tree, PhyloTree *toRemove );

void simplifyTree( PhyloTree *tree );

PhyloTree *extract_tree( PhyloTree *tree, char **labels, int count );

PhyloTree *random_tree( MultiFastaSeq &seqs );

PhyloTree *random_tree( char **fasta_ptrs, int numSeqs );

PhyloTree **get_leaves( PhyloTree *tree );

void get_leaves( PhyloTree *tree, PhyloTree **leaves, int &num );

double get_dist( PhyloTree *first, PhyloTree *sec );

double **get_dists( PhyloTree *tree, char **&labels, int num_leaves );

double **get_ordered_dists( PhyloTree *tree, char **fastaLines, int num_seqs = 0 );

double **get_ordered_dists( PhyloTree *tree, MultiFastaSeq &seqs );

double **get_ordered_dists( PhyloTree *tree, FastaAlign al );

void add_cond_matrices( PhyloTree *tree );

#endif
