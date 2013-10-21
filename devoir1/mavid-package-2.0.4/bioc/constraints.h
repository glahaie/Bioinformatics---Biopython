#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include "fasta.h"
#include "seq.h"
#include "tree.h"
#include "alignments.h"
#include "common.h"
#include "match.h"
#include <stdio.h>


#define MAX_MISSING 3

struct Constraint {

  int first_index, second_index;
  int first_start, first_end, second_start, second_end;

};

struct ConstraintsList {

  Constraint cons;
  ConstraintsList *next;

};

struct ConstraintsTree {  

  // node data
  char *label;
  int index;
  ConstraintsTree *parent;
  ConstraintsList *list;

  // children data
  int numChild;
  ConstraintsTree **child;

};


struct Range {

  int start, end;

};

struct Map {

  int num, width;
  Range **segment;

};

enum WeakType {LEQ, GEQ};

struct Weak {

  int i,j;
  int i_index, j_index;
  WeakType type;

};

struct WeakList {

  Weak weak;
  WeakList *next;

};


void delete_constraints_tree( ConstraintsTree *tree );

ConstraintsTree *copy_tree( PhyloTree *tree );

void ancestralize_constraints( ConstraintsTree *cons, WeakList *weak,
			       FastaAlign al0, int *shift0, int len0,
			       FastaAlign al1, int *shift1, int len1,
			       MultiFastaSeq &seqs,
			       ConstraintsList *&a_cons, WeakList *&a_weak );

ConstraintsTree *phylo_to_constraints( PhyloTree *tree, ConstraintsTree *parent = NULL );

ConstraintsTree *make_constraintsTree( FILE *in, PhyloTree *tree );

Map *read_map_file( char *file_name );

WeakList *get_weak( Map *map, PhyloTree *tree, bool &non_overlapping,
		    bool &left_first );

void filter_by_weak( matchlist *&head, WeakList *weaks );

int select_consistent( matchlist *maskingMatches,
		       int* &outBaseLocs, int* &outSecLocs, int* &outLens );

void unionize( ConstraintsList *&list );

#endif
