#ifndef MAVID_REFINE_H
#define MAVID_REFINE_H

#include "alignments.h"
#include "malign.h"
#include "tree.h"
#include "constraints.h"


FastaAlign refine_alignment( FastaAlign al, PhyloTree *edge, MultiFastaSeq &seqs, PhyloTree *tree,
                             double **&profile, ConstraintsTree *cons, Map *map );

FastaAlign refine_leaves( FastaAlign al, MultiFastaSeq &seqs, PhyloTree *tree,
                          double **&profile, ConstraintsTree *cons, Map *map );

#endif
