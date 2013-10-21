#ifndef PRUNE_H
#define PRUNE_H

#include "alignments.h"
#include "tree.h"

double *column_likelihoods( FastaAlign al, PhyloTree *tree );

double column_likelihood( bp *col, PhyloTree *tree );

int prune( bp *col, PhyloTree *tree, double *&message );

#endif
