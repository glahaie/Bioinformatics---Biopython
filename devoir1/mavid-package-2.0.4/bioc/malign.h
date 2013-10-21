#ifndef MAVID_MALIGN_H
#define MAVID_MALIGN_H

#include "common.h"
#include "tree.h"
#include "alignments.h"
#include "constraints.h"


FastaAlign multiAlign( MultiFastaSeq &seqs, PhyloTree *tree,
		       ConstraintsTree *cons, Map *map );

FastaAlign mergeAlign( FastaAlign al0, FastaAlign al1, MultiFastaSeq &seqs, PhyloTree *tree,
		       ConstraintsTree *cons, Map *map );

void reconstruct( FastaAlign al, PhyloTree *tree, MultiFastaSeq &seqs,
                  bp *&seq, char *&mask, int &len, int *&shift,
                  double **&a_profile );

FastaAlign glue_alignments( FastaAlign al0, int *img0, int *shift0, int len0,
                            FastaAlign al1, int *img1, int *shift1, int len1 );

FastaAlign glue_alignments( FastaAlign al0, int *img0, int *shift0, int len0, float dist0,
                            FastaAlign al1, int *img1, int *shift1, int len1, float dist1 );

FastaAlign make_trivial_alignment( MultiFastaSeq &seqs, int index );

FastaAlign make_nonoverlapping_alignment( FastaAlign al0, FastaAlign al1, bool left_first );

void infer_profile( double *prof0, double **matrix0, double *prof1, double **matrix1, double *dest );

#endif
