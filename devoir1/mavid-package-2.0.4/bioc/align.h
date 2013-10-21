
#ifndef MAVID_ALIGN_H
#define MAVID_ALIGN_H

#include "common.h"
#include "match.h"
#include "constraints.h"


double globalAlign( bp *baseSeq, bp *secSeq, char *baseMask, char *secMask,
		    int baseLen, int secLen,
		    bool noMasking,
                    double **scores, double gap_open, double gap_extend,
		    double **baseProf, double **secProf,
		    int **baseImg, int **secImg,
		    ConstraintsList *cons, WeakList *weak );


double globalAlign( bp *baseSeq, bp *secSeq, char *baseMask, char *secMask,
		    int baseLen, int secLen,
		    bool noMasking,
                    double **scores, double gap_open, double gap_extend,
		    double **baseProf, double **secProf,
		    int **baseImg, int **secImg,
		    WeakList *weak );

double matchedAlign( bp *baseSeq, bp *secSeq, char *baseMask, char *secMask,
		     int baseLen, int secLen,
		     int baseFrom, int baseTo, int secFrom, int secTo,
		     matchlist *maskedMatches, matchlist *unmaskedMatches,
		     double **scores, double gap_open, double gap_extend,
		     double **baseProf, double **secProf,
		     int *baseImg, int *secImg );

double basicAlign( bp *baseSeq, bp *secSeq, char *baseMask, char *secMask,
		   int baseFrom, int baseTo, int secFrom, int secTo,
		   double **scores, double gap_open, double gap_extend,
		   matchlist *maskedMatches, matchlist *unmaskedMatches,
		   double **baseProf, double **secProf,
		   int *baseImg, int *secImg );

double gappedAlign( bp *first, bp *second, int fromfirst, int tofirst,
		    int fromsecond, int tosecond,
                    double **scores, double gap_open, double gap_extend,
		    double **baseProf, double **secProf,
		    int *imagefirst, int *imagesecond );

bool verifyAlignment( int baseFrom, int baseTo, int secFrom, int secTo,
                      int *baseImg, int *secImg );

bool verifyUnwritten( int from, int to, int *img );

#endif
