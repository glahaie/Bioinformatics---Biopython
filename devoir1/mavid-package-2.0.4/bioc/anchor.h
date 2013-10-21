
#ifndef MAVID_ANCHOR
#define MAVID_ANCHOR

#include "match.h"
#include "common.h"

int select_anchors( bp *baseSeq, int baseLen, int baseFrom, int baseTo,
                    bp *secSeq, int secLen, int secFrom, int secTo,
                    matchlist* maskingMatches, matchlist* nomaskingMatches,
                    bool isShort,
                    int* &outBaseLocs, int* &outSecLocs, int* &outLens );

int findUniqueOrder( int *locs, int *lens, int num, int *outLocs, 
		     int *outLens, int *order );

int flankAlign( bp *baseSeq, bp *secSeq, int len );

void fixAnchor( int basePos, int secPos, int len,
		int *baseImg, int *secImg );

#endif
