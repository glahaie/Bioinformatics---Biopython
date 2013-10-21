
#ifndef MAVID_MATCH
#define MAVID_MATCH

#include "sufftree.h"
#include "common.h"

struct matchlist {
  int basePos, secPos, len;
  matchlist *next;
};


void getMatchList( bp *baseSeq, char *baseMask, int baseLen,
		   bp *secSeq, char *secMask, int secLen,
		   int minLen, bool useMasking,
		   matchlist *&head, matchlist *&tail );

matchlist **filterMatches( matchlist *&head, int *baseStarts, int *baseEnds,
			   int *secStarts, int *secEnds,
			   int numAnchors );

void sortMatches( matchlist *&head, matchlist *&tail );

void clipMask( bp *seq, char *mask, int len,
	       bp *&preppedSeq, char *&cutMask, int *&shifts,
	       int &preppedLen );


void prepSeqs( bp *baseSeq, char *baseMask, int baseLen,
	       bp *secSeq, char *secMask, int secLen,
	       int minLen, bool useMasking,
	       bp *&preppedSeq, int *&nextMask, int *&baseShifts,
	       int *&secShifts, int &breakPt, int &preppedLen );

void findMatches( bp *baseSeq, char *baseMask, int baseLen,
		  bp *secSeq, char *secMask, int secLen, 
		  int minLen );

bool is_maximal_match( bp *baseSeq, char *baseMask, int baseLen,
		       bp *secSeq, char *secMask, int secLen,
                       bool useMasking, int basePos, int secPos, int len );

#endif
