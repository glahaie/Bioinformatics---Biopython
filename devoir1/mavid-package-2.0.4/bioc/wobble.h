
#ifndef WOBBLE_H
#define WOBBLE_H

#include "common.h"

void joinMatches( bp *baseSeq, char *baseMask, int baseLen,
		  bp *secSeq, char *secMask, int secLen,
		  matchlist *head );

void makeWobbleSeq( bp *seq, char *mask, int len,
		    bp *&wseq, char *&wmask, int &wlen );

void getWobbleMatchList( bp *baseSeq, char *baseMask, int baseLen,
			 bp *secSeq, char *secMask, int secLen,
			 int minWLen,
			 matchlist *&head, matchlist *&tail );

#endif
