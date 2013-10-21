#include "common.h"
#include "match.h"
#include <stdio.h>


#define MAX_DIST 30
//#define MAX_DIST -1

int matchcomp( const void *a, const void *b )
{

  if( (*((matchlist **)a))->basePos > (*((matchlist **)b))->basePos )
    return 1;
  else
    return -1;

}

void joinMatches( bp *baseSeq, char *baseMask, int baseLen,
		  bp *secSeq, char *secMask, int secLen,
		  matchlist *&head, matchlist *&tail )
{

  matchlist **diagheads = new matchlist* [baseLen + secLen + 1],
    **diagtails = new matchlist* [baseLen + secLen + 1];

  for( int i = 0; i < baseLen + secLen + 1; i++ )
    diagheads[i] = diagtails[i] = NULL;
  diagheads += baseLen;
  diagtails += baseLen;

  // separate into diagonals
  for( matchlist *temp = head; temp != NULL; ){
    int d = temp->secPos - temp->basePos;

    if( diagheads[d] ){
      diagtails[d]->next = temp;
      diagtails[d] = temp;
      temp = temp->next;
      diagtails[d]->next = NULL;
    }
    else {
      diagheads[d] = diagtails[d] = temp;
      temp = temp->next;
      diagtails[d]->next = NULL;
    }
  }

  // sort each diagonal
  matchlist **scratch = new matchlist* [baseLen];

  for( int i = -baseLen; i < secLen + 1; i++ ){
    int dcount = 0;

    for( matchlist *temp = diagheads[i]; temp != NULL; temp = temp->next )
      scratch[dcount++] = temp;

    if( dcount == 0 )
      continue;

    qsort( scratch, dcount, sizeof( matchlist* ), matchcomp );
    for( int j = 0; j < dcount; j++ ){
      if( j < dcount - 1 &&
	  scratch[j + 1]->basePos - (scratch[j]->basePos +
				     scratch[j]->len) <= MAX_DIST ){
	scratch[j + 1]->len += scratch[j + 1]->basePos - scratch[j]->basePos;
	scratch[j + 1]->basePos = scratch[j]->basePos;
	scratch[j + 1]->secPos = scratch[j]->secPos;
	delete scratch[j];
	scratch[j] = NULL;
      }
    }

    // rejoin the matches
    diagheads[i] = diagtails[i] = NULL;
    for( int j = 0, last = -1; j < dcount; j++ ){
      if( !scratch[j] )
	continue;
      diagtails[i] = scratch[j];
      if( last != -1 )
	scratch[last]->next = scratch[j];
      else
	diagheads[i] = scratch[j];

      last = j;
    }
    diagtails[i]->next = NULL;
  }

  // rejoin the diagonals
  head = tail = NULL;
  for( int i = -baseLen; i < secLen + 1; i++ ){
    if( !head ){
      head = diagheads[i];
      tail = diagtails[i];
    }
    else if( diagheads[i] ){
      tail->next = diagheads[i];
      tail = diagtails[i];
    }
  }

  diagheads -= baseLen;
  diagtails -= baseLen;
  delete[] diagheads;
  delete[] diagtails;
  delete[] scratch;

}


void makeWobbleSeq( bp *seq, char *mask, int len,
		    bp *&wseq, char *&wmask, int &wlen )
{

  wseq = new bp [2 * len + 2];
  wmask = new char [2 * len + 2];

  wlen = 0;
  for( int f = 0; f < 3; f++ ){
    for( int i = f; i < len; i++ ){
      if( (i - f) % 3 == 2 )
	continue;
      wseq[wlen] = seq[i];
      wmask[wlen++] = mask[i];
    }
    wseq[wlen] = BASE_N;
    wmask[wlen++] = MASKED;
  }

}


void getWobbleMatchList( bp *baseSeq, char *baseMask, int baseLen,
			 bp *secSeq, char *secMask, int secLen,
			 int minWLen,
			 matchlist *&head, matchlist *&tail )
{

  bp *baseWSeq, *secWSeq;
  char *baseWMask, *secWMask;
  int baseWLen, secWLen;

  makeWobbleSeq( baseSeq, baseMask, baseLen,
		 baseWSeq, baseWMask, baseWLen );
  makeWobbleSeq( secSeq, secMask, secLen,
		 secWSeq, secWMask, secWLen );

  getMatchList( baseWSeq, baseWMask, baseWLen,
		secWSeq, secWMask, secWLen,
		minWLen, true,
		head, tail );

  delete[] baseWSeq;
  delete[] secWSeq;
  delete[] baseWMask;
  delete[] secWMask;

  int baseBreak1, baseBreak2, secBreak1, secBreak2;
  matchlist *heads[3][3], *tails[3][3];

  for( int i = 0; i < 3; i++ )
    for( int j = 0; j < 3; j++ )
      heads[i][j] = tails[i][j] = NULL;

  baseBreak1 = 2 * (baseLen/3) + (baseLen % 3) + 1;
  baseBreak2 = baseBreak1 + 2 * ((baseLen - 1)/3) + ((baseLen - 1) % 3) + 1;
  secBreak1 = 2 * (secLen/3) + (secLen % 3) + 1;
  secBreak2 = secBreak1 + 2 * ((secLen - 1)/3) + ((secLen - 1) % 3) + 1;
  for( matchlist *temp = head; temp != NULL; ){
    int baseF, secF;

    if( temp->basePos >= baseBreak2 ){
      temp->basePos -= baseBreak2;
      baseF = 2;
    }
    else if( temp->basePos >= baseBreak1 ){
      temp->basePos -= baseBreak1;
      baseF = 1;
    }
    else {
      baseF = 0;
    }

    if( temp->secPos >= secBreak2 ){
      temp->secPos -= secBreak2;
      secF = 2;
    }
    else if( temp->secPos >= secBreak1 ){
      temp->secPos -= secBreak1;
      secF = 1;
    }
    else {
      secF = 0;
    }

    if( (temp->basePos + temp->len - 1) % 2 == 0 )
      temp->len--;
    if( temp->basePos % 2 == 1 ){
      temp->len--;
      temp->basePos++;
      temp->secPos++;
    }

    if( temp->secPos % 2 == 1 ){
      matchlist *next = temp->next;
      //      if( last != NULL )
      //	last->next = next;
      //      else
      //	head = next;
      //      if( next == NULL )
      //	tail = last;
      delete temp;
      temp = next;
    }
    else {
      temp->basePos = 3 * (temp->basePos / 2) + baseF;
      temp->secPos = 3 * (temp->secPos / 2) + secF;
      temp->len = 3 * (temp->len / 2);
      if( (temp->basePos + temp->len > baseLen) ||
	  (temp->secPos + temp->len > secLen ) ){
	temp->len--;
      }

      //      printf( "%u %u %u\n", temp->basePos, temp->secPos, temp->len );
      // comment out
      //      last = temp;
      if( heads[baseF][secF] ){
	tails[baseF][secF]->next = temp;
	tails[baseF][secF] = temp;
	temp = temp->next;
	tails[baseF][secF]->next = NULL;
      }
      else {
	heads[baseF][secF] = tails[baseF][secF] = temp;
	temp = temp->next;
	tails[baseF][secF]->next = NULL;
      }

      //      last = temp;
      //      temp = temp->next;
    }

  }


  head = tail = NULL;
  for( int i = 0; i < 3; i++ ){
    for( int j = 0; j < 3; j++ ){
      if( !heads[i][j] )
	continue;

      joinMatches( baseSeq, baseMask, baseLen,
		   secSeq, secMask, secLen,
		   heads[i][j], tails[i][j] );
      if( head ){
	tail->next = heads[i][j];
	tail = tails[i][j];
      }
      else {
	head = heads[i][j];
	tail = tails[i][j];
      }
    }
  }

  // REMOVE ME
  for( matchlist *temp = head; temp != NULL; temp = temp->next ){
    if( baseSeq[temp->basePos] != secSeq[temp->secPos] ){
      temp->basePos++;
      temp->secPos++;
      temp->len--;
    }
    if( baseSeq[temp->basePos + temp->len - 1] != secSeq[temp->secPos + temp->len - 1] )
      temp->len--;
  }

  sortMatches( head, tail );


#ifdef DEBUG
  for( matchlist *temp = head; temp != NULL; temp = temp->next ){
    assert( temp->basePos >= 0 );
    assert( temp->basePos + temp->len <= baseLen );
    assert( temp->secPos >= 0 );
    assert( temp->secPos + temp->len <= secLen );
    /*
    for( int i = 0; i < temp->len; i++ ){
      assert( ((i + 1) % 3 == 0) ||
	      (baseSeq[temp->basePos + i] == secSeq[temp->secPos + i]) );
	      }*/
  }
#endif

}
