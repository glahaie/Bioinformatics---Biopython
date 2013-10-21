
#include "common.h"
#include "match.h"
#include "sufftree.h"

int *sortvals;
int byVal(const void *i, const void *j)
{

  if( sortvals[*((int *)i)] > sortvals[*((int *)j)] )
    return 1;
  if( sortvals[*((int *)i)] < sortvals[*((int *)j)] )
    return -1;

  return 0;
}





// given a list of matches and three arrays corresponding to a set of anchors
// this function filters those matches into an array of lists where each
// list corresponds to all the matches which falls between a pair of anchors
// (or before the first anchors or after the last). it returns that array
// of lists.
matchlist **filterMatches( matchlist *&head, int *baseStarts, int *baseEnds,
			   int *secStarts, int *secEnds,
			   int numAnchors )
{

  // there will be numAnchors + 1 regions between the anchors so we will
  // need numAnchors + 1 lists of matches.
  matchlist **heads = new matchlist* [numAnchors + 1];
  matchlist **tails = new matchlist* [numAnchors + 1];
  
  // initialize the lists to be empty
  for( int i = 0; i < numAnchors + 1; i++ ){
    heads[i] = NULL;
    tails[i] = NULL;
  }

  // now go through one by one and filter out each match.
  while( head != NULL ){
    // save the next match pointer for later.
    matchlist *next = head->next;
    int loc;

    // locate the limiting anchor(the next anchor which comes after this
    // position) in baseStarts.
    // this should be replaced by binary search later.
    for( loc = 0; loc < numAnchors && baseStarts[loc] < head->basePos; loc++ )
      {}
    //    loc = binSearch( baseStarts, numAnchors, head->basePos );

    // does the match given by head fit between the limiting anchors?
    bool isBeforeFirst, isBetween, isAfterLast;

    isBeforeFirst = (loc == 0) &&
      (head->basePos <= baseStarts[0]) && (head->secPos <= secStarts[0]);
    // if loc == numAnchors then this position came after the last anchor.
    isAfterLast = (loc == numAnchors);
    isBetween = (loc != 0) && (loc != numAnchors) &&
      (baseStarts[loc - 1] <= head->basePos) &&
      (head->basePos <= baseStarts[loc]) &&
      (secStarts[loc - 1] <= head->secPos) &&
      (head->secPos <= secStarts[loc]);

    // first let's chop the match as needed
    if( isBeforeFirst )
      head->len = MIN( head->len, baseStarts[0] - head->basePos,
		       secStarts[0] - head->secPos );

    if( isAfterLast ){
      int startChop = MAX( 0, baseEnds[loc - 1] - head->basePos + 1,
			   secEnds[loc - 1] - head->secPos + 1 );

      head->basePos += startChop;
      head->secPos += startChop;
      head->len -= startChop;
    }

    if( isBetween ){
      int startChop = MAX( 0, baseEnds[loc - 1] - head->basePos + 1,
			   secEnds[loc - 1] - head->secPos + 1 );

      head->basePos += startChop;
      head->secPos += startChop;
      head->len -= startChop;
      head->len = MIN( head->len, baseStarts[loc] - head->basePos,
		       secStarts[loc] - head->secPos );
    }

    // if any of the conditions are satisfied and the match is still long
    // enough then we should save the match in the appropriate list.
    if( (isBeforeFirst || isBetween || isAfterLast) && head->len >= 5 ){
      // is it the first entry in this list?
      if( heads[loc] == NULL )
	heads[loc] = tails[loc] = head;
      else {
	tails[loc]->next = head;
	tails[loc] = head;
      }
      // terminate the list.
      tails[loc]->next = NULL;
    }
    // if the match is not consistent with the anchors, just throw it away.
    else
      delete head;

    // advance to the next match.
    head = next;

  }

  for( int i = 0; i < numAnchors + 1; i++ )
    sortMatches( heads[i], tails[i] );

  delete[] tails;
  return heads;

}

// performs a very simple sort-by-length(a bucket sort?) on a list of matches.
void sortMatches( matchlist *&head, matchlist *&tail )
{

  int maxLen = 0;

  if( head == NULL )
    return;

  // first let's prep all the lists we'll need.
  for( matchlist *temp = head; temp != NULL; temp = temp->next )
    if( temp->len > maxLen )
      maxLen = temp->len;

  matchlist **heads = new matchlist* [maxLen + 1],
    **tails = new matchlist* [maxLen + 1];

  for( int i = 0; i <= maxLen; i++ )
    heads[i] = tails[i] = NULL;

  // now we seperate our big list out into a bunch of smaller lists.
  for( matchlist *temp = head; temp != NULL; ){
    matchlist *next = temp->next;

    if( tails[temp->len] == NULL )
      heads[temp->len] = temp;
    else
      tails[temp->len]->next = temp;
    tails[temp->len] = temp;
    tails[temp->len]->next = NULL;
    
    temp = next;
  }

  // now let's link them all together.
  int i = maxLen;

  head = heads[i];
  tail = tails[i];

  // go to the next non-empty list.
  while( i > 0 && heads[--i] == NULL )
    {}


  // we use '>' instead of '>=' because i gets decremented first thing in the
  // loop.
  while( i > 0 ){
    tail->next = heads[i];
    tail = tails[i];

    // go to the next non-empty list.
    while( heads[--i] == NULL && i > 0 )
      {}

  }
  tail->next = NULL;

  delete[] heads;
  delete[] tails;

}


void prepSeqs( bp *baseSeq, char *baseMask, int baseLen,
	       bp *secSeq, char *secMask, int secLen,
	       int minLen, bool useMasking,
	       bp *&preppedSeq, int *&nextMask, int *&baseShifts,
	       int *&secShifts, int &breakPt, int &preppedLen )
{

  preppedSeq = new bp [baseLen + secLen + 1];
  nextMask = new int [baseLen + secLen + 1];
  baseShifts = new int [baseLen];
  secShifts = new int [secLen];
  preppedLen = 0;

  for( int i = 0, shift = 0; i < baseLen; ){
    // skip over any masked parts
    if( (useMasking && baseMask[i] == MASKED) || baseSeq[i] == BASE_N ){
      preppedSeq[preppedLen] = BASE_N;
      baseShifts[preppedLen] = shift;
      nextMask[preppedLen] = preppedLen;
      preppedLen++; i++;
    }
    while( i < baseLen && ((useMasking && baseMask[i] == MASKED) ||
			   baseSeq[i] == BASE_N) ){
      shift++;
      i++;
    }
    // copy unmasked part, but record where we begin in case it ends up being
    // too short.
    int unmaskedLen = 0;

    while( i < baseLen && (!useMasking || baseMask[i] == UNMASKED) &&
	   baseSeq[i] != BASE_N ){
      preppedSeq[preppedLen] = baseSeq[i];
      baseShifts[preppedLen] = shift;
      preppedLen++; i++; unmaskedLen++;
    }

    for( int j = preppedLen - unmaskedLen; j < preppedLen; j++ )
      nextMask[j] = preppedLen;

    if( unmaskedLen < minLen ){
      preppedLen -= unmaskedLen ;
      shift += unmaskedLen;
    }
    
  }

  preppedSeq[preppedLen] = BASE_N;
  nextMask[preppedLen] = preppedLen;
  breakPt = preppedLen;
  preppedLen++;

  for( int i = 0, shift = 0; i < secLen; ){
    // skip over any masked parts
    if( (useMasking && secMask[i] == MASKED) || secSeq[i] == BASE_N ){
      preppedSeq[preppedLen] = BASE_N;
      secShifts[preppedLen - breakPt - 1] = shift;
      nextMask[preppedLen] = preppedLen;
      preppedLen++; i++;
    }
    while( i < secLen && ((useMasking && secMask[i] == MASKED) ||
			  secSeq[i] == BASE_N) ){
      shift++;
      i++;
    }
    // copy unmasked part, but record where we begin in case it ends up being
    // too short.
    int unmaskedLen = 0;

    while( i < secLen && (!useMasking || secMask[i] == UNMASKED) &&
	   secSeq[i] != BASE_N ){
      preppedSeq[preppedLen] = secSeq[i];
      secShifts[preppedLen - breakPt - 1] = shift;
      preppedLen++; i++; unmaskedLen++;
    }

    for( int j = preppedLen - unmaskedLen; j < preppedLen; j++ )
      nextMask[j] = preppedLen;

    if( unmaskedLen < minLen ){
      preppedLen -= unmaskedLen;
      shift += unmaskedLen;
    }
  }

#ifdef DEBUG
  // sanity check the results
  assert( preppedLen > 0 );
  for( int i = 0; i < preppedLen; i++ ){
    assert( 0 <= preppedSeq[i] && preppedSeq[i] <= 4 );
    //    assert( i == 0 || (preppedSeq[i] != BASE_N || preppedSeq[i - 1] != BASE_N) );
  }
#endif

}


int matchcompare( const void *m1, const void *m2 )
{

  match mm1 = *((match *)m1), mm2 = *((match *)m2);

  if( mm1.basePos > mm2.basePos )
    return 1;
  else if( mm1.basePos < mm2.basePos )
    return -1;
  else {
    if( mm1.secPos > mm2.secPos )
      return 1;
    else
      return -1;
  }
}


void getMatchList( bp *baseSeq, char *baseMask, int baseLen,
		   bp *secSeq, char *secMask, int secLen,
		   int minLen, bool useMasking,
		   matchlist *&head, matchlist *&tail )
{

  bp *preppedSeq;
  int *nextMask, *baseShifts, *secShifts;
  int breakPt, preppedLen;

  prepSeqs( baseSeq, baseMask, baseLen,
	    secSeq, secMask, secLen,
	    minLen, useMasking,
	    preppedSeq, nextMask, baseShifts, secShifts,
	    breakPt, preppedLen );

  char *transSeq = new char [preppedLen + 2];

  for( int i = 0; i < preppedLen; i++ )
    transSeq[i] = ALPHA[(int)preppedSeq[i]];
  transSeq[preppedLen] = '$';
  transSeq[preppedLen + 1] = '\0';

  ilist **junk;
  mlist matches;
  SuffixTree *t = new SuffixTree( transSeq, "ACGTN$" );

  junk = getMatches( *t->root, preppedSeq, nextMask, breakPt, minLen, 0, true,
		     matches );

  delete t;
  delete[] transSeq;
  delete[] preppedSeq;
  delete[] nextMask;

  int numMatch = 0;

  for( mitem *temp = matches.head; temp != NULL; temp = temp->next ){
    numMatch++;
    temp->data.basePos += baseShifts[temp->data.basePos];
    temp->data.secPos += secShifts[temp->data.secPos];
  }

#ifdef DEBUG

  for( mitem *temp = matches.head; temp != NULL; temp = temp->next ){
    assert( is_maximal_match( baseSeq, baseMask, baseLen, secSeq, secMask, secLen,
			      useMasking, temp->data.basePos, temp->data.secPos, temp->data.len ) );
  }

#endif

  head = tail = NULL;
  for( mitem *temp = matches.head; temp != NULL; ){
    int bP = temp->data.basePos, sP = temp->data.secPos, len = temp->data.len;
    matchlist *newMatch = new matchlist;
    mitem *next = temp->next;

    newMatch->basePos = bP;
    newMatch->secPos = sP;
    newMatch->len = len;
    newMatch->next = NULL;
    if( tail == NULL )
      head = tail = newMatch;
    else {
      tail->next = newMatch;
      tail = newMatch;
    }
    delete temp;
    temp = next;

    //    cout << newMatch->basePos << " " << newMatch->secPos << " "
    //	 << newMatch->len << endl;

  }

  sortMatches( head, tail );

  delete[] baseShifts;
  delete[] secShifts;

}


bool is_maximal_match( bp *baseSeq, char *baseMask, int baseLen,
		       bp *secSeq, char *secMask, int secLen,
		       bool useMasking, int basePos, int secPos, int len )
{

  // check that it's a match
  for( int k = 0; k < len; k++ ){
    if( baseSeq[basePos + k] != secSeq[secPos + k] ||
	(useMasking && baseMask[basePos + k] != UNMASKED) ||
	(useMasking && secMask[secPos + k] != UNMASKED) ||
	baseSeq[basePos + k] == BASE_N )
      return false;
  }
  // check left maximality
  if( (basePos != 0) && (secPos != 0) &&
      (baseSeq[basePos - 1] != BASE_N) &&
      (baseSeq[basePos - 1] == secSeq[secPos - 1]) &&
      (!useMasking || (baseMask[basePos - 1] == UNMASKED)) &&
      (!useMasking || (secMask[secPos - 1] == UNMASKED)) )
    return false;
  // check right maximality
  if( (basePos + len != baseLen) && (secPos + len != secLen) &&
      (baseSeq[basePos + len] != BASE_N) &&
      (baseSeq[basePos + len] == secSeq[secPos + len]) &&
      (!useMasking || (baseMask[basePos + len] == UNMASKED)) &&
      (!useMasking || (secMask[secPos + len] == UNMASKED)) )
    return false;

  return true;

}

