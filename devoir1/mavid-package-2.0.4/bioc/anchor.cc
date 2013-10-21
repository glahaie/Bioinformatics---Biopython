
#include "common.h"
#include "anchor.h"
#include "align.h"

#define MATCH 0
#define BASESEQ_GAP 1
#define SECSEQ_GAP 2

#define MAVID_MAX_NUM_MATCH 900
#define MAVID_MIN_ANCHOR_RATIO .50
#define MAVID_MIN_LEN_RATIO 0.2

#define FLANK_LEN 8


int by_end_basepos(const void *a, const void *b)
{

  if( (*(matchlist **)a)->basePos + (*(matchlist **)a)->len > (*(matchlist **)b)->basePos + (*(matchlist **)b)->len )
    return 1;
  else
    return -1;

}


// given a list of matches, this function finds the set of consistent matches
// of greatest total length among the first (up to) MAVID_MAX_NUM_MATCH matches
// in that list(so we should have the list sorted according to decreasing
// length). it returns the number of anchors it chooses, and sets and fills
// outBaseLocs, outSecLocs, and outLens with the information about the 
// matches. if isShort is true, then it requires that the total length of the
// anchor set it find is at least MAVID_MIN_ANCHOR_RATIO times the length of
// the shorter sequence. if it is not, then no anchors are returned.
int select_anchors( bp *baseSeq, int baseLen, int baseFrom, int baseTo,
		    bp *secSeq, int secLen, int secFrom, int secTo,
		    matchlist* maskingMatches, matchlist* nomaskingMatches,
		    bool isShort,
		    int* &outBaseLocs, int* &outSecLocs, int* &outLens )
{

  // If we don't have any matches, return 0.
  if( maskingMatches == NULL && nomaskingMatches == NULL )
    return 0;

  // Grab as many as possible up to MAVID_MAX_NUM_MATCH
  matchlist **best = new matchlist* [MAVID_MAX_NUM_MATCH];
  int num = 0;
  matchlist *temp;

  temp = maskingMatches;
  while( num < MAVID_MAX_NUM_MATCH && temp != NULL ){
    best[num] = temp;
    num++;
    temp = temp->next;
  }

  // if we've run out of masking matches and still don't have enough, move
  // on to the nomasking matches.
  temp = nomaskingMatches;
  while( num < MAVID_MAX_NUM_MATCH && temp != NULL ){
    best[num] = temp;
    num++;
    temp = temp->next;
  }


  /*  for( int i = 0; i < num; i++ )
    cout << best[i]->basePos << " " << best[i]->secPos << " " << best[i]->len << endl;
  cout << "done" << endl;
  exit(0);*/


  // The first thing we'll do is extend the matches.
  for( int i = 0; i < num; i++ ){
    while( (best[i]->basePos > baseFrom) && (best[i]->secPos > secFrom) &&
	   (baseSeq[best[i]->basePos - 1] == secSeq[best[i]->secPos - 1]) ){
      best[i]->basePos--;
      best[i]->secPos--;
      best[i]->len++;
    }
    while( (best[i]->basePos + best[i]->len < baseTo) &&
	   (best[i]->secPos + best[i]->len < secTo) &&
	   (baseSeq[best[i]->basePos + best[i]->len] ==
	    secSeq[best[i]->secPos + best[i]->len]) )
      best[i]->len++;
  }

  qsort( best, num, sizeof(matchlist *), by_end_basepos );

  /*  for( int i = 0; i < num; i++ )
    cout << best[i]->basePos << " " << best[i]->secPos << " " << best[i]->len << endl;
    exit(0);*/

  // calculate the score for each match
  double *scores = new double [num];

  for( int i = 0; i < num; i++ ){
    double leftBonus, rightBonus;
    int leftLen, rightLen;

    leftLen = MIN( best[i]->basePos, best[i]->secPos, FLANK_LEN );
    rightLen = MIN( baseLen - best[i]->basePos - best[i]->len,
		    secLen - best[i]->secPos - best[i]->len,
		    FLANK_LEN );

    leftBonus = flankAlign( baseSeq + best[i]->basePos - leftLen,
			    secSeq + best[i]->secPos - leftLen,
			    leftLen );
    rightBonus = flankAlign( baseSeq + best[i]->basePos + best[i]->len,
			     secSeq + best[i]->secPos + best[i]->len,
			     rightLen );

    scores[i] = 
      ((double)best[i]->len * (double)best[i]->len) +
      (double)(leftBonus + rightBonus);
  }

  // now find the maximal score of a consistent anchor set ending at each match.
  // it should be noted that this code depends rather heavily on the score of each match
  // being positive (that's not such an unreasonable requirement though).
  double *max_score = new double [num], best_score = 0.0;
  int *prev_ptr = new int [num], best_ptr = -1;

  for( int i = 0; i < num; i++ ){
    double best_extended_score = 0.0;

    // intialize the prev_ptr to be nothing and initialize the score to be the score of this match by itself.
    prev_ptr[i] = -1;
    max_score[i] = scores[i];

    // scan back and look for anchor sets it can extend
    for( int j = 0; j < i; j++ ){
      // can we extend this one? and if so is it better than the current best?
      if( best[j]->basePos + best[j]->len <= best[i]->basePos &&
	  best[j]->secPos + best[j]->len <= best[i]->secPos &&
	  max_score[j] >= best_extended_score ){
	prev_ptr[i] = j;
	best_extended_score = max_score[j];
      }
    }
    max_score[i] += best_extended_score;

    if( max_score[i] >= best_score ){
      best_ptr = i;
      best_score = max_score[i];
    }
  }


  int num_matches = 0;

  // find maximum anchor length and save only those anchors whose length
  // is above the required ratio. num_matches will be too high after this
  // but it gives us an upper bound.
  int max_len = 0, min_len = 0;

  for( int i = best_ptr; i != -1; i = prev_ptr[i] ){
    if( best[i]->len > max_len )
      max_len = best[i]->len;
    num_matches++;
  }

  min_len = (int)(((float)max_len)*MAVID_MIN_LEN_RATIO);
  outBaseLocs = new int [num_matches];
  outSecLocs = new int [num_matches];
  outLens = new int [num_matches];
  num_matches = 0;
  for( int i = best_ptr; i != -1; i = prev_ptr[i] ){
    if( best[i]->len >= min_len ){
      // REMOVE!
      //	cout << "accepted: " << baseLocs[i] << " " << secLocs[j] << " " <<  baseLens[i] << endl;
      outBaseLocs[num_matches] = best[i]->basePos;
      outSecLocs[num_matches] = best[i]->secPos;
      outLens[num_matches] = best[i]->len;
      num_matches++;
    }
  }

  // make sure the total length of our anchors is above the minimum.
  int totalLen;

  totalLen = 0;
  for( int i = 0; i < num_matches; i++ )
    totalLen += outLens[i];

  delete[] best;
  delete[] max_score;
  delete[] prev_ptr;

  if( (isShort == true) &&
      (totalLen < MAVID_MIN_ANCHOR_RATIO *
       (double)MIN( baseTo - baseFrom, secTo - secFrom )) ){
    delete[] outBaseLocs;
    delete[] outSecLocs;
    delete[] outLens;
    return 0;
  }

  // reverse the arrays so that they are in ascending order.
  int mid = (num_matches >> 1);
  for( int i = 0; i < mid; i++ ){
    swap( outBaseLocs[i], outBaseLocs[num_matches - i - 1] );
    swap( outSecLocs[i], outSecLocs[num_matches - i - 1] );
    swap( outLens[i], outLens[num_matches - i - 1] );
  }

  return num_matches;

}


int *posvals;
int *lenvals;
int byPosAndLen(const void *i, const void *j)
{

  if( posvals[*((int *)i)] > posvals[*((int *)j)] )
    return 1;
  if( posvals[*((int *)i)] < posvals[*((int *)j)] )
    return -1;
  if( lenvals[*((int *)i)] > lenvals[*((int *)j)] )
    return 1;

  return -1;

}


#define FLANK_MATCH 1
#define FLANK_MISMATCH -1
#define FLANK_GAP -1

int *flank_this_col = NULL, *flank_last_col = NULL;

int flankAlign( bp *baseSeq, bp *secSeq, int len )
{

  int *temp;

  if( flank_this_col == NULL ){
    flank_this_col = new int [FLANK_LEN + 1];
    flank_last_col = new int [FLANK_LEN + 1];
  }

  flank_last_col[0] = 0;
  for( int i = 1; i < FLANK_LEN + 1; i++ ){
    flank_last_col[i] = flank_last_col[i - 1] + FLANK_GAP;
  }

  // adjust the sequence ptrs to be one-based
  baseSeq--; secSeq--;

  for( int i = 1; i <= len; i++ ){
    flank_this_col[0] = flank_last_col[0] + FLANK_GAP;
    for( int j = 1; j <= len; j++ ){
      int score;

      score = MAX( flank_last_col[j] + FLANK_GAP,
		   flank_this_col[j - 1] + FLANK_GAP,
		   flank_last_col[j - 1] + ((baseSeq[i] == secSeq[j]) ?
				      FLANK_MATCH: FLANK_MISMATCH ));
      flank_this_col[j] = score;
    }
    temp = flank_this_col;
    flank_this_col = flank_last_col;
    flank_last_col = temp;
  }

  return flank_last_col[len];

}


void fixAnchor( int basePos, int secPos, int len,
		int *baseImg, int *secImg )
{

  for( int k = 0; k < len; k++ ){
#ifdef DEBUG
    assert( baseImg[basePos + k] == UNWRITTEN && 
	    secImg[secPos + k] == UNWRITTEN );
#endif
    baseImg[basePos + k] = secPos + k;
    secImg[secPos + k] = basePos + k;
  }

}
