
#include "align.h"
#include "common.h"
#include "anchor.h"
#include "match.h"
#include "wobble.h"
#include "constraints.h"


#define MIN_MASKING_LEN 8
#define MIN_NOMASKING_LEN 10
#define MIN_WOBBLE_LEN 10

#define MAVID_EXTEND 0

#define LARGE_MIN_LEN 30
#define LARGE_LEN 800000

#define MAX_GAP 2000000

double globalAlign( bp *baseSeq, bp *secSeq, char *baseMask, char *secMask,
		    int baseLen, int secLen,
		    bool noMasking,
		    double **scores, double gap_open, double gap_extend,
		    double **baseProf, double **secProf,
		    int **baseImg, int **secImg,
		    ConstraintsList *cons, WeakList *weak )
{


  *baseImg = new int [baseLen];
  *secImg = new int [secLen];

#ifdef DEBUG
  for( int i = 0; i < baseLen; i++ )
    (*baseImg)[i] = UNWRITTEN;
  for( int i = 0; i < secLen; i++ )
    (*secImg)[i] = UNWRITTEN;
#endif

  int num_anchors = 2;
  double score = 0;

  ConstraintsList *last = NULL;
  for( ConstraintsList *temp = cons; temp != NULL;
       last = temp, temp = temp->next ){
    num_anchors += 2;

    int gap = 0;
    if( last ){
      gap = MAX( temp->cons.first_start - last->cons.first_start - 1,
		 temp->cons.second_start - last->cons.second_start - 1 );
    } else {
      gap = MAX( temp->cons.first_start, temp->cons.second_start );
    }
    num_anchors += gap/MAX_GAP;
  }

  if( last ){
    num_anchors += MAX( baseLen - last->cons.first_start - 1,
			secLen - last->cons.second_start - 1 )/MAX_GAP;
  }
  else {
    num_anchors += MAX( baseLen, secLen )/MAX_GAP;
  }

  int *first_pos = new int [num_anchors];
  int *second_pos = new int [num_anchors];

  first_pos[0] = -1;
  second_pos[0] = -1;
  num_anchors = 1;
  for( ConstraintsList *temp = cons; temp != NULL; temp = temp->next ){
    int first_gap = temp->cons.first_start - first_pos[num_anchors - 1],
      second_gap = temp->cons.second_start - second_pos[num_anchors - 1];
    int gap_count = MAX( first_gap, second_gap )/MAX_GAP;

    for( int i = 0; i < gap_count; i++ ){
      first_pos[num_anchors] = first_pos[num_anchors - 1] + first_gap/(gap_count + 1);
      second_pos[num_anchors++] = second_pos[num_anchors - 1] + second_gap/(gap_count + 1);
    }
    first_pos[num_anchors] = temp->cons.first_start;
    second_pos[num_anchors++] = temp->cons.second_start;
    first_pos[num_anchors] = temp->cons.first_end - 1;
    second_pos[num_anchors++] = temp->cons.second_end - 1;
  }
  
  int first_gap = baseLen - first_pos[num_anchors - 1],
    second_gap = secLen - second_pos[num_anchors - 1];
  int gap_count = MAX( first_gap, second_gap )/MAX_GAP;

  for( int i = 0; i < gap_count; i++ ){
    first_pos[num_anchors] = first_pos[num_anchors - 1] + first_gap/(gap_count + 1);
    second_pos[num_anchors++] = second_pos[num_anchors - 1] + second_gap/(gap_count + 1);
  }

  first_pos[num_anchors] = baseLen;
  second_pos[num_anchors] = secLen;


  // get consistent set
  matchlist *thig = NULL;

  for( int i = 1; i < num_anchors; i++ ){
    matchlist *to_add = new matchlist;

    to_add->basePos = first_pos[i];
    to_add->secPos = second_pos[i];
    to_add->len = 1;
    to_add->next = thig;
    thig = to_add;
  }

  int *temp_first_pos, *temp_second_pos, *temp_lens, num_consistent;
  num_consistent = select_consistent( thig, temp_first_pos, temp_second_pos, temp_lens );

  for( int i = 0; i < num_consistent; i++ ){
    first_pos[i + 1] = temp_first_pos[i];
    second_pos[i + 1] = temp_second_pos[i];
  }
  num_anchors = num_consistent + 1;
  first_pos[num_anchors] = baseLen;
  second_pos[num_anchors] = secLen;

  for( int i = 1; i <= num_anchors; i++ ){
    int *temp_base_img, *temp_sec_img;

    printf( "Aligning [%i,%i] to [%i,%i]\n", first_pos[i - 1] + 1,
	    first_pos[i] - 1, second_pos[i - 1] + 1, second_pos[i] - 1 );

    // correct score when gapping!
    if( first_pos[i] - first_pos[i - 1] - 1 <= 0 ){
            for( int p = second_pos[i - 1] + 1; p <= second_pos[i] - 1; p++ ){
#ifdef DEBUG
	assert( (*secImg)[p] == UNWRITTEN );
#endif
	(*secImg)[p] = -1;
      }
    }
    else if( second_pos[i] - second_pos[i - 1] - 1 <= 0 ){
      for( int p = first_pos[i - 1] + 1; p <= first_pos[i] - 1; p++ ){
#ifdef DEBUG
	assert( (*baseImg)[p] == UNWRITTEN );
#endif
	(*baseImg)[p] = -1;
      }
    }
    else {
      for( WeakList *temp = weak; temp != NULL; temp = temp->next ){
	weak->weak.i -= first_pos[i - 1] + 1;
	weak->weak.j -= second_pos[i - 1] + 1;
      }
      score += globalAlign( baseSeq + first_pos[i - 1] + 1,
			    secSeq + second_pos[i - 1] + 1,
			    baseMask + first_pos[i - 1] + 1,
			    secMask + second_pos[i - 1] + 1,
			    first_pos[i] - first_pos[i - 1] - 1,
			    second_pos[i] - second_pos[i - 1] - 1,
			    noMasking, scores, gap_open, gap_extend,
			    baseProf + first_pos[i - 1] + 1, secProf + second_pos[i - 1] + 1,
			    &temp_base_img, &temp_sec_img, weak );
      for( WeakList *temp = weak; temp != NULL; temp = temp->next ){
	weak->weak.i += first_pos[i - 1] + 1;
	weak->weak.j += second_pos[i - 1] + 1;
      }
      for( int k = 0; k < first_pos[i] - first_pos[i - 1] - 1; k++ ){
#ifdef DEBUG
	assert( first_pos[i - 1] + 1 + k < baseLen &&
		(*baseImg)[first_pos[i - 1] + 1 + k] == UNWRITTEN );
#endif
	(*baseImg)[first_pos[i - 1] + 1 + k] = (temp_base_img[k] == -1) ? -1 :
	  second_pos[i - 1] + 1 + temp_base_img[k];
      }
      for( int k = 0; k < second_pos[i] - second_pos[i - 1] - 1; k++ ){
#ifdef DEBUG
	assert( second_pos[i - 1] + 1 + k < secLen &&
		(*secImg)[second_pos[i - 1] + 1 + k] == UNWRITTEN );
#endif
	(*secImg)[second_pos[i - 1] + 1 + k] = (temp_sec_img[k] == -1) ? -1 :
	  first_pos[i - 1] + 1 + temp_sec_img[k];
      }

      delete[] temp_base_img;
      delete[] temp_sec_img;

    }

    if( i != num_anchors ){
      (*baseImg)[first_pos[i]] = second_pos[i];
      (*secImg)[second_pos[i]] = first_pos[i];
      score += scores[baseSeq[first_pos[i]]][secSeq[second_pos[i]]];
    }
  }

#ifdef DEBUG
  for( int i = 0; i < baseLen; i++ )
    assert( (*baseImg)[i] != UNWRITTEN );
  for( int i = 0; i < secLen; i++ )
    assert( (*secImg)[i] != UNWRITTEN );
#endif

  delete[] first_pos;
  delete[] second_pos;

  return score;

}





// this function takes two sequences along with masks, and fills baseImg and
// secImg with a global alignment of those two sequences. it returns the
// score of that alignment. this function is really just a shell for the
// matchedAlign function which does the real work of the alignment.
double globalAlign( bp *baseSeq, bp *secSeq, char *baseMask, char *secMask,
		    int baseLen, int secLen,
		    bool noMasking,
		    double **scores, double gap_open, double gap_extend,
		    double **baseProf, double **secProf,
		    int **baseImg, int **secImg, WeakList *weak )
{

  *baseImg = new int [baseLen];
  *secImg = new int [secLen];

  // for debugging purposes, initialize all images to be UNWRITTEN
  // and then use this to make sure we never rewrite over part of our alignment
  // and also to check that the entire image gets filled.
#ifdef DEBUG
  for( int i = 0; i < baseLen; i++ )
    (*baseImg)[i] = UNWRITTEN;
  for( int i = 0; i < secLen; i++ )
    (*secImg)[i] = UNWRITTEN;
#endif


  // get the wobble matches
  matchlist *wobbleHead, *wobbleTail;
  bool useWobble = true;

  if( useWobble ){
    getWobbleMatchList( baseSeq, baseMask, baseLen,
			secSeq, secMask, secLen,
			MIN_WOBBLE_LEN,
			wobbleHead, wobbleTail );
  }
  else {
    wobbleHead = wobbleTail = NULL;
  }

  if( baseLen >= LARGE_LEN || secLen >= LARGE_LEN ){
    matchlist *bigMatchesHead, *bigMatchesTail;
    int numAnchors;
    int *baseStarts, *secStarts, *lens;
    double score = 0;
    int *tempBaseImg, *tempSecImg;

    getMatchList( baseSeq, baseMask, baseLen,
		  secSeq, secMask, secLen,
		  LARGE_MIN_LEN, true,
		  bigMatchesHead, bigMatchesTail );

    // if we had any wobble matches, stick them on the front of the match list.
    if( wobbleHead ){
      wobbleTail->next = bigMatchesHead;
      bigMatchesHead = wobbleHead;
    }

    filter_by_weak( bigMatchesHead, weak );
    numAnchors = select_anchors( baseSeq, baseLen, 0, baseLen,
				 secSeq, secLen, 0, secLen,
				 bigMatchesHead, NULL,
				 false,
				 baseStarts, secStarts, lens );

    // if we couldn't get any anchors(or rather, decided not to use them for
    // whatever reason) then do a "basic" alignment of the two regions.
    if( numAnchors == 0 ){
      score = basicAlign( baseSeq, secSeq, baseMask, secMask,
			  0, baseLen, 0, secLen,
			  scores, gap_open, gap_extend,
			  bigMatchesHead, NULL,
			  baseProf, secProf,
			  *baseImg, *secImg );

      return score;
    }

    // now let's extend the anchors and fix them
    int *baseEnds = new int [numAnchors],
      *secEnds = new int [numAnchors];

    for( int i = 0; i < numAnchors; i++ ){
      baseEnds[i] = baseStarts[i] + lens[i] - 1;
      secEnds[i] = secStarts[i] + lens[i] - 1;
    }

    for( int i = 0; i < numAnchors; i++ ){
      int rightExt, leftExt;

      if( i == 0 )
	leftExt = MIN( baseStarts[i], secStarts[i] );
      else
	leftExt = MIN( baseStarts[i] - baseEnds[i - 1] - 1,
		       secStarts[i] - secEnds[i - 1] - 1 );
      leftExt = MIN( leftExt, MAVID_EXTEND );

      if( i == numAnchors - 1 )
	rightExt = MIN( baseLen - baseEnds[i] - 1, secLen - secEnds[i] - 1 );
      else
	rightExt = MIN( baseStarts[i + 1] - baseEnds[i] - 1,
			secStarts[i + 1] - secEnds[i] - 1 );
      rightExt = MIN( rightExt, MAVID_EXTEND );

      // fix the anchor
      fixAnchor( baseStarts[i], secStarts[i], lens[i],
		 *baseImg, *secImg );
      for( int p = 0; p < lens[i]; p++ )
	score += scores[baseSeq[baseStarts[i] + p]][secSeq[secStarts[i] + p]];


      // now align the extensions and adjust the Starts and Ends
      // accordingly.
      baseStarts[i] -= leftExt; secStarts[i] -= leftExt;
      baseEnds[i] += rightExt; secEnds[i] += rightExt;
      gappedAlign( baseSeq, secSeq, baseStarts[i], baseStarts[i] + leftExt - 1,
		   secStarts[i], secStarts[i] + leftExt - 1,
		   scores, gap_open, gap_extend,
		   NULL, NULL,
		   *baseImg, *secImg );
      gappedAlign( baseSeq, secSeq, baseEnds[i] - rightExt + 1, baseEnds[i],
		   secEnds[i] - rightExt + 1, secEnds[i],
		   scores, gap_open, gap_extend,
		   NULL, NULL,
		   *baseImg, *secImg );

      // the ends of the extensions may be gapped in which case we want to
      // unalign them
      for( int j = baseStarts[i]; (*baseImg)[j] == GAPPED; j++ ){
#ifdef DEBUG
	(*baseImg)[j] = UNWRITTEN;
#endif
	baseStarts[i]++;
      }
      for( int j = baseEnds[i]; (*baseImg)[j] == GAPPED; j-- ){
#ifdef DEBUG
	(*baseImg)[j] = UNWRITTEN;
#endif
	baseEnds[i]--;
      }
      for( int j = secStarts[i]; (*secImg)[j] == GAPPED; j++ ){
#ifdef DEBUG
	(*secImg)[j] = UNWRITTEN;
#endif
	secStarts[i]++;
      }
      for( int j = secEnds[i]; (*secImg)[j] == GAPPED; j-- ){
#ifdef DEBUG
	(*secImg)[j] = UNWRITTEN;
#endif
	secEnds[i]--;
      }

    }

#ifdef DEBUG
    cout << "anchoring from " << 0 << " (" << 0 << ")  to  "
	 << baseLen << " (" << secLen << ")\n";
    for( int i = 0; i < numAnchors; i++ ){
      cout << "anchor: (" << baseStarts[i] << ", " << baseEnds[i] << ") x ("
	   << secStarts[i] << ", " << secEnds[i] << ")\t";
      cout << "length " << baseEnds[i] - baseStarts[i] + 1 << " x "
	   << secEnds[i] - secStarts[i] + 1 << endl;
    }
    cout << endl;
#endif
    

    // now start recursing into the gaps.
    for( int i = 0; i < numAnchors; i++ ){
      // now recurse into the region in-between the anchors
      int baseOffset, secOffset;
      int *tempBaseImg, *tempSecImg;

      if( i == 0 ){
	baseOffset = 0;
	secOffset = 0;
      }
      else {
	baseOffset = baseEnds[i - 1] + 1;
	secOffset = secEnds[i - 1] + 1;
      }

      score += globalAlign( baseSeq + baseOffset, secSeq + secOffset,
			    baseMask + baseOffset, secMask + secOffset,
			    baseStarts[i] - baseOffset,
			    secStarts[i] - secOffset,
			    noMasking,
			    scores, gap_open, gap_extend,
			    baseProf + baseOffset, secProf + secOffset,
			    &tempBaseImg, &tempSecImg, NULL );
      for( int k = 0; k < baseStarts[i] - baseOffset; k++ ){
#ifdef DEBUG
	assert( (*baseImg)[baseOffset + k] == UNWRITTEN );
#endif
	(*baseImg)[baseOffset + k] = (tempBaseImg[k] == -1) ? -1 : secOffset + tempBaseImg[k];
      }
      for( int k = 0; k < secStarts[i] - secOffset; k++ ){
#ifdef DEBUG
	assert( (*secImg)[secOffset + k] == UNWRITTEN );
#endif
	(*secImg)[secOffset + k] = (tempSecImg[k] == -1) ? -1 : baseOffset + tempSecImg[k];
      }
      delete[] tempBaseImg;
      delete[] tempSecImg;
    }

    // and finally recurse into the region after the last anchor
    int baseOffset = baseEnds[numAnchors - 1] + 1;
    int secOffset = secEnds[numAnchors - 1] + 1;


    score += globalAlign( baseSeq + baseOffset, secSeq + secOffset,
			  baseMask + baseOffset, secMask + secOffset,
			  baseLen - baseOffset, secLen - secOffset,
			  noMasking,
			  scores, gap_open, gap_extend,
			  baseProf + baseOffset, secProf + secOffset,
			  &tempBaseImg, &tempSecImg, NULL );
    for( int i = 0; i < baseLen - baseOffset; i++ ){
#ifdef DEBUG
      assert( (*baseImg)[baseOffset + i] == UNWRITTEN );
#endif
      (*baseImg)[baseOffset + i] = (tempBaseImg[i] == -1) ? -1 : secOffset + tempBaseImg[i];
    }
    for( int i = 0; i < secLen - secOffset; i++ ){
#ifdef DEBUG
      assert( (*secImg)[secOffset + i] == UNWRITTEN );
#endif
      (*secImg)[secOffset + i] = (tempSecImg[i] == -1) ? -1 : baseOffset + tempSecImg[i];
    }

    delete[] tempBaseImg;
    delete[] tempSecImg;
    delete[] baseStarts;
    delete[] baseEnds;
    delete[] secStarts;
    delete[] secEnds;
    delete[] lens;

    for( matchlist *temp = bigMatchesHead; temp != NULL; ){
      matchlist *next = temp->next;
      delete temp;
      temp = next;
    }
    return score;

  }

  // we get the list of matches both with and without masking and append
  // the unmasked list to the end of the masking on.
  matchlist *maskingHead, *maskingTail, *nomaskingHead, *nomaskingTail;

  getMatchList( baseSeq, baseMask, baseLen,
		secSeq, secMask, secLen,
		MIN_MASKING_LEN, true,
		maskingHead, maskingTail );
  if( !noMasking ){
    getMatchList( baseSeq, baseMask, baseLen,
		  secSeq, secMask, secLen,
		  MIN_NOMASKING_LEN, false,
		  nomaskingHead, nomaskingTail );
  }
  else
    nomaskingHead = nomaskingTail = NULL;


  /*while( wobbleHead != NULL ){
    cout << wobbleHead->basePos << " " << wobbleHead->secPos << " " << wobbleHead->len << endl;
    wobbleHead = wobbleHead->next;
  }
  cout << endl << endl;

  while( maskingHead != NULL ){
    cout << maskingHead->basePos << " " << maskingHead->secPos << " " << maskingHead->len << endl;
    maskingHead = maskingHead->next;
  }
  cout << endl << endl;
  exit(0);*/

  // if we had any wobble matches, stick them on the front of the match list.
  if( wobbleHead ){
    wobbleTail->next = maskingHead;
    maskingHead = wobbleHead;
  }

  // and then have matchedAlign do the real work.
  filter_by_weak( maskingHead, weak );
  filter_by_weak( nomaskingHead, weak );
  return matchedAlign( baseSeq, secSeq, baseMask, secMask,
		       baseLen, secLen,
		       0, baseLen, 0, secLen,
		       maskingHead, nomaskingHead,
		       scores, gap_open, gap_extend,
		       baseProf, secProf,
		       *baseImg, *secImg );

}


// this takes two sequences along with a linked list of matches between them
// and computes a global alignment of the regions [baseFrom, baseTo) and
// [secFrom, secTo), filling the corresponding regions of baseImg and secImg
// with that alignment. from the list of matches, it chooses a set of "anchors"
// and then calls itself recursively on the regions between the anchors. if
// anchors cannot be used for whatever reason then it calls the basicAlign
// function to complete the alignment.
double matchedAlign( bp *baseSeq, bp *secSeq, char *baseMask, char *secMask,
		     int baseLen, int secLen,
		     int baseFrom, int baseTo, int secFrom, int secTo,
		     matchlist *maskingMatches, matchlist *nomaskingMatches,
		     double **scores, double gap_open, double gap_extend,
		     double **baseProf, double **secProf,
		     int *baseImg, int *secImg )
{

  // Call the select_anchors function to get the anchor locations and lengths.
  double score = 0;
  int numAnchors;
  int *baseLocs, *secLocs, *lens;

  bool isShort = ((baseTo - baseFrom + 2) <= 2048)
    && ((secTo - secFrom + 2) <= 2048);
  numAnchors = select_anchors( baseSeq, baseLen, baseFrom, baseTo,
			       secSeq, secLen, secFrom, secTo,
			       maskingMatches, nomaskingMatches,
			       isShort,
			       baseLocs, secLocs, lens );

  // if we couldn't get any anchors(or rather, decided not to use them for
  // whatever reason) then do a "basic" alignment of the two regions.
  if( numAnchors == 0 ){
    score = basicAlign( baseSeq, secSeq, baseMask, secMask,
			baseFrom, baseTo, secFrom, secTo,
			scores, gap_open, gap_extend,
			maskingMatches, nomaskingMatches,
			baseProf, secProf,
			baseImg, secImg );

    return score;
  }

  // now let's extend the anchors and fix them
  int *baseEnds = new int [numAnchors],
    *secEnds = new int [numAnchors];

  for( int i = 0; i < numAnchors; i++ ){
    baseEnds[i] = baseLocs[i] + lens[i] - 1;
    secEnds[i] = secLocs[i] + lens[i] - 1;
  }

  for( int i = 0; i < numAnchors; i++ ){
    int rightExt, leftExt;

    if( i == 0 )
      leftExt = MIN( baseLocs[i] - baseFrom, secLocs[i] - secFrom );
    else
      leftExt = MIN( baseLocs[i] - baseEnds[i - 1] - 1,
		     secLocs[i] - secEnds[i - 1] - 1 );
    leftExt = MIN( leftExt, MAVID_EXTEND );

    if( i == numAnchors - 1 )
      rightExt = MIN( baseTo - baseEnds[i] - 1, secTo - secEnds[i] - 1 );
    else
      rightExt = MIN( baseLocs[i + 1] - baseEnds[i] - 1,
		      secLocs[i + 1] - secEnds[i] - 1 );
    rightExt = MIN( rightExt, MAVID_EXTEND );

    // fix the anchor
    fixAnchor( baseLocs[i], secLocs[i], lens[i],
	       baseImg, secImg );
    for( int p = 0; p < lens[i]; p++ )
      score += scores[baseSeq[baseLocs[i] + p]][secSeq[secLocs[i] + p]];


    // now align the extensions and adjust the Locs and Ends
    // accordingly.
    baseLocs[i] -= leftExt; secLocs[i] -= leftExt;
    baseEnds[i] += rightExt; secEnds[i] += rightExt;
    gappedAlign( baseSeq, secSeq, baseLocs[i], baseLocs[i] + leftExt - 1,
		 secLocs[i], secLocs[i] + leftExt - 1,
		 scores, gap_open, gap_extend,
		 NULL, NULL,
		 baseImg, secImg );
    gappedAlign( baseSeq, secSeq, baseEnds[i] - rightExt + 1, baseEnds[i],
		 secEnds[i] - rightExt + 1, secEnds[i],
		 scores, gap_open, gap_extend,
		 NULL, NULL,
		 baseImg, secImg );

    // the ends of the extensions may be gapped in which case we want to
    // unalign them
    for( int j = baseLocs[i]; baseImg[j] == GAPPED; j++ ){
#ifdef DEBUG
      baseImg[j] = UNWRITTEN;
#endif
      baseLocs[i]++;
    }
    for( int j = baseEnds[i]; baseImg[j] == GAPPED; j-- ){
#ifdef DEBUG
      baseImg[j] = UNWRITTEN;
#endif
      baseEnds[i]--;
    }
    for( int j = secLocs[i]; secImg[j] == GAPPED; j++ ){
#ifdef DEBUG
      secImg[j] = UNWRITTEN;
#endif
      secLocs[i]++;
    }
    for( int j = secEnds[i]; secImg[j] == GAPPED; j-- ){
#ifdef DEBUG
      secImg[j] = UNWRITTEN;
#endif
      secEnds[i]--;
    }

  }


#ifdef DEBUG
  cout << "anchoring from " << baseFrom << " (" << secFrom << ")  to  " << baseTo << " (" << secTo << ")\n";
  for( int i = 0; i < numAnchors; i++ ){
    cout << "anchor: (" << baseLocs[i] << ", " << baseEnds[i] << ") x ("
	 << secLocs[i] << ", " << secEnds[i] << ")\t";
    cout << "length " << baseEnds[i] - baseLocs[i] + 1 << " x "
	 << secEnds[i] - secLocs[i] + 1 << "\t";
    for( int p = baseLocs[i]; p <= baseEnds[i]; p++ )
      cout << DNA2alpha(baseSeq[p]);
    cout << endl;
  }
  cout << endl;
#endif
    

  // once the anchors have been extended, we can filter the matches into
  // numAnchors + 1 lists according to which anchors they're between.
  matchlist **maskingMatcheses, **nomaskingMatcheses;

  maskingMatcheses = filterMatches( maskingMatches, baseLocs, baseEnds,
				    secLocs, secEnds,
				    numAnchors );
  nomaskingMatcheses = filterMatches( nomaskingMatches, baseLocs, baseEnds,
				      secLocs, secEnds,
				      numAnchors );


  // now start recursing into the gaps.
  for( int i = 0; i < numAnchors; i++ ){
    // now recurse into the region in-between the anchors
    int baseOffset, secOffset;

    if( i == 0 ){
      baseOffset = baseFrom;
      secOffset = secFrom;
    }
    else {
      baseOffset = baseEnds[i - 1] + 1;
      secOffset = secEnds[i - 1] + 1;
    }

    score += matchedAlign( baseSeq, secSeq, baseMask, secMask,
			   baseLen, secLen,
			   baseOffset, baseLocs[i], secOffset, secLocs[i],
			   maskingMatcheses[i], nomaskingMatcheses[i],
			   scores, gap_open, gap_extend,
			   baseProf, secProf,
			   baseImg, secImg );

  }

  // and finally recurse into the region after the last anchor
  int baseOffset = baseEnds[numAnchors - 1] + 1;
  int secOffset = secEnds[numAnchors - 1] + 1;

  score += matchedAlign( baseSeq, secSeq, baseMask, secMask,
			 baseLen, secLen,
			 baseOffset, baseTo, secOffset, secTo,
			 maskingMatcheses[numAnchors],
			 nomaskingMatcheses[numAnchors],
			 scores, gap_open, gap_extend,
			 baseProf, secProf,
			 baseImg, secImg );

  delete[] baseLocs;
  delete[] secLocs;
  delete[] baseEnds;
  delete[] secEnds;
  delete[] lens;
  delete[] maskingMatcheses;
  delete[] nomaskingMatcheses;


#ifdef DEBUG
  // verify output
  assert( verifyAlignment( baseFrom, baseTo, secFrom, secTo,
			   baseImg, secImg ) );
#endif

  return score;

}


// this function takes two sequences and computes a global alignment of the
// regions [baseFrom, baseTo) and [secFrom, secTo). if the regions are "short"
// then it uses a Smith-Waterman procedure to compute an optimal(with respect
// to the values match, mismatch, gapOpen and gap) alignment. otherwise, it
// computes a decidely non-optimal alignment(the reasoning being that if we
// have reached this point with "long" regions then those regions are very
// non-homologous so the optimality of an "optimal alignment" is meaningless
// anyways.
double basicAlign( bp *baseSeq, bp *secSeq, char *baseMask, char *secMask,
		   int baseFrom, int baseTo, int secFrom, int secTo,
		   double **scores, double gap_open, double gap_extend,
		   matchlist *maskingMatches, matchlist *nomaskingMatches,
		   double **baseProf, double **secProf,
		   int *baseImg, int *secImg )
{

#ifdef DEBUG
  // verify input
  assert( verifyUnwritten( baseFrom, baseTo, baseImg ) );
  assert( verifyUnwritten( secFrom, secTo, secImg ) );
#endif

  double score;

  // right now we just delete the matches passed to this function. one day
  // maybe we'll do something cooler.
  while( maskingMatches != NULL ){
    matchlist *next = maskingMatches->next;
    delete maskingMatches;
    maskingMatches = next;
  }
  while( nomaskingMatches != NULL ){
    matchlist *next = nomaskingMatches->next;
    delete nomaskingMatches;
    nomaskingMatches = next;
  }

  // if one of our regions is empty then just gap the other region
  if( baseFrom == baseTo ){
    for( int i = secFrom; i < secTo; i++ ){
      secImg[i] = -1;
    }
    return gap_open + gap_extend * (secTo - secFrom - 1);
  }

  if( secFrom == secTo ){
    for( int i = baseFrom; i < baseTo; i++ ){
      baseImg[i] = -1;
    }
    return gap_open + gap_extend * (baseTo - baseFrom - 1);
  }

  // if the regions are "short" but non-empty then use gappedAlign
  if( baseFrom != baseTo && secFrom != secTo &&
      (baseTo - baseFrom + 1 < 30000) && (secTo - secFrom + 1 < 30000) &&
      ((baseTo - baseFrom + 1) * (secTo - secFrom + 1) <= 4048 * 4048) ){
    score = gappedAlign( baseSeq, secSeq, baseFrom, baseTo - 1,
			 secFrom, secTo - 1,
			 scores, gap_open, gap_extend,
			 baseProf, secProf,
			 baseImg, secImg );
  }
  else {
    // otherwise let's do this garbage alignment.
    for( int i = baseFrom; i < baseTo; i++ ){
      baseImg[i] = GAPPED;
    }
    for( int i = secFrom; i < secTo; i++ ){
      secImg[i] = GAPPED;
    }

    if( baseFrom == baseTo ){
      score = gap_open + (secTo - secFrom - 1) * gap_extend;
    }
    else if( secFrom == secTo ){
      score = gap_open + (baseTo - baseFrom - 1) * gap_extend;
    }
    else {
      score = 2 * gap_open + (secTo + baseTo - secFrom - baseFrom - 2) * gap_extend;
    }
  }

#ifdef DEBUG
  // verify output
  assert( verifyAlignment( baseFrom, baseTo, secFrom, secTo,
			   baseImg, secImg ) );
#endif  

  return score;

}



#define BIG_INT 0x3FFFFFFF

double scoreMat[][5] = {{4.0,-4.0,-4.0,-4.0,-4.0},{-4.0,4.0,-4.0,-4.0,-4.0},{-4.0,-4.0,4.0,-4.0,-4.0},
		    {-4.0,-4.0,-4.0,4.0,-4.0}, {-4.0,-4.0,-4.0,-4.0,-4.0}};

double gappedAlign( bp *first, bp *second, int fromfirst, int tofirst,
		    int fromsecond, int tosecond,
		    double **scores, double gap_open, double gap_extend,
		    double **baseProf, double **secProf,
		    int *imagefirst, int *imagesecond )
{
  // Performs a global alignment on two sequences. M contains the scores,
  //    P contains the pointers. The encoding is 0=up, 1=left, 2=upleft.
  //    The input is two sequences, and positions within them to align.
  //    The answer is returned in 'imagefirst' and 'imagesecond '(zero based).

  register int i = 0, j = 0;

  enum alignPtrs {toD0, toD1, toD2, toU0, toU1, toU2, toL0, toL1, toL2};

  int lengthfirst = tofirst - fromfirst + 1,
    lengthsecond = tosecond - fromsecond + 1;

  double score = 0;

  if( lengthfirst == 0 && lengthsecond == 0 ){
    return score;
  }

  double *thisCol0, *lastCol0, *thisCol1, *lastCol1, *thisCol2, *lastCol2;
  thisCol0 = new double [lengthfirst + 1];
  thisCol1 = new double [lengthfirst + 1];
  thisCol2 = new double [lengthfirst + 1];
  lastCol0 = new double [lengthfirst + 1];
  lastCol1 = new double [lengthfirst + 1];
  lastCol2 = new double [lengthfirst + 1];

  char *MatrixP = new char [3 * (lengthfirst + 1) * (lengthsecond + 1)];

  // Boundary Conditions:

  lastCol0[0] = 0;
  lastCol1[0] = -BIG_INT;
  lastCol2[0] = -BIG_INT;

  int pPtr = 3;

  double gapScore = gap_open;
  for( i = 1; i < lengthfirst + 1; i++, pPtr += 3 ){
    lastCol0[i] = -BIG_INT;
    lastCol1[i] = gapScore;
    lastCol2[i] = gapScore + gap_open;
    MatrixP[pPtr] = toU1;
    MatrixP[pPtr + 1] = toU1;
    MatrixP[pPtr + 2] = toU1;
    gapScore += gap_extend;
  }
  MatrixP[0] = 0;
  MatrixP[1] = toU0;
  MatrixP[2] = toU0;

  pPtr = 3 * (lengthfirst + 1);
  for( j = 1; j < lengthsecond + 1; j++, pPtr += 3*(lengthfirst + 1) ){
    MatrixP[pPtr] = toL2;
    MatrixP[pPtr + 1] = toL2;
    MatrixP[pPtr + 2] = toL2;
  }
  MatrixP[3*(lengthfirst + 1)] = toL0;

// End of boundary conditions.

  first += fromfirst;
  second += fromsecond;
  if( baseProf != NULL )
    baseProf += fromfirst;
  if( secProf != NULL )
    secProf += fromsecond;

  int i_1, j_1; // these will be i-1 and j-1 in the next loop.
  pPtr = 3 * (lengthfirst + 1) + 3;

  thisCol0[0] = -BIG_INT;
  thisCol1[0] = 2 * gap_open;
  thisCol2[0] = gap_open;
  double score_col[5];

  for( j = 1, j_1 = 0; j < lengthsecond + 1; j++, j_1++ ){
    char c = second[j_1];

    if( secProf ){
      for( int k = 0; k < 4; k++ ){
	// FIX! ADD BACK IN GAP PENALTY! MAYBE!
	score_col[k] = secProf[j_1][DNA_A]*scores[k][DNA_A] + secProf[j_1][DNA_C]*scores[k][DNA_C] +
	  secProf[j_1][DNA_G]*scores[k][DNA_G] + secProf[j_1][DNA_T]*scores[k][DNA_T];
	  //+ secProf[j_1][4]*scores[k][DNA_GAP];
      }
      score_col[4] = secProf[j_1][DNA_A]*scores[DNA_GAP][DNA_A] + secProf[j_1][DNA_C]*scores[DNA_GAP][DNA_C] +
	secProf[j_1][DNA_G]*scores[DNA_GAP][DNA_G] + secProf[j_1][DNA_T]*scores[DNA_GAP][DNA_T];
    }
    for( i = 1, i_1 = 0; i < lengthfirst + 1; i++, i_1++, pPtr += 3 ){
      double U0 = thisCol0[i_1],
	U1 = thisCol1[i_1],
	U2 = thisCol2[i_1],
	L0 = lastCol0[i],
	L1 = lastCol1[i],
	L2 = lastCol2[i],
	D0 = lastCol0[i_1],
	D1 = lastCol1[i_1],
	D2 = lastCol2[i_1];

      double d;
      if( baseProf ){
	d = score_col[DNA_A]*baseProf[i_1][DNA_A] + score_col[DNA_C]*baseProf[i_1][DNA_C]
	  + score_col[DNA_G]*baseProf[i_1][DNA_G] + score_col[DNA_T]*baseProf[i_1][DNA_T]; // + score_col[4]*baseProf[i_1][4];
      }
      else {
	d = scores[first[i_1]][c];
      }

      if( D0 >= D1 ){
	if( D2 >= D0 ){
	  thisCol0[i] = D2 + d;
	  MatrixP[pPtr] = toD2;
	}
	else {
	  thisCol0[i] = D0 + d;
	  MatrixP[pPtr] = toD0;
	}
      }
      else {
	if( D2 >= D1 ){
	  thisCol0[i] = D2 + d;
	  MatrixP[pPtr] = toD2;
	}
	else {
	  thisCol0[i] = D1 + d;
	  MatrixP[pPtr] = toD1;
	}
      }

      if( U0 >= U2 ){
	if( U0 + gap_open >= U1 + gap_extend ){
	  thisCol1[i] = U0 + gap_open;
	  MatrixP[pPtr + 1] = toU0;
	}
	else {
	  thisCol1[i] = U1 + gap_extend;
	  MatrixP[pPtr + 1] = toU1;
	}
      }
      else {
	if( U2 + gap_open >= U1 + gap_extend ){
	  thisCol1[i] = U2 + gap_open;
	  MatrixP[pPtr + 1] = toU2;
	}
	else {
	  thisCol1[i] = U1 + gap_extend;
	  MatrixP[pPtr + 1] = toU1;
	}
      }

      if( L0 >= L1 ){
	if( L0 + gap_open >= L2 + gap_extend ){
	  thisCol2[i] = L0 + gap_open;
	  MatrixP[pPtr + 2] = toL0;
	}
	else {
	  thisCol2[i] = L2 + gap_extend;
	  MatrixP[pPtr + 2] = toL2;
	}
      }
      else {
	if( L1 + gap_open >= L2 + gap_extend ){
	  thisCol2[i] = L1 + gap_open;
	  MatrixP[pPtr + 2] = toL1;
	}
	else {
	  thisCol2[i] = L2 + gap_extend;
	  MatrixP[pPtr + 2] = toL2;
	}
      }
    }
    pPtr += 3;
    swap( thisCol0, lastCol0 );
    swap( thisCol1, lastCol1 );
    swap( thisCol2, lastCol2 );
    thisCol0[0] = -BIG_INT;
    thisCol1[0] = lastCol1[0] + gap_extend;
    thisCol2[0] = lastCol2[0] + gap_extend;
  }

  // Find the optimal alignment and fill in the arrays 'imagefirst'
  //    and 'imagesecond'.

  i = lengthfirst;
  j = lengthsecond;
  int l;

  if( lastCol0[i] >= lastCol1[i] ){
    if( lastCol0[i] >= lastCol2[i] ){
      l = 0;
    }
    else {
      l = 2;
    }
  }
  else {
    if( lastCol1[i] >= lastCol2[i] ){
      l = 1;
    }
    else {
      l = 2;
    }
  }

  if( l == 0 )
    score = lastCol0[i];
  else if( l == 1 )
    score = lastCol1[i];
  else
    score = lastCol2[i];

  while( i > 0 || j > 0 ){
    switch( MatrixP[3*((lengthfirst + 1)*j + i) + l] ){
    case toU0:
#ifdef DEBUG
      assert( imagefirst[fromfirst + i - 1] == UNWRITTEN );
#endif
      imagefirst[fromfirst + --i] = -1;
      l = 0;
      break;
    case toU1:
#ifdef DEBUG
      assert( imagefirst[fromfirst + i - 1] == UNWRITTEN );
#endif
      imagefirst[fromfirst + --i] = -1;
      l = 1;
      break;
    case toU2:
#ifdef DEBUG
      assert( imagefirst[fromfirst + i - 1] == UNWRITTEN );
#endif
      imagefirst[fromfirst + --i] = -1;
      l = 2;
      break;
    case toL0:
#ifdef DEBUG
      assert( imagesecond[fromsecond + j - 1] == UNWRITTEN );
#endif
      imagesecond[fromsecond + --j] = -1;
      l = 0;
      break;
    case toL1:
#ifdef DEBUG
      assert( imagesecond[fromsecond + j - 1] == UNWRITTEN );
#endif
      imagesecond[fromsecond + --j] = -1;
      l = 1;
      break;
    case toL2:
#ifdef DEBUG
      assert( imagesecond[fromsecond + j - 1] == UNWRITTEN );
#endif
      imagesecond[fromsecond + --j] = -1;
      l = 2;
      break;
    case toD0:
#ifdef DEBUG
      assert( imagefirst[fromfirst + i - 1] == UNWRITTEN );
      assert( imagesecond[fromsecond + j - 1] == UNWRITTEN );
#endif
      imagefirst[fromfirst + i - 1] = fromsecond + j - 1;
      imagesecond[fromsecond + --j] = fromfirst + --i;
      l = 0;
      break;
    case toD1:
#ifdef DEBUG
      assert( imagefirst[fromfirst + i - 1] == UNWRITTEN );
      assert( imagesecond[fromsecond + j - 1] == UNWRITTEN );
#endif
      imagefirst[fromfirst + i - 1] = fromsecond + j - 1;
      imagesecond[fromsecond + --j] = fromfirst + --i;
      l = 1;
      break;
    case toD2:
#ifdef DEBUG
      assert( imagefirst[fromfirst + i - 1] == UNWRITTEN );
      assert( imagesecond[fromsecond + j - 1] == UNWRITTEN );
#endif
      imagefirst[fromfirst + i - 1] = fromsecond + j - 1;
      imagesecond[fromsecond + --j] = fromfirst + --i;
      l = 2;
      break;
#if DEBUG
    default:  assert(0);
#endif
    }
  }

  delete[] thisCol0;
  delete[] thisCol1;
  delete[] thisCol2;
  delete[] lastCol0;
  delete[] lastCol1;
  delete[] lastCol2;
  delete[] MatrixP;
  
  return score;

}


bool verifyAlignment( int baseFrom, int baseTo, int secFrom, int secTo,
		      int *baseImg, int *secImg )
{

  for( int i = baseFrom; i < baseTo; i++ ){
    if( baseImg[i] != GAPPED ){
      if( secFrom > baseImg[i] || baseImg[i] >= secTo )
	return false;
      if( i != 0 && baseImg[i - 1] >= baseImg[i] )
	return false;
    }
  }
  for( int i = secFrom; i < secTo; i++ ){
    if( secImg[i] != GAPPED ){
      if( baseFrom > secImg[i] || secImg[i] >= baseTo )
	return false;
      if( i != 0 && secImg[i - 1] >= secImg[i] )
	return false;
    }
  }

  return true;

}


bool verifyUnwritten( int from, int to, int *img )
{

  for( int i = from; i < to; i++ )
    if( img[i] != UNWRITTEN )
      return false;

  return true;

}
