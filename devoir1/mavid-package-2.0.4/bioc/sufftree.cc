
#include "align.h"
#include "sufftree.h"
#define COPYHUGETUPLEENDS

Node * SuffixTree::addNode(int t,int b,Node * p) {
  myNodes[nNodes].children=myNodeLinks+(nNodes*ALPHASIZE);

  myNodes[nNodes].setParams(t,b,p,nNodes);

  return(&myNodes[nNodes++]);
}

/* McCreight's Algorithm for suffix tree construction
 * Note that the variable names mostly follow the terminology in the original
 * paper:
 * [McCreight76] E. M. McCreight. A Space-Economical Suffix Tree Construction
 *               Algorithm. _Journal_of_the_ACM, 23(2):262-272, 1976.
 */
SuffixTree::SuffixTree(const char *strg, const char *alphabet) {
  alpha=alphabet;
  buildAlphaHash();

  s=strg;
  sLen=strlen(s);
  if(sLen<2) {
    cerr << "ERROR: Suffix tree corner case [string length < 2]\n";
    exit(1);
  };

  // Hog up all space we may ever need
  // NOTE: The storage requirement for this is upper bounded
  //       by 2 x actual usage, and in reality probably much
  //       smaller than the usage from any solution based
  //       on dynamic node allocation. This is a GOOD thing =)
  myNodes     = new Node[2*sLen-1]();
  myNodeLinks = new Node * [(2*sLen-1)*ALPHASIZE];
  nNodes=0;

  root=addNode(0,EMPTY,NULL);
  SUFFDEBUG { cout << "MAKE root: " << root << "\n"; };

  Node * newnode = addNode(0,sLen-1,root);
  SUFFDEBUG { cout << "MAKE first: " << newnode << "\n"; };
  root->children[alphaIndex(s[0])]=newnode;

  head_bot=EMPTY;
  head_prev_cl=root; // contracted locus
  // head_prev_el=root; // extended locus
  head_prev=root;

  // first char indices of the head segments, *ON* the head of the
  // current suffix
  int chi,alpha,beta,gamma;
  Node *c, *d, *f, *fprev;

  for(int i=1;i<sLen;i++) {
    // insert the suffix starting with s[i]

    // Substep A
    chi=i-1;
    if(head_bot==EMPTY) {
      alpha=chi;
      beta=chi;
      gamma=i;
    } else {
      alpha=i;
      beta=i-1+head_prev_cl->total;
      if(beta<alpha) beta=alpha; // this is iff alpha is empty
      gamma=head_bot+1;
    };
    // found chi/alpha/beta

    if(alpha==beta) {
       c=root;
    } else {
      // assert(alpha<beta);
      c=head_prev->slink;
    };

    if(!c) c=root;

    SUFFDEBUG {
      cout << "i: " << i << ", alpha: " << alpha << ", beta: " << beta
	   << ", gamma: " << gamma << " / c: " << c << endl;
    };
    

    // Substep B: Rescanning
    int previdx;
    fprev=c;
    previdx=alphaIndex(strg[c->total+i]);
    f=fprev->children[previdx];

    // REMOVE!!!
    int count = 0;

    while(f && f->total<gamma-i) {
      fprev=f;
      previdx=alphaIndex(strg[f->total+i]);
      f=f->children[previdx];
      count++;
    };

    // Assertion:
    // At this point, EITHER beta is empty, so fprev=c. Then either f=NULL
    // if we fall out of the tree already or f is whatever node will be the
    // the starting point for scanning for gamma. OR beta is not empty in
    // which case, since alpha+beta is in the tree already, its extended
    // locus is f (so alpha+beta ends between f->top and f->bot,
    // inclusive), and f's parent is fprev.
    
    // now insert alpha+beta node if necessary between fprev and f
    // int gamma_is_empty=0;
    int found_old_d=0;
    int newnodereqd=1; // for alpha+beta+gamma

    Node *slink_to;
    if(/*alpha!=beta &&*/ f && f->total>gamma-i && fprev->total<gamma-i) {
      // (so it's necessary)
      d=addNode(i+fprev->total,gamma-1,fprev);
      SUFFDEBUG { cout << "MAKE d: " << d << " / fprev: "
		   << fprev << " / f: " << f << "\n"; };

      f->top+=(gamma-i-(fprev->total));
      d->children[alphaIndex(strg[f->top])]=f;
      previdx=alphaIndex(strg[fprev->total+i]);
      fprev->children[previdx]=d;
      head_prev_cl=fprev;
      // gamma_is_empty=1;
    } else if(f && f->total==gamma-i) {
      // or we're there already and ready to start scanning
      d=f;
      found_old_d=1;
      head_prev_cl=f;
    } else if(fprev->total==gamma-i) {
      // we fell out of the tree at fprev.
      // gamma's gotta be empty
      // OR we're at root
      d=fprev;
      head_prev_cl=fprev;
      // gamma_is_empty=1;
    } else {
      cerr << "Can't figure out locus of alpha+beta at i=" << i << ".\n";
      exit(1);
    };

    // d is now definitely the locus of alpha+beta.

    // Substep C: Scanning
    slink_to=d; // (found_old_d?fprev:d);
    if(head_prev && head_prev->total==1 && !head_prev->slink) {
      head_prev->slink=root;
    } else if(head_prev && head_prev!=root && !head_prev->slink) {
      SUFFDEBUG {
	if(slink_to && slink_to->total+1!=head_prev->total) {
	  cerr << "Bzzt. Something's screwy.\n";
	  exit(1);
	};
      };
      head_prev->slink=slink_to;
    };
    
    int j,k;

    f=d->children[alphaIndex(strg[gamma])];

    if(f) {
      // if we *might* actually be inserting a new node *between* two nodes

      // if(!gamma_is_empty) {
      // if gamma is not empty (ie alpha+beta locus existed already), we
      // need to scan down to make a new locus for alpha+beta+gamma

      // fprev=(found_old_d?fprev:((gamma==i)?root:d));
      // d=fprev->children[alphaIndex(strg[gamma])];

      k = f->top; // d->top + gamma - (fprev->total+i);

      for(j=gamma;j<sLen;j++) {
	// this loop will *ALWAYS* terminate by "break", since the
	// search has to "fall out of the tree
      
	if(strg[k]!=strg[j]) break; // the tail diverged in mid-edge

	// advance down the tree
	if(k==f->bot) {
	  previdx=alphaIndex(strg[j+1]);
	  newnode=f->children[previdx];
	  if(newnode==NULL) {
	    newnodereqd=0;
	    j++;
	    break; // the tail diverged at a node
	  };
	  d=f;
	  f=newnode;
	  k=f->top;
	} else {
	  k++;
	};
      };

      // at this point, original meaning of "d" is corrupted:
      // now, if alpha+beta+gamma has a locus, it's f
      //      if not, it has to be inserted between d and f

      // shouldn't need this:
      // if(k==0) newnodereqd=0;

      Node *attach_to=d; // ((d==f)?fprev:d);
      if(newnodereqd) {
	// make a new node if the tail diverged in mid-edge
	// assert(fprev!=NULL);
	previdx=alphaIndex(strg[j-k+f->top]); // redundant code(?)
	newnode=addNode(j-k+f->top,j-1,attach_to);
	SUFFDEBUG { cout << "MAKE newnodereqd: " << newnode << " / j: "
		     << j << ", k: " << k << "\n"; };

	f->top=k;
	attach_to->children[previdx]=newnode;

	newnode->children[alphaIndex(strg[k])]=f;
	previdx=alphaIndex(strg[j]);
      } else {
	newnode=f;
      };
      /* } else {
	 // so gamma is empty, we only need to build a terminal node
	 j=gamma;
	 newnode=d;
	 previdx=alphaIndex(strg[gamma]);
	 };
      */
    } else {
      // else f=NULL and we're just building a terminal node
      j=gamma;
      newnode=d;
      // fprev=root;
      // d=root;
    };

    // make a new terminal node
    previdx=alphaIndex(strg[j]); // potentially redundant code!
    newnode->children[previdx]=addNode(j,sLen-1,newnode);
    SUFFDEBUG { cout << "MAKE terminal: " << newnode->children[previdx] << "\n"; };

    if(newnode==root) { head_bot=EMPTY; }
    else { head_bot=i+newnode->total-1; };
    head_prev=newnode;
    // head_prev_el=d;

    // and we're done with step i.
  };

  ready=1; // whoopee!

};


/*******************************************************/
/* Low-level query functions to SuffixTrees            */

int SuffixTree::containsP(char *query,int length) {
  // is a substring included in the tree?
  // (start with <*query> and check <length> characters)

  Node *current=root;
  int posQ; // position in query string
  int posT; // position in tree's master string

  if(!current) return(0); // root==NULL?! this ain't no tree...

  for(posQ=0,posT=0; posQ<length; posQ++,posT++) {
    if(current->bot<posT) {
      // switch to next node
      current=current->children[alphaIndex(query[posQ])];
      if(!current) return(0);
      posT=current->top;
    };
    if(s[posT]!=query[posQ]) return(0);
  };

  return(1);
};

/*******************************************************/
/* Functions for integrity verification                */

int Node::countTerminals(int where[]) {
  int amIterminal=1;
  int res=0;
  for(int i=0;i<ALPHASIZE;i++) {
    amIterminal=amIterminal&&!children[i];
    if(children[i]) res+=children[i]->countTerminals(where);
  };
  if(where && amIterminal) where[this->total]=this->num;
  return(amIterminal?1:res);
};

int Node::countTerminals() {
  return(this->countTerminals(NULL));
};



/*******************************************************/
/* External query/wrapper functions for SuffixTrees    */


char alpha[]="acgtnN#$";

SuffixTree * getSuffixTree (int *hseq, char *hmask, int hseqlen,
			    int *mseq, char *mmask, int mseqlen);

SuffixTree * getSuffixTree (int *hseq, int hseqlen, int *mseq, int mseqlen) {
  return(getSuffixTree(hseq,NULL,hseqlen,mseq,NULL,mseqlen));
};

SuffixTree * getSuffixTree (int *hseq, char *hmask, int hseqlen,
			    int *mseq, char *mmask, int mseqlen) {
  SuffixTree *t=NULL;
  char *bigstr;

  if(!(bigstr=(char *)malloc((hseqlen+mseqlen+3)*sizeof(char)))) {
    cerr << "ERROR: Out of memory\n";
    exit(1);
  };
  
  for( int i=0;i<hseqlen;i++)
    bigstr[i]=((hmask&&hmask[i])?'n':alpha[hseq[i]]);
  bigstr[hseqlen]='#';
  for( int i=0;i<mseqlen;i++)
    bigstr[hseqlen+1+i]=((mmask&&mmask[i])?'N'
			 :((mseq[i]==4)?'N'
			   :alpha[mseq[i]]));
  bigstr[hseqlen+mseqlen+1]='$';
  bigstr[hseqlen+mseqlen+2]='\0';
  
  t=new SuffixTree(bigstr,alpha); // whoohoo. build the big tree. *poof*
  
  return(t);
};


void markNSeqsMasked(Node *n,SuffixTree *t,int regs,int boundaries[],
		     char *masks[],int tlength,int *markdata) {
  // NOTE: if tlength<=0, no masking is done, so each node has at
  // least one bit set. consequently, the "masks" array is NEVER accessed,
  // and may thus be set to NULL

  if(n->bot==boundaries[regs-1]) {
    // if terminal node
    int myreg;
    for(myreg=0;myreg<regs;myreg++) 
      if(n->top<=boundaries[myreg]) break;
    
    markdata[n->num]=1<<myreg;
    for(int i=0;i<tlength;i++)
      if((masks[myreg])[boundaries[regs-1]-n->total-
		       ((myreg==0)?-1:boundaries[myreg-1])+i] ||
	 t->alphaIndex(t->s[t->sLen-n->total+i]) > 3 )
	// if either base is masked, or an invalid base is found in the
	// match (like "n")
	{ markdata[n->num]=0; break; };

    return;
  };

  // non-terminal
  markdata[n->num]=0;
  int realChildren=0;
  for(int i=0;i<ALPHASIZE;i++)
    if(n->children[i]) {
      markNSeqsMasked(n->children[i],t,regs,boundaries,masks,tlength,markdata);
      markdata[n->num]|=markdata[n->children[i]->num];
      realChildren++;
    };

  if(realChildren==0) {
    cerr << "ERROR: Childless non-terminal node; suffix tree corrupted.\n";
    exit(1);
  };

  return;
}


void gatherMatchesFromNode(Node *n, int markdata[], int regs,
			   int boundaries[], int *idxs[], int lens[]) {
  if(!markdata[n->num]) return;

  if(n->bot==boundaries[regs-1]) {
    // terminal
    for(int i=0,m=markdata[n->num];i<regs;i++,m>>=1)
      if(m&0x1) {
	(idxs[i])[lens[i]]=boundaries[regs-1]-n->total-
	  ((i==0)?-1:boundaries[i-1]);
	lens[i]++;
      };
  } else {
    // non-terminal
    for(int i=0;i<ALPHASIZE;i++)
      if(n->children[i])
	gatherMatchesFromNode(n->children[i],markdata,regs,boundaries,
			      idxs,lens);
  };
};

int *idx;
int byIdx(const void *a,const void *b) {
  return(idx[*((int *)a)]-idx[*((int *)b)]);
};


void printTree( Node n, char *s, int indent )
{

  for( int i = 0; i < indent; i++ )
    cout << " ";
  for( int i = n.top; i <= n.bot; i++ )
    cout << s[i];
  cout << endl;

  for( int i = 0; i < ALPHASIZE; i++ )
    if( n.children[i] != NULL )
      printTree( *(n.children[i]), s, indent + n.bot - n.top + 1);

}


void mergeMatches( ilist **toRet, ilist ***childLists )
{

  for( int i = 0; i < 5; i++ ){
    for( int j = 0; j < 5; j++ ){
      if( childLists[j] ){
	toRet[0][i].join( childLists[j][0][i] );
	toRet[1][i].join( childLists[j][1][i] );
      }
    }
  }

}


void addMatches( mlist &matches, ilist ***childLists, int curDepth,
		 bool edgeMasked, int matchLen )
{

  if( edgeMasked ){
    match temp;

    temp.len = matchLen;
    for( int i_1 = 0; i_1 < 5; i_1++ ){
      if( childLists[i_1] == NULL )
	continue;
      for( int j_1 = 0; j_1 < 5; j_1++ ){
	iitem *first = childLists[i_1][0][j_1].head;

	while( first != NULL ){
	  temp.basePos = first->data;
	  for( int i_2 = 0; i_2 < 5; i_2++ ){
	    if( childLists[i_2] == NULL )
	      continue;
	    for( int j_2 = 0; j_2 < 5; j_2++ ){
	      if( j_2 == j_1 && j_1 != BASE_N )
		continue;
	      iitem *second = childLists[i_2][1][j_2].head;
	      
	      while( second != NULL ){
		temp.secPos = second->data;
		matches.add( temp );
		second = second->next;
	      }
	    }
	  }
	  first = first->next;
	}
      }
    }
  }
  else {
    match temp;

    temp.len = curDepth;
    for( int i_1 = 0; i_1 < 5; i_1++ ){
      if( childLists[i_1] == NULL )
	continue;
      for( int j_1 = 0; j_1 < 5; j_1++ ){
	iitem *first = childLists[i_1][0][j_1].head;

	while( first != NULL ){
	  temp.basePos = first->data;
	  for( int i_2 = 0; i_2 < 5; i_2++ ){
	    if( (i_2 == i_1 && i_1 != BASE_N) || childLists[i_2] == NULL )
	      continue;
	    for( int j_2 = 0; j_2 < 5; j_2++ ){
	      if( j_2 == j_1 && j_1 != BASE_N )
		continue;
	      iitem *second = childLists[i_2][1][j_2].head;
	      
	      while( second != NULL ){
		temp.secPos = second->data;
		matches.add( temp );
		second = second->next;
	      }
	    }
	  }
	  first = first->next;
	}
      }
    }
  }
}



ilist **getMatches( Node n, bp *seq, int *nextMask, int breakPt,
		    int minDepth, int curDepth,
                    bool isMatch, mlist &matches )
{

  // is the edge coming into this node masked? if so, note this and find the
  // amount of unmasked material in the edge.
  bool edgeMasked = false;
  int matchLen;

  matchLen = curDepth;
  if( n.bot >= nextMask[n.top] ){
    edgeMasked = true;
    matchLen -= n.bot - nextMask[n.top] + 1;
  }

  // if the edge is masked then this is as long a match as we're going to get.
  // if it's less than the minimum length, give up.
  if( edgeMasked && (matchLen < minDepth) )
    return NULL;

  // by the end of this function, toRet[i][j] will be a list of the starting
  // positions for every suffix going through this node which starts in
  // sequence in which is preceeded by either character j.
  ilist **toRet = new ilist* [2];
  // childLists[i] will be the toRet from the ith child.
  ilist ***childLists = new ilist** [5];

  toRet[0] = new ilist [5];
  toRet[1] = new ilist [5];

  // now get the toRets from the children. also note whether this is a
  // terminal node.
  bool isTerminal = true;

  for( int i = 0; i < 5; i++ ){
    if( n.children[i] != NULL ){
      isTerminal = false;
      childLists[i] = getMatches( *n.children[i], seq, nextMask, breakPt,
                                  minDepth, curDepth + n.children[i]->bot -
                                  n.children[i]->top + 1,
                                  !edgeMasked && isMatch && (i != 4),
				  matches );
    }
    else
      childLists[i] = NULL;
  }

  // there are two instances when we need to add a suffix to our lists. if
  // this node has a $-leaf and if this node is terminal.
  if( n.children[5] != NULL ){
    int pos = (n.children[5])->top - curDepth;
    char left = seq[pos - 1];
    char which_seq = 0;

    if( pos >= breakPt ){
      which_seq = 1;
      pos -= breakPt + 1;
    }
    if( childLists[4] == NULL ){
      childLists[4] = new ilist* [2];
      childLists[4][0] = new ilist [5];
      childLists[4][1] = new ilist [5];
    }
    childLists[4][which_seq][left].add( pos );
    //    cout << pos << "," << (which_seq ? '1' : '0') << " created" << endl;
  }

  if( isTerminal ){
    int pos = n.bot - curDepth + 1;
    char left = (pos == 0) ? 4 : seq[pos - 1];
    char which_seq = 0;

    if( pos >= breakPt ){
      which_seq = 1;
      pos -= breakPt + 1;
    }
#ifdef DEBUG
    assert( 0 <= left && left < 5 );
#endif
    toRet[which_seq][left].add( pos );
    //    cout << pos << "," << (which_seq ? '1' : '0') << " created" << endl;
  }
  else {
    if( isMatch && matchLen >= minDepth ){
      addMatches( matches, childLists, curDepth, edgeMasked, matchLen );
    }
    mergeMatches( toRet, childLists );
  }

  for( int i = 0; i < 5; i++ ){
    if( childLists[i] != NULL ){
      delete[] childLists[i][0];
      delete[] childLists[i][1];
      delete[] childLists[i];
    }
  }
  delete[] childLists;

  if( curDepth <= minDepth ){
    for( int i = 0; i < 5; i++ ){
      iitem *temp = toRet[0][i].head;

      while( temp != NULL ){
	iitem *next = temp->next;
	//	cout << temp->data << ",0 deleted" << endl;
	delete temp;
	temp = next;
      }
      temp = toRet[1][i].head;

      while( temp != NULL ){
	iitem *next = temp->next;
	//	cout << temp->data << ",1 deleted" << endl;
	delete temp;
	temp = next;
      }
    }
    delete[] toRet[0];
    delete[] toRet[1];
    delete[] toRet;
    toRet = NULL;
  }

  return toRet;

}
