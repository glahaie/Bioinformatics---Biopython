// Main SuffixTree class
#ifndef SUFFTREE_H
#define SUFFTREE_H

#include <iostream>
#include <strings.h>
#include <stdlib.h>

using namespace std;

#include "list.h"

struct match {
  int basePos, secPos, len;
};

typedef list<int> ilist;
typedef listelem<int> iitem;
typedef list<match> mlist;
typedef listelem<match> mitem;


#define EMPTY -1
#define STREE_VERSION 2 // to mark stored-tree files
#define SUFFDEBUG if(0)

// CHANGE ME!!!!!!!!!!!
#define ALPHASIZE 6

class Node {
public:

  Node() { };

  inline void setParams(int t,int b,Node * p,int n) {
    // ASSUMPTION: children is allocated already
    top=t;
    bot=b;
    parent=p;
    num=n;

    // this works if EMPTY==-1, as is:
    total=(p?(p->total+(b-t+1)):0);

    if(p) {
      if(p->bot<0) {
	init=t;
      } else {
	init=p->init;
      };
    } else {
      init=-1;
    };

    for(int i=0;i<ALPHASIZE;i++) { children[i]=NULL; };
    slink=NULL;
    
    SUFFDEBUG if(top>bot && bot!=EMPTY) {
      cerr << "Strange node request: " << top << "-" << bot
	   << " / parent " << parent << "\n";
      exit(1);
    };

    return;
  };

  ~Node() {};

  // human-readable output
  void printReadable(ostream & out);
  friend ostream& operator<<(ostream&, const Node&);

  // verification:
  int countTerminals();
  int countTerminals(int where[]);

  // top and bottom of the substring the edge coming ***INTO*** the
  // node corresponds to, as indices on the array
  int top,bot;
  int total,num,init;
  
  Node ** children; // as indices into the node array
  Node * slink; // suffix link
  Node * parent;

};

typedef struct NodeList {
  Node *n;
  NodeList *nxt;
} NodeList;


class SuffixTree {
private:
  Node * addNode(int t,int b,Node * p);
  void buildAlphaHash() {
    for(int i=0;i<ALPHASIZE;i++) alphaHash[(int)alpha[i]]=i;
  };

public:
  int ready;
  SuffixTree() { ready=0; };
  SuffixTree(const char *strg, const char *alphabet);
  ~SuffixTree() {
    //delete root;
    delete [] myNodes;
    delete [] myNodeLinks;
  };

  int  containsP(char *query,int length);

  // human-readable output
  void printReadable(ostream & out);
  friend ostream& operator<<(ostream&, const SuffixTree&);

  // machine-readable output
  int saveTo(ostream & out);

  // input from saved tree
  int loadFrom(istream & in);

  // this should be rewritten for arbitrary alphabets: 
  inline int alphaIndex(char c) {
    return alphaHash[(unsigned int)c];
  };

  Node * root;
  int nNodes;

  // previous head
  Node * head_prev; 

  // contracted/extended loci of prev. heads in pre-prev. tree
  Node * head_prev_cl;
  // Node * head_prev_el;
  int head_top, head_bot;

  const char * s;
  int sLen;

  // for arbirary alphabets:
  const char * alpha;
  int alphaHash[256];

  // Node and child-pointer storage
  Node  * myNodes;
  Node ** myNodeLinks;
};


ilist **getMatches( Node n, bp *seq, int *nextMask, int breakPt,
		    int minDepth, int curDepth,
                    bool isMatch, mlist &matches );

#endif // SUFFTREE_H

