
#ifndef BIOC_AXT
#define BIOC_AXT

class AXTBlock;
class AXTFile;

class AXTBlock {

 public:
  int index;
  char baseName[256], secName[256];
  int baseStart, baseEnd, secStart, secEnd;
  char strand;
  float score;

  int linelen;
  char *baseLine, *secLine;

  AXTBlock *next;

  AXTBlock();
  ~AXTBlock();
  void readBlock( FILE *in );
  void writeBlock( FILE *out );
  void printBlock();
  void printSubBlock( int startI, int stopI, int baseOffset, int secOffset,
		      int score );
  AXTFile splitBlock( int startIdx, float **scores, float gapOpen,
		      float gapExtend, float theshold );

};


class AXTFile {

 public:
  AXTBlock *head, *tail;

  AXTFile();
  ~AXTFile();
  void addBlock( AXTBlock *toAdd );
  void readAXTFile( char *fileName );
  void writeAXTFile( char *fileName);
  void printAXTFile();
  AXTFile splitAXT( float **scores, float gapOpen, float gapExtend,
		    float threshold );
};

#endif
