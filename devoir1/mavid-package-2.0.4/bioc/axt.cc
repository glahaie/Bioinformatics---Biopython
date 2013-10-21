
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "axt.h"


AXTBlock::AXTBlock()
{

  baseLine = secLine = NULL;

}


void AXTBlock::readBlock( FILE *in )
{

  int lineStart;

  if( baseLine != NULL ){
    delete baseLine;
    delete secLine;
  }

  fscanf( in, "%u %255s %u %u %255s %u %u %c", &index, baseName, &baseStart,
	  &baseEnd, secName, &secStart, &secEnd, &strand );
  while( getc( in ) != '\n' );
  lineStart = ftell( in );
  linelen = 0;
  while( getc( in ) != '\n' )
    linelen++;
  baseLine = new char [linelen];
  secLine = new char [linelen];
  fseek( in, lineStart, SEEK_SET );
  linelen = 0;
  for( char c = getc( in ); c != '\n'; c = getc ( in ) )
    baseLine[linelen++] = c;
  linelen = 0;
  for( char c = getc( in ); c != '\n' && !feof(in) ; c = getc ( in ) )
    secLine[linelen++] = c;
  // advance to next block or end of file
  char c;

  while( !isdigit(c = getc( in )) && !feof( in ));
  if( isdigit( c ) )
    ungetc( c, in );

  next = NULL;

}


void AXTBlock::writeBlock( FILE *out )
{

  fprintf( out, "%u %s %u %u %s %u %u %c %f\n", index, baseName, baseStart, baseEnd,
	   secName, secStart, secEnd, strand, score );
  for( int i = 0; i < linelen; i++ )
    fprintf( out, "%c", baseLine[i] );
  fprintf( out, "\n" );
  for( int i = 0; i < linelen; i++ )
    fprintf( out, "%c", secLine[i] );
  fprintf( out, "\n\n" );

}


void AXTBlock::printBlock()
{

  printf( "%u %s %u %u %s %u %u %c %f\n", index, baseName, baseStart, baseEnd,
	  secName, secStart, secEnd, strand, score );
  for( int i = 0; i < linelen; i++ )
    printf( "%c", baseLine[i] );
  printf( "\n" );
  for( int i = 0; i < linelen; i++ )
    printf( "%c", secLine[i] );
  printf( "\n\n" );

}


char char2idx( char c )
{

  switch( c ){
  case 'a':
  case 'A':
    return 0;

  case 'c':
  case 'C':
    return 1;

  case 'g':
  case 'G':
    return 2;

  case 't':
  case 'T':
    return 3;

  default:
    return 4;
  }

}


AXTFile AXTBlock::splitBlock( int startIdx, float **scores, float gapOpen,
			      float gapExtend, float threshold )
{

  AXTFile toRet;
  AXTBlock *toAdd = new AXTBlock;
  bool inBaseGap = false, inSecGap = false;
  float maxScore = 0, curScore = 0;

  for( int curI = 0, startI = 0, maxI = 0,
	 curBP = baseStart, curSP = secStart, startBP = 0, startSP = 0,
	 maxBP = 0, maxSP = 0; curI <= linelen; curI++ ){
    if( curScore == 0 || curI == linelen ){
      if( maxScore >= threshold ){
	toAdd->index = startIdx;
	strcpy( toAdd->baseName, baseName );
	strcpy( toAdd->secName, secName );
	toAdd->baseStart = startBP;
	toAdd->secStart = startSP;
	toAdd->baseEnd = maxBP;
	toAdd->secEnd = maxSP;
	toAdd->strand = strand;
	toAdd->score = maxScore;
	toAdd->linelen = maxI - startI + 1;
	toAdd->baseLine = new char [toAdd->linelen];
	strncpy( toAdd->baseLine, baseLine + startI, toAdd->linelen );
	toAdd->secLine = new char [toAdd->linelen];
	strncpy( toAdd->secLine, secLine + startI, toAdd->linelen );
	toRet.addBlock( toAdd );

	toAdd = new AXTBlock;
	startIdx++;
	curI = maxI + 1;
	curBP = maxBP;
	curSP = maxSP;
	if( baseLine[curI - 1] != '-' ){
	  curBP++;
	}
	if( secLine[curI - 1] != '-' ){
	  curSP++;
	}
	curScore = maxScore = 0;

      }
      else if( curI == linelen ){
	break;
      }
      startBP = curBP;
      startSP = curSP;
      startI = curI;
    }
    // time to update the score
    if( baseLine[curI] == '-' ){
      if( inBaseGap )
	curScore += gapExtend;
      else
	curScore += gapOpen;
      inBaseGap = true;
      inSecGap = false;
    }
    else if( secLine[curI] == '-' ){
      if( inSecGap )
	curScore += gapExtend;
      else
	curScore += gapOpen;
      inBaseGap = false;
      inSecGap = true;
    }
    else {
      inBaseGap = inSecGap = false;
      curScore += scores[char2idx(baseLine[curI])][char2idx(secLine[curI])];
    }
    curScore = (curScore < 0) ? 0 : curScore;
    // done updating score
    if( maxScore < curScore ){
      maxScore = curScore;
      maxI = curI;
      maxBP = curBP;
      maxSP = curSP;
    }
    if( baseLine[curI] != '-' ){
      curBP++;
    }
    if( secLine[curI] != '-' ){
      curSP++;
    }    
  }

  return toRet;

}


//
//  AXTFile functions
//

AXTFile::AXTFile()
{

  head = tail = NULL;

}

AXTFile::~AXTFile()
{

}

void AXTFile::addBlock( AXTBlock *toAdd )
{

  if( head == NULL ){
    head = tail = toAdd;
    toAdd->next = NULL;
  }
  else {
    tail->next = toAdd;
    tail = toAdd;
    tail->next = NULL;
  }

}


void AXTFile::readAXTFile( char *fileName )
{

  FILE *in;

  in = fopen( fileName, "r" );

  if( !in )
    return;

  while( !feof( in ) ){
    AXTBlock *newBlock = new AXTBlock;

    newBlock->readBlock( in );
    addBlock( newBlock );
  }

}


void AXTFile::printAXTFile()
{

  for( AXTBlock *curBlock = head; curBlock != NULL; curBlock = curBlock->next )
    curBlock->printBlock();

}


AXTFile AXTFile::splitAXT( float **scores, float gapOpen, float gapExtend,
			   float threshold )
{

  AXTFile toRet;
  int curIdx = 1;

  for( AXTBlock *curBlock = head; curBlock != NULL;
       curBlock = curBlock->next ){
    AXTFile temp;

    temp = curBlock->splitBlock( curIdx, scores, gapOpen, gapExtend,
				 threshold );
    if( temp.head == NULL )
      continue;
    if( toRet.head == NULL )
      toRet.head = temp.head;
    else
      toRet.tail->next = temp.head;
    toRet.tail = temp.tail;
    curIdx = toRet.tail->index;

  }

  return toRet;

}

