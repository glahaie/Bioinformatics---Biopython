
#include <string.h>
#include "seq.h"

char DNA2alpha( char c )
{

  if( c < 6 )
    return DNAHASH[(int)c];
  else
    return 'N';

}

char Prot2alpha( char c )
{

  if( c < 2 )
    return PROTHASH[(int)c];
  else
    return 'X';

}

char alpha2DNA( char c )
{

  switch( c ){
  case 'A':
  case 'a':
    return DNA_A;

  case 'C':
  case 'c':
    return DNA_C;

  case 'G':
  case 'g':
    return DNA_G;

  case 'T':
  case 't':
    return DNA_T;

  case '-':
    return DNA_GAP;

  default:
    return DNA_N;

  }

}


char alpha2Prot( char c )
{

  return PROT_A;

}


bool isDNA( char c )
{

  return (strchr( "ACGTNacgtn", c ) != NULL );

}


bool isProt( char c )
{

  return (strchr( "Aa", c ) != NULL );

}
