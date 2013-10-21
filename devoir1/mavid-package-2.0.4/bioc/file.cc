
#include <stdio.h>

// This function reads in a file and stores it(null-terminated) at the
// address returned.
char *readFile( char *fileName )
{

  FILE *in;
  char *toReturn;
  int len;

  in = fopen( fileName, "r" );
  if( !in )
    return NULL;

  // get the file length
  fseek( in, 0, SEEK_END );
  len = ftell( in );
  fseek( in, 0, SEEK_SET );

  // read in the file
  toReturn = new char [len + 1];
  fread( toReturn, 1, len, in );
  toReturn[len] = 0;

  return toReturn;

}
