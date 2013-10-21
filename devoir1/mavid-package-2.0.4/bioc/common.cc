

#include "common.h"


bp char2bp( char c )
{

  switch( c ){
  case 'a':
  case 'A':
    return BASE_A;
  
  case 'c':
  case 'C':
    return BASE_C;
  
  case 'g':
  case 'G':
    return BASE_G;
  
  case 't':
  case 'T':
      return BASE_T;
  
  default:
    return BASE_N;
  }

}


char bp2char( bp b )
{
  
  switch( b ){
  case BASE_A:
    return 'a';
  case BASE_C:
    return 'c';
  case BASE_G:
    return 'g';
  case BASE_T:
    return 't';
  case BASE_N:
    return 'n';
  default:
#if DEBUG
    assert(0);
#endif
    return 'n';
  }

}


bp complement( bp b )
{

  switch( b ){
  case BASE_A:
    return BASE_T;
  case BASE_C:
    return BASE_G;
  case BASE_G:
    return BASE_C;
  case BASE_T:
    return BASE_A;
  case BASE_N:
    return BASE_N;
  default:
#if DEBUG
    assert(0);
#endif
    return BASE_N;
  }

}


void arrayZero( char *arr, int sz )
{

  memset( arr, 0, sz );

}


void arrayZero( int *arr, int sz )
{

  for( register int i = 0; i < sz; arr[i++] = 0 );

}


void arrayZero( double *arr, int sz )
{

  for( register int i = 0; i < sz; arr[i++] = 0 );

}


void array_reverse( char *array, char *rev_array, int length )
{

  for( register int i = 0; i < length; i++ ) 
    rev_array[length - i - 1] = array[i];

}


void reverse_complement( bp *seq_array, int *rev_array, int array_len )
{

  int i = 0, j = array_len - 1;

  while( i < array_len ){
    rev_array[i] = complement( seq_array[j] );
    i++;
    j--;
  }

}

////////////////////////////////////////////////////////////////////////////////

// 'Quicksort':  Sorts the array 'vals' in the index range [p, q] and permutes
//             the corresponding array 'indices', in exactly the same way.
//             (Often, the indices array will be initialized (0, 1, 2, ..., n)
//             or something like that.)
// Note:  Corresponding values from 'vals' and 'array' (an originally parallel
//      array to 'vals') will be:  'vals[i]' <--> 'array[indices[i]]' because
//      'array' still has the original order.

void quicksort( int *vals, int *indices, int p, int r )
{

  if( p < r ){
    int q;
    q = partition( vals, indices, p, r );
    quicksort( vals, indices, p, q );
    quicksort( vals, indices, q + 1, r );
  }

}

////////////////////////////////////////////////////////////////////////////////

// 'Partition':  The function which is the essence of 'Quicksort'.

int partition( int *vals, int *indices, int p, int r )
{

  int x = vals[p];
  int i = p - 1, j = r + 1;

  while (1){
    do {
      j--;
    } while (vals[j] > x);
    do {
      i++;
    } while (vals[i] < x);
    if( i < j ){
      int tempval = vals[i];
      vals[i] = vals[j];
      vals[j] = tempval;
      int tempidx = indices[i];
      indices[i] = indices[j];
      indices[j] = tempidx;
    }
    else {    
      return j;
    }
  }

}


// binary searchs the array arr.
int binSearch( int *arr, int sz, int key )
{

  int high, low;

  for( low = -1, high = sz; high - low > 1; ){
    int i;

    i = (high + low) >> 1;
    if ( key <= arr[i] )
      high = i;
    else
      low  = i;

  }

  return high;

}
