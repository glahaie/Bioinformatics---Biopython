

#ifndef MAVID_COMMON_H
#define MAVID_COMMON_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <sstream>

using namespace std;

#define FALSE false
#define TRUE true
#define SMALL_VALUE 0.0000000000001     // Twelve 0's and then 1;
#define Infinity HUGE_VAL
#define OUT_OF_RANGE -1
#define BAD_CHAR_SIZE 256

#define GAPPED -1
#define UNWRITTEN -37

#define UNMASKED 0
#define MASKED 1

#define BASE_A 0
#define BASE_C 1
#define BASE_G 2
#define BASE_T 3
#define BASE_N 4
#define BASE_GAP 5

const char ALPHA[] = {'A','C','G','T','N','-'};

typedef char bp;

bp char2bp( char c );

char bp2char( bp b );

bp complement( bp b );

void arrayZero( char *arr, int sz );

void arrayZero( int *arr, int sz );

void arrayZero( double *arr, int sz );

void array_reverse( char *array, char *rev_array, int length );

void reverse_complement( int *seq_array, int *rev_array, int array_len );

void quicksort( int *vals, int *indices, int p, int r );

int partition( int *vals, int *indices, int p, int r );

int binSearch( int *arr, int sz, int key );

void warn( char *buff );

/*
#ifndef __GLIBCPP_INTERNAL_ALGOBASE_H
template <class T> inline T swap( T t1, T t2 )
{

  T temp;

  temp = t1;
  t1 = t2;
  t2 = temp;

}
#endif
*/

template <class T> inline T MIN( const T t1, const T t2 )
{

  return (t1 < t2) ? t1 : t2;

}

template <class T> inline T MAX( const T t1, const T t2 )
{

  return (t1 > t2) ? t1 : t2;

}

template <class T> inline T MIN( const T t1, const T t2, const T t3 )
{

  return (t1 < t2) ? ((t1 < t3) ? t1 : t3) : ((t2 < t3) ? t2 : t3);

}

template <class T> inline T MAX( const T t1, const T t2, const T t3 )
{

  return (t1 > t2) ? ((t1 > t3) ? t1 : t3) : ((t2 > t3) ? t2 : t3);

}

template <class T> inline T ABSDIFF( const T t1, const T t2 )
{

  if (t1 > t2)
    return t1 - t2;
  else
    return t2 - t1;

}

template <class T> inline void arrayInit( T val, T *arr, int sz )
{

  int i;
  
  for( i = 0; i < sz; arr[i++] = val );

}

template <class T> T arrayMax( T *array, int arraySize )
{

  T max = array[0];
  
  if( !array || arraySize <= 0 )
    return 0;
  for( register int i = 1; i < arraySize; i++ )
    if( array[i] > max) max = array[i];

  return max;

}    

template <class T> T arrayMin(T *array, int arraySize)
{
  T value, min = array[0];                          
  
  if (!array || arraySize <= 0) return 0;
  for (register int i = 1; i < arraySize; i++)
    if ((value = array[i]) < min) min = value;
  return min;
}



#endif
