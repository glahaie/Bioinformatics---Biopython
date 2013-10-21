
#include <stdlib.h>

template <class T>
class listelem {

 public:
  T data;
  listelem<T> *next;
  listelem(T dat) { data = dat; next = NULL; };
  listelem() { next = NULL; };
  
};

template <class T>
class list {

 public:
  listelem<T> *head, *tail;

  list(void) { head = tail = NULL; };
  list(T init) { head = tail = new listelem<T>(init); };
  void add( T data )
    {
      if( tail == NULL )
	head = tail = new listelem<T>(data);
      else {
	tail->next = new listelem<T>(data);
	tail = tail->next;
      }
    };
  void join( list<T> toJoin )
    {
      if( tail == NULL ){
	head = toJoin.head;
	tail = toJoin.tail;
      }
      else {
	if( toJoin.head == NULL )
	  return;
	tail->next = toJoin.head;
	tail = toJoin.tail;
      }
    };
};

