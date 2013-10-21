
#ifndef BIOC_SEQ
#define BIOC_SEQ

/*
** Enums, defines, and whatnot
*/

enum { DNA_A, DNA_C, DNA_G, DNA_T, DNA_N, DNA_GAP };
enum { PROT_A };

#define DNAHASH "ACGTN-"
#define PROTHASH "A"

/*
** Structs, classes, and whatnot
*/

/*
** Function prototypes
*/

char alpha2DNA( char c );
char alpha2Prot( char c );
char DNA2alpha( char c );
char Prot2alpha( char c );

bool isDNA( char alpha );
bool isProt( char alpha );

#endif
