
#ifndef BIOC_MATRICES
#define BIOC_MATRICES

/*
** Enums, defines, and whatnot
*/

/*
** Structs, classes, and whatnot
*/

/*
** Function prototypes
*/

double **hoxd_scores( );

double **simple_scores( );

void sc_fix(double *s,double t);

void sc_fix_gapped(double *s,double t);

void f_fix(double *s,double t);

void p_fix(double *s,double t);

void f_fix_gapped(double *s,double t);

void readScores( char *fileName, float **&scores, float &gapOpen,
                 float &gapExtend );

void scoreMatrix( double **&mat, double t );

void transition_matrix_gapped( double **&mat, double t );

void mutationMatrix( double **mat, double t );

void score_matrix(double **&scores,double t);

void gapped_score_matrix(double **&scores,double t);

#endif
