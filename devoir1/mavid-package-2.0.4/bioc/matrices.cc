#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
/* #include "nrutil.h" */

#define ALIGN 3000000
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


/* ***************** */		
/* Memory allocation */
/* ***************** */		


double **hoxd_scores( )
{

  double **to_ret;

  to_ret = new double* [6];
  for( int i = 0; i < 6; i++ )
    to_ret[i] = new double [6];

  to_ret[0][0] =   91; to_ret[0][1] = -114; to_ret[0][2] =  -31; to_ret[0][3] = -123; to_ret[0][4] = -123; to_ret[0][5] = -123;
  to_ret[1][0] = -114; to_ret[1][1] =  100; to_ret[1][2] = -125; to_ret[1][3] =  -31; to_ret[1][4] = -125; to_ret[1][5] = -125;
  to_ret[2][0] =  -31; to_ret[2][1] = -125; to_ret[2][2] =  100; to_ret[2][3] = -114; to_ret[2][4] = -114; to_ret[2][5] = -114;
  to_ret[3][0] = -123; to_ret[3][1] =  -31; to_ret[3][2] = -114; to_ret[3][3] =   91; to_ret[3][4] = -123; to_ret[3][5] = -123;
  to_ret[4][0] = -123; to_ret[4][1] = -125; to_ret[4][2] = -125; to_ret[4][3] = -123; to_ret[4][4] = -125; to_ret[4][5] = -125;
  to_ret[5][0] = -123; to_ret[5][1] = -125; to_ret[5][2] = -125; to_ret[5][3] = -123; to_ret[5][4] = -125; to_ret[5][5] = -125;

  return to_ret;

}


double **simple_scores( )
{

  double **to_ret;

  to_ret = new double* [6];
  for( int i = 0; i < 6; i++ )
    to_ret[i] = new double [6];

  to_ret[0][0] =  1; to_ret[0][1] = -1; to_ret[0][2] = -1; to_ret[0][3] = -1; to_ret[0][4] = -1; to_ret[0][5] = -1;
  to_ret[1][0] = -1; to_ret[1][1] =  1; to_ret[1][2] = -1; to_ret[1][3] = -1; to_ret[1][4] = -1; to_ret[1][5] = -1;
  to_ret[2][0] = -1; to_ret[2][1] = -1; to_ret[2][2] =  1; to_ret[2][3] = -1; to_ret[2][4] = -1; to_ret[2][5] = -1;
  to_ret[3][0] = -1; to_ret[3][1] = -1; to_ret[3][2] = -1; to_ret[3][3] =  1; to_ret[3][4] = -1; to_ret[3][5] = -1;
  to_ret[4][0] = -1; to_ret[4][1] = -1; to_ret[4][2] = -1; to_ret[4][3] = -1; to_ret[4][4] =  1; to_ret[4][5] = -1;
  to_ret[5][0] = -1; to_ret[5][1] = -1; to_ret[5][2] = -1; to_ret[5][3] = -1; to_ret[5][4] = -1; to_ret[5][5] =  1;

  return to_ret;

}


/* Dynamic memory allocation: doubles */
double *dmalloc(long n)
{
  double *x;

  if((x = (double *)malloc((unsigned)(n * sizeof(double)))) == NULL) {
    printf("Error: Couldn't get %ld doubles.\n", n);
    exit(1);
  }
  return x;
}


/* Dynamic memory allocation: integers */
int *dmalloci(int n) 
{
  int *x;

  if((x = (int *)malloc((unsigned)(n * sizeof(int)))) == NULL) {
    printf("Error: Couldn't get %d integers.\n", n);
    exit(1);
  }
  return x;
}


/* Dynamic memory allocation: long integers */
long *dmallocs(long n) 
{
  long *x;

  if((x = (long *)malloc((unsigned)(n * sizeof(long)))) == NULL) {
    printf("Error: Couldn't get %ld long integers.\n", n);
    exit(1);
  }
  return x;
}


/* Dynamic memory allocation: characters */
char *dmallocc(int n) 
{
  char *x;

  if((x = (char *)malloc(n * sizeof(char))) == NULL) {
    printf("Error: Couldn't get %d characters.\n", n);
    exit(1);
  }
  return x;
}


/* Dynamic memory allocation: characters pointers */
char **dmalloccp(int n) 
{
  char **x;

  if((x = (char **)malloc(n * sizeof(char *))) == NULL) {
    printf("Error: Couldn't get %d character pointers.\n", n);
    exit(1);
  }
  return x;
}


/* Print multiple alignment */
void prmal(FILE *fl, long left, long right, char *x1, char *x2, char *x3)
{
  long i, j, nrow;
  
  nrow = (right-left+1)/70;
  for (i = left; i < left + nrow*70; i+= 70) {
    for (j = i; j < i+70; j++) {
      fprintf(fl, "%c", x1[j]);
    }
    fprintf(fl, "\n");
    for (j = i; j < i+70; j++) {
      fprintf(fl, "%c", x2[j]);
    } 
    fprintf(fl, "\n");
    for (j = i; j < i+70; j++) {
      fprintf(fl, "%c", x3[j]);
    }
    fprintf(fl, "\n\n");
  }    
  for (j = left + nrow*70; j <= right; j++) {
    fprintf(fl, "%c", x1[j]);
  }
  fprintf(fl, "\n");
  for (j = left + nrow*70; j <= right; j++) {
    fprintf(fl, "%c", x2[j]);
  }
  fprintf(fl, "\n");
  for (j = left + nrow*70; j <= right; j++) {
    fprintf(fl, "%c", x3[j]);
  }
  fprintf(fl, "\n\n\n");
}


/* *********************** */
/* Random number generator */
/* *********************** */

#define M 714025
#define IA 1366
#define IC 150889

/* Generate a sample from U(0,1) */
double runiform()
{
  static long iy,ir[98];
  static int iff=0;
  int j;
  static time_t now;
  
  if (iff==0) now=time(0);
  
  if (now <0 || iff ==0) {
    iff=1;
    if ((now=(IC-(now)) % M) < 0) now=-(now);
    for (j=1;j<=97;j++) {
      now=(IA*(now)+IC) % M;
      ir[j]=(now);
    }
    now=(IA*(now)+IC) % M;
    iy=(now);
  }
  j = 1 + (int)(97.0*(double)iy/(double)M);
  if (j>97 || j < 1) printf("RUNI: This cannot happen");
  iy=ir[j];
  now=(IA*(now)+IC) %M;
  ir[j]=(now);
  return (double) iy/M;
}


/* Generate n samples from U(0,1) */ 
double *rdm(long n)
{
  long i;
  double *x, *dmalloc(long), runiform();

  x = dmalloc(n);
  for (i=0;i<n;i++)
    x[i] = runiform();
  return x; 
}


/* Generate a sample from Multinomial(N,pr) */
/* ct: count vector of length N             */
/* pr: probability vector of length n       */
void multinomial(int N, int *ct, double *pr, int n)
{
  int i, cell, step;
  double *cf, *x, *dmalloc(long), *rdm(long);

  cf = dmalloc(n);
  x = rdm(N);
  cf[0] = pr[0];
  for (i = 1; i < n; i++) 
    cf[i] = cf[i-1] + pr[i];
  for (i = 0; i < n; i++) 
    ct[i] = 0; 
  for (i = 0; i < N; i++) {
    cell = 0;
    step = 0;
    while (x[i] > cf[step]) {
      cell++;
      step++;
    }
    ct[cell] += 1; 
  }
  free(cf); 
  free(x);
}


/* Generate n independent samples from X = Geometric(p) */
/* p is the probability of success */
/* X is the number of independent trials needed to get the first success */
/* P(X=1) = p, P(X=2) = (1-p) p, P(X=3) = (1-p)^2 p, etc */
int *geometric(int n, double p)
{
  int *geo, i, *dmalloci(int);
  double *x, *rdm(long);

  x = rdm(n);
  geo = dmalloci(n);
  for (i = 0; i < n; i++)
    geo[i] = (int) (1 + log(1-x[i]) / log(1-p) );
  return geo;
}



/* ************* */
/* Miscellaneous */
/* ************* */


/* Print multiple alignment */
void pma_mfa(FILE *fl, char **name, char **seq, int left, int right, int n)
{
  char *x, **p;
  int ct, gap, w = 70, len;
  
  len = right - left + 81;
  if ((x = (char *)malloc((unsigned)(n*len*sizeof(char)))) == NULL) {
     printf("Not enough memory\n");
     exit(1);
  }
  if ((p = (char **)malloc((unsigned)(n*sizeof(char *)))) == NULL) {
     printf("Not enough memory\n");
     exit(1);
  }
  for( int j = 0; j < n; j++ )
    p[j] = x + j*len;
  ct = 0;
  for( int i = left; i <= right; i++ ){
    gap = 0;
    for( int j = 0; j < n; j++ ){
      if (seq[j][i] == '-')
        gap++;
    }
    if( gap < n ){
      ct++;
      for( int j = 0; j < n; j++)
        *(p[j]++) = seq[j][i];
      if (ct == w) {
        for( int j = 0; j < n; j++ )
          *(p[j]++) = '\n';
        ct = 0;
      } 
    }
  }
  for( int j = 0; j < n; j++ ){
    fprintf(fl, "%s\n", name[j]);
    *p[j] = '\0';
    fprintf(fl, "%s\n", x+j*len);
  }
}


/* Print multiple alignment */
void oldpma_mfa(FILE *fl, char **name, int left, int right, char *x[], int n)
{
  int i, j, k, nrow, w = 70;
  
  nrow = (right-left+1)/w;
  for (j = 0; j < n; j++) {
    fprintf(fl, "%s\n", name[j]);
    for (i = left; i < left + nrow*w; i+= w) {
      for (k = i; k < i+w; k++) 
        fprintf(fl, "%c", x[j][k]);
      fprintf(fl, "\n");
    }
    for (k = left + nrow*w; k <= right; k++)
      fprintf(fl, "%c", x[j][k]);
    fprintf(fl, "\n");
  }    
}


/* Print multiple alignment */
void pma_phy(FILE *fl, char **name, int left, int right, char *x[], int n)
{
  int i, j, k, nrow, w = 70;
  
  nrow = (right-left+1)/w;
  for (i = left; i < left + nrow*w; i+= w) {
    for (j = 0; j < n; j++) {
      for (k = i; k < i+w; k++) 
        fprintf(fl, "%c", x[j][k]);
      fprintf(fl, "\n");
    }
    fprintf(fl, "\n");
  }
  for (j = 0; j < n; j++) {
    for (k = left + nrow*w; k <= right; k++)
      fprintf(fl, "%c", x[j][k]);
    fprintf(fl, "\n");
  }    
  fprintf(fl, "\n\n");
}


/* Print multiple alignment */
void prmual(FILE *fl, int left, int right, char *x[], int n)
{
  int i, j, k, nrow;
  
  nrow = (right-left+1)/70;
  for (i = left; i < left + nrow*70; i+= 70) {
    for (j = 0; j < n; j++) {
      for (k = i; k < i+70; k++) 
        fprintf(fl, "%c", x[j][k]);
      fprintf(fl, "\n");
    }
    fprintf(fl, "\n");
  }    
  for (j = 0; j < n; j++) {
    for (k = left + nrow*70; k <= right; k++)
      fprintf(fl, "%c", x[j][k]);
    fprintf(fl, "\n");
  }
  fprintf(fl, "\n\n");
}


/* Get number of lines of a file */
int n_line(FILE *fl)
{
  int c, num = 0;
  
  while ((c = fgetc(fl)) != EOF)
    if (c == '\n')
      num++;
  return num;
}
  

/* Get number of sequences of a mfa file */
int n_seq(FILE *fl)
{
  char x[81];
  int num = 0;
  
  while (!feof(fl)) {
    fgets(x, 81, fl);
    if (x[0] == '>')
      num++;
  }
  return num;
}
  

/* *************************** */
/* Basic operations on doubles */
/* *************************** */

/* Copy an array to another */
void dblcp(double *des, double *sou, long n) 
{
  while (n-- > 0)
    *des++ = *sou++;
}

/* Add array x2 to x1 */
void dbladd(double *x1, double *x2, long n) 
{
  while (n-- > 0)
    *x1++ += *x2++;
}

/* Multiply x1 by a scalar s */
void dblmul(double *x1, double s, long n) 
{
  while (n-- > 0)
    *x1++ *= s;
}

/* Divide x1 by a scalar s */
void dbldiv(double *x1, double s, long n) 
{
  while (n-- > 0)
    *x1++ /= s;
}

/* vector sum */
double vec_sum(double *x, long n, long skip)
{
  double s = 0;

  for (; n > 0; n -= skip, x += skip)
    s += *x;
  return s;
}
  
/* Row sum */
void row_sum(double *rs, double *x, long rn, long cn)
{
  long t_rn = rn, n = rn * cn;
  double vec_sum(double *, long, long);
 
  while (t_rn-- > 0)
    *rs++ = vec_sum(x++, n, rn); 
}
     

/* ************************ */
/* Basic operations on ints */
/* ************************ */

/* Copy an array of ints of length n to another */
void intcp(int *des, int *sou, int n) 
{
  while (n-- > 0)
    *des++ = *sou++;
}

/* Add an array of ints x2 to x1 */
void intadd(int *x1, int *x2, int n) 
{
  while (n-- > 0)
    *x1++ += *x2++;
}

/* vector sum */
int vec_isum(int *x, int n, int skip)
{
  int s = 0;

  for (; n > 0; n -= skip, x += skip)
    s += *x;
  return s;
}
  
/* Row sum */
void row_isum(int *rs, int *x, int rn, int cn)
{
  int t_rn = rn, n = rn*cn, vec_isum(int *, int, int);
 
  while (t_rn-- > 0)
    *rs++ = vec_isum(x++, n, rn); 
}
     


/* Calculate sum c(i)*log(f(i)) */
double multi_llh(int *c, double *f, int n) 
{
  double l = 0;

  while (n-- > 0)
    l += *c++ * log(*f++);
  return l;
}   


/* Transpose */
void transpose(double *xt, double *x, long nrow, long ncol)
{
  long i, j;

  for (i = 0; i < ncol; i++)
    for (j = 0; j < nrow; j++)
      xt[j*ncol+i] = x[i*nrow+j];
}


/* l1 distance between two vectors */
double dist_l1(double *x, double *y, long n)
{
  double d = 0;

  while (n-- > 0)
    d += fabs(*x++ - *y++);
  return d;
} 


/* Desktop Calculator Algorithm */
void desk(double *x, long n, double *mean, double *var)
{
  long t_n = n;

  for (*mean = *var = 0; t_n > 0; t_n--, x++) {
    *mean += *x;
    *var += *x * *x;
  }
  *mean /= n;
  *var = *var / n - *mean * *mean ; 
}


/* Correlation coefficient */
double corrcoef(double *x, double *y, long n)
{
  long i;
  double corr = 0, mean_x, mean_y, var_x, var_y, *dmalloc(long);

  desk(x, n, &mean_x, &var_x);
  desk(y, n, &mean_y, &var_y);
  
  for (i = 0; i < n; i++)
    corr += (x[i] - mean_x) * (y[i] - mean_y);
  corr /= n * sqrt(var_x) * sqrt(var_y);
  return corr;
}     


/* Maximum */
double getmax(double *x, long n, long *ind)
{
  long i;
  double max; 

  max = x[0];
  *ind = 0;
  for (i = 1; i < n; i++) {
    if (x[i] > max) {
      max = x[i];
      *ind = i;
    }
  }
  return max;
}


double gmax(double *x, int n, int *ind)
{
  int i;
  double max; 

  max = x[0];
  *ind = 0;
  for (i = 1; i < n; i++) {
    if (x[i] > max) {
      max = x[i];
      *ind = i;
    }
  }
  return max;
}


/* Minimum */
double getmin(double *x,long n,long *ind)
{
  long i;
  double min; 

  min = x[0];
  *ind = 0;
  for (i = 1; i < n; i++) {
    if (x[i] < min) {
      min = x[i];
      *ind = i;
    }
  }
  return min;
}


double gmin(double *x, int n, int *ind)
{
  int i;
  double min; 

  min = x[0];
  *ind = 0;
  for (i = 1; i < n; i++) {
    if (x[i] < min) {
      min = x[i];
      *ind = i;
    }
  }
  return min;
}


/* printing matrix of integers */
void printm_int(int *x, int n, int p)
{
  int i, j;

  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++)
      printf("%d ", x[j*n+i]);
    printf("\n");
  }
}
 

/* printing matrix */
void printm(double *x, long n, long p)
{
  long i, j;

  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++)
      printf("%lf ", x[j*n+i]);
    printf("\n");
  }
}
 

/* printing matrix to a file*/
void fprintm(FILE *fl, double *x, long n, long p)
{
  long i, j;

  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++)
      fprintf(fl, "%lf ", x[j*n+i]);
    fprintf(fl,"\n");
  }
}
 

/* scanning matrix */
void scanm(double *x, long n, long p)
{
  long i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < p; j++)
      scanf("%lf", x+j*n+i);
}


/* scanning matrix from a file*/
void fscanm(FILE *fl, double *x, long n, long p)
{
  long i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < p; j++)
      fscanf(fl, "%lf", x+j*n+i);
}


/* load a sequence */
void read_seq(FILE *fl, char *x)
{
  int c, on;
  on = 1;
  while ((c = fgetc(fl)) != EOF) {
    if (c == '>')
      on = 0;
    if (c != ' ' && c != '\n' && on == 1)
      *x++ = c;
    if (on == 0 && c == '\n')
      on = 1;
  }
  *x = '\0';
}


/* kick gaps */
void kickgap(char *seq, int nl, char *nseq)
{
  char **ps, **pns, **dmalloccp(int);
  int len, nlen, is_gap; 

  len = strlen(seq) / nl;
  ps = dmalloccp(nl);
  pns = dmalloccp(nl);
  ps[0] = seq;
  pns[0] = nseq;
  for( int i = 1; i < nl; i++ ){
    ps[i] = ps[i-1] + len;
    pns[i] = pns[i-1] + len;
  }
  nlen = 0;
  while (*ps[nl-1] != '\0') {
    is_gap = 0;
    for( int i = 0; i < nl; i++){
      if (*ps[i] == '-') {
        is_gap = 1;
        break;
      }
    }
    if (is_gap == 0) {
      for( int i = 0; i < nl; i++ )
        *pns[i]++ = *ps[i];
      nlen++;
    }
    for( int i = 0; i < nl; i++ )
      ps[i]++;
  }
  free(ps); free(pns);
}


/* Scalar Product */
double dot(double *v, double *w, long n, long iv, long iw)
{
  double sum = 0;  

  for (; n > 0; n--, v += iv, w += iw)
    sum += *v * *w;
  return sum;
}


/* Matrix X'X*/
void psd(double *x, long n, long p, double *xtx)
{
  long i, j;
  double dot(double *, double *, long, long, long);

  for (i = 0; i < p; i++)
    for (j = 0; j < i+1; j++)
      xtx[i*p+j] = xtx[j*p+i] = dot(x+i*n, x+j*n, n, 1, 1);
}


/* Matrix Multiplication */
void matmul(double *x, double *y, long l, long m, long n, double *z)
{
  long i, j;
  double dot(double *, double *, long, long, long);
 
  for (i = 0; i < l; i++)
    for (j = 0; j < n; j++)
      z[j*l+i] = dot(x+i, y+j*m, m, l, 1);
}


/* Cholesky Decomposition */
void cholesky(double *x, long n)
{
  long i, j, k;
  double tmp, e;

  e = pow(10., -8.);
  for (i = 0; i < n; i++) { 
    tmp = 0;
    for (k = 0; k < i; k++) 
      tmp += x[i*n+k] * x[i*n+k];
    x[i*n+i] = sqrt(x[i*n+i] - tmp);
    if (fabs(x[i*n+i]) < e) {
      printf("Division by %lf in column %ld\n", x[i*n+i], i);
      break;
    }
    for (j = i + 1; j < n; j++) {
      tmp = 0;
      for (k = 0; k < i; k++) 
        tmp += x[i*n+k] * x[j*n+k];
      x[j*n+i] = (x[j*n+i] - tmp) / x[i*n+i];
      x[i*n+j] = 0;
    }
  }
}


/* Beaton Sweep: 1 column */
void bs(double *x, long n, long p, long k)
{
  long i, j;
  double d, b, e;

  e = exp(-8. * log(10));
  d = x[k*p+k];
  if (fabs(d) < e){
    printf("Beaton sweep: Division by %lf in %ldth column \n", d, k+1);
    exit(1);
  }
  for (j = 0; j < p; j++)
    x[j*n+k] /= d;
  for (i = 0; i < n; i++) {
    if (i == k) 
      continue;
    b = x[k*n+i];
    for (j = 0; j < p; j++)
      x[j*n+i] -= b * x[j*n+k];
    x[k*n+i] = -b / d;
  }
  x[k*n+k] = 1 / d;
}


/* Inverse via Beaton sweep */
void bs_all(double *x, long n)
{
  long j;

  for (j = 0; j < n; j++)
    bs(x ,n, n, j);
}


/* Sorting */
void sort(double *z, long n)
{
  long nodd,i,j,j1,j2,k;
  double v;

  nodd = n - 1;

  for(i=0;i<nodd;i++)
     {if(z[i] > z[i+1])
        {v = z[i + 1];
         if(z[0] >= v)k = -1;
         else 
           {j1 = 0;
            j2 = i;
            k = i / 2;
            while(z[k] != v) 
              {k = (j1 + j2) / 2;
               if(k == j1)break;
               else
                 {if(z[k] > v)j2 = k;
                  else if(z[k] < v)j1 = k; }
              } 
            }
         for(j=i;j>k;j--)z[j + 1] = z[j];
         z[k + 1] = v; 
        }
     }
}

void isort(double *z, int *ind, int n)
{
  int nodd, i, j, j1, j2, k, iv;
  double v;

  nodd = n - 1;

  for (i = 0; i < n; i++)
    ind[i] = i;
  for (i = 0; i < nodd; i++) {
    if(z[i] > z[i+1]) {
      v = z[i+1]; 
      iv = ind[i+1];
      if (z[0] >= v) 
        k = -1;
      else {
        j1 = 0;
        j2 = i;
        k = i / 2;
        while (z[k] != v) {
          k = (j1 + j2) / 2;
          if (k == j1) 
            break;
          else {
            if (z[k] > v) 
              j2 = k;
            else if (z[k] < v) 
              j1 = k; 
          }
        } 
      }
      for (j = i; j > k; j--) {
        z[j+1] = z[j]; 
        ind[j+1] = ind[j];
      }
      z[k+1] = v; 
      ind[k+1] = iv;
    }
  }
}


/* QR decomposition */
void QR(double *x, long n, long p, long s, double *r)
{
  long i, j, k;
  double e, tmp, dot(double *,double *,long,long,long);
  
  if (s > p) 
    s = p; 
  e = exp(-8.*log(10));
  for (j = 0; j < s; j++) {
    r[j*s+j] = sqrt(dot(x+j*n, x+j*n, n, 1, 1));
    if (fabs(r[j*s+j]) < e) {
      printf("QR decomposition: Division by %lf in column %ld\n", r[j*s+j], j);
      break;
    }
    for (i = 0; i < n; i++) 
      x[j*n+i] /= r[j*s+j]; 
    for (k = j+1; k < s; k++) {
      r[j*s+k] = 0;
      r[k*s+j] = dot(x+j*n, x+k*n, n, 1, 1);
      for (i = 0; i < n; i++) 
        x[k*n+i] -= x[j*n+i] * r[k*s+j];
    }
    for (k = s; k < p; k++) {
      tmp = dot(x+j*n, x+k*n, n, 1, 1);
      for (i = 0; i < n; i++) 
        x[k*n+i] -= x[j*n+i] * tmp;
    }
  }
}


/* Inverting an Upper Triangular Matrix */
void upper(double *t, long p, double *u)
{
  long i, j, k;
  double tmp, e;

  e = exp(-8.*log(10));
  for (i = 0; i < p; i++) {
    if (fabs(t[i*p+i]) < e)
      printf("Inverting upper triangular matrix: Division by %lf in column %ld\n", t[i*p+i], i);
    u[i*p+i] = 1 / t[i*p+i];
    for (j = i+1; j < p; j++) {
      tmp = 0;
      for (k = 0; k < j; k++) 
        tmp += u[k*p+i] * t[j*p+k];
      u[j*p+i] = -tmp / t[j*p+j];
      u[i*p+j] = 0;
    }
  }
} 



/* Householder transformations to tridiagonal form */
void tred2(double *a, long n, double *d, double *e)
{
  long l, k, j, i;
  double scale, hh, h, g, f;

  for (i = n-1; i >= 1; i--) {
    l = i-1;
    h = scale = 0.0;
    if (l > 0) {
      for (k = 0; k <= l; k++) 
        scale += fabs(a[k*n+i]); 
      if (scale == 0.0) 
        e[i] = a[l*n+i]; 
      else { 
        for (k = 0; k <= l; k++) { 
          a[k*n+i] /= scale; 
          h += a[k*n+i] * a[k*n+i]; 
        } 
        f = a[l*n+i];
        g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i] = scale * g;
        h -= f * g;
        a[l*n+i] = f - g;
        f = 0.0;  
        for (j = 0; j <= l; j++) {
          a[i*n+j] = a[j*n+i] / h;
	  g = 0.0;
 	  for (k = 0; k <= j; k++) 
            g += a[k*n+j] * a[k*n+i];
	  for (k = j+1; k <= l; k++) 
            g += a[j*n+k] * a[k*n+i];
	  e[j] = g / h;
	  f += e[j] * a[j*n+i];
	}
	hh = f / (h+h);
	for (j = 0; j <= l; j++) {
	  f = a[j*n+i];
	  e[j] = g = e[j] - hh*f;
	  for (k = 0; k <= j; k++) 
            a[k*n+j] -= f*e[k] + g * a[k*n+i];
	}
      }
    }
    else 
      e[i] = a[l*n+i];
    d[i] = h;
  }
  d[0] = 0.0;  
  e[0] = 0.0;
  for (i = 0; i < n; i++) {
    l = i-1;
    if (d[i]) { 
      for (j = 0; j <= l; j++) {
        g = 0.0;
	for (k = 0; k <= l; k++) 
          g += a[k*n+i] * a[j*n+k];
	for (k = 0; k <= l; k++) 
          a[j*n+k] -= g * a[i*n+k];
      }
    } 
    d[i] = a[i*n+i];
    a[i*n+i] = 1.0;
    for (j = 0; j <= l; j++) 
      a[i*n+j] = a[j*n+i] = 0.0;
  }
} 


/* Pythagoras */
double pythag(double a, double b)
{
  return sqrt(a*a + b*b);
}


/* Tridiagonal QL implicit */
void tqli(double *d, double *e, long n, double *z)
{
  long m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b, pythag(double a,double b);

  for (i = 1; i < n; i++) 
    e[i-1] = e[i];
  e[n-1] = 0.0;
  for (l = 0; l < n; l++) {
    iter = 0;
    do {
      for (m = l; m < n-1; m++) {
        dd = fabs(d[m]) + fabs(d[m+1]);
        if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
/*        if (iter++ = 30) nrerror("Too many iterations in tqli");*/
        g = (d[l+1] - d[l]) / (2.0 * e[l]);
        r = pythag(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + SIGN(r,g));
        s = c = 1.0;
        p = 0.0;
        for (i = m-1; i >= l; i--) {
          f = s * e[i];
	  b = c * e[i];
	  e[i+1] = r = pythag(f, g);
          if (r == 0.0) {
	    d[i+1] -= p;
	    e[m] = 0.0;
	    break;
  	  }
	  s = f / r;
	  c = g / r;
	  g = d[i+1] - p;
	  r = (d[i] - g) * s + 2.0 * c * b;
	  d[i+1] = g + (p = s * r);
	  g = c * r - b;
	  for (k = 0; k < n; k++) {
	    f = z[(i+1)*n+k];
	    z[(i+1)*n+k] = s * z[i*n+k] + c * f;
	    z[i*n+k] = c * z[i*n+k] - s * f;
	  }
        }
	if (r == 0.0 && i >= l) 
          continue;
	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    } 
    while (m != l)
      ;
  }
}


/* Stationary distribution from arbitrary rate matrix */
void stadis(double *sta, double *x, long n)
{
  long i, j;
  double *m, *rhs, *dmalloc(long);

  m = dmalloc((n-1)*(n-1));
  rhs = dmalloc(n-1);

  for (i = 0; i < n-1; i++) 
    rhs[i] = -x[(i+1)*n];
  for (i = 0; i < n-1; i++)
    for (j = 0; j < n-1; j++)
      m[i*(n-1)+j] = x[(i+1)*n+j+1] + rhs[i]; 
  bs_all(m, n-1);
  matmul(rhs, m, 1, n-1, n-1, sta+1);
  sta[0] = 0;
  for (i = 1; i < n; i++) 
    sta[0] -= sta[i];
  sta[0] += 1;
  free(m); 
  free(rhs);
}


/* Nondegenerate null vector */
void get_nv(double *nv, double *x, int n)
{
  long i, j, k;
  double f, vec_sum(double *, long, long);

  for (i = 0; i < n-1; i++)
    for (j = i+1; j < n; j++) {
      f = x[i+i*n] / x[j+i*n];
      x[j+i*n] = 0; 
      for (k = i+1; k < n; k++)
        x[j+k*n] = x[i+k*n] - x[j+k*n] * f; 
    }
  nv[n-1] = 1;
  for (i = n-2; i >= 0; i--) {
    nv[i] = 0;
    for (j = i+1; j < n; j++)
      nv[i] -= nv[j] * x[i+j*n]; 
    nv[i] /= x[i+i*n];
  }
  f = vec_sum(nv, n, 1);
  for (i = 0; i < n; i++)
    nv[i] /= f;
}


/* *********************************** */
/* Functions for genomics applications */
/* *********************************** */


/* Parsing annotation for genes and (exons, utrs) */
void par_anno(FILE *fl,long *gene, long *gene_n, long *exon, long *exon_n, long *utr, long *utr_n) {
  char x1[81], z[10]; 
  long p1, p2;
  
  *gene_n = 0; *exon_n = 0; *utr_n = 0;
  while (!feof(fl)) {
    fgets(x1, 81, fl);
    if (*x1 == '>' || *x1 == '<') {
      sscanf(x1+2,"%ld %ld", gene, gene+1); 
      gene += 2; 
      (*gene_n)++; 
    }
    else {
      sscanf(x1,"%ld %ld %s", &p1, &p2, z); 
      if (z[0] == 'e') {
        *exon++ = p1;
        *exon++ = p2;
        (*exon_n)++; 
      }
      if (z[0] == 'u') {
        *utr++ = p1;
        *utr++ = p2;
        (*utr_n)++; 
      }
    } 
  } 
}


/* Percentage overlap with regions */
void overlap(double *per, long *apos, long p, long *gene, long gene_n) 
{
  long i, j, left, right, found;

  for (j = 0; j < p; j++) {
    per[j] = 0;
    found = 0;
    i = 0;
    while (found == 0 && i < gene_n) {
      if (apos[6*j] <= gene[2*i+1] && apos[6*j+3] >= gene[2*i]) {
        left = apos[6*j];
        right = apos[6*j+3];
        if (apos[6*j+3] > gene[2*i+1]) 
          right = gene[2*i+1];
        if (apos[6*j] < gene[2*i]) 
          left = gene[2*i];
        per[j] = (right - left + 1.) / (apos[6*j+3] - apos[6*j] + 1); 
        found = 1;
      }
      i++;
    }
  }
}        


/* Stationary distribution from reversible rate matrix */
void qtosta(double *sta, double *x, long n)
{
  long i;
  double s = 0;

  for (i = 0; i < n; i++) {
    sta[i] = x[i*n] / x[i];
    s += sta[i];
  }
  for (i = 0; i < n; i++) 
    sta[i] /= s; 
} 
  

/* Stationary distribution from any stochastic matrix */
void ptosta(double *sta, double *p, int n)
{
  int i;
  double *tp, *dmalloc(long);

  tp = dmalloc(n*n);
  transpose(tp, p, n, n);
  for (i = 0; i < n; i++)
    tp[i+i*n] -= 1;
  get_nv(sta, tp, n);
  free(tp);
} 
  

/* Reversed transition matrix */
void reversed(double *rp, double *p, int n)
{
  int i, j;
  double *v;

  v = dmalloc(n);
  ptosta(v, p, n);
  for (i = 0; i < n; i++) 
    for (j = 0; j < n; j++) 
      rp[i+j*n] = v[j] * p[j+i*n] / v[i];
  free(v);
}


/* Simplest rate matrix */
void flat(double *x, double *sta, long n)
{
  long i, j;
  double exp_ch = 0;

  j = 0;
  for (i = 0; i < n; i++)
    if (sta[i] < 1./(10*n))
      j++; 
  if (j > 0)
    for (i = 0; i < n; i++)
      sta[i] = 1. / n;
  for (i = 0; i < n; i++) 
    exp_ch += sta[i] * (1 - sta[i]);
  exp_ch *= 100;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i == j) 
        x[i*n+j] = (sta[i] - 1) / exp_ch;
      else 
        x[i*n+j] = sta[i] / exp_ch;
    }
  }
}
      

/* Length and percent identity */
void summ(double *per, int *len, int *ct, int n, int p)
{
  int n2 = n*n, vec_isum(int *, int, int);

  for (; p-- > 0; ct += n2, per++, len++) {
    *len = vec_isum(ct, n2, 1);
    *per = (double)vec_isum(ct, n2, n+1) / (*len);
  }
}  


/* GC content */
void gccnt(double *gc1, double *gc2, double *ct, long p)
{
  double length, vec_sum(double *, long, long);

  for (; p-- > 0; ct += 16) {
    length = vec_sum(ct, 16, 1); 
    *gc1++ = (vec_sum(ct+4, 4, 1) + vec_sum(ct+8, 4, 1)) / length;
    *gc2++ = (vec_sum(ct+1, 16, 4) + vec_sum(ct+2, 16, 4)) / length;
  }
}
     

/* Composition */
void comp(double *c1, double *c2, double *length, double *ct, long n, long p)
{
  long i, n2 = n*n;
  double vec_sum(double *, long, long);

  for (; p-- > 0; ct += n2, c1 += n, c2 += n) {
    for (i = 0; i < n; i++) {
      c1[i] = vec_sum(ct+i, n2, n);
      c2[i] = vec_sum(ct+i*n, n, 1);
    }
    *length++ = vec_sum(c1, n, 1);
  }
}
     

/* Selecting alignments with GC content within a range */
long slct(double *sct, double *att, double low, double high, double *ct, long n, long p)
{
  long n2 = n*n, sp = 0;

  for (; p-- > 0; ct += n2, att++) {
    if (*att > low && *att <= high) {
      dblcp(sct, ct, n2);
      sct += n2;
      sp++; 
    }
  }
  return sp;
}


/* Collapsing counts */
void coll(int *nct, double *nt, int *np, int *ct, double *t, int n, int p)
{
  int n2 = n*n, *ind, *iptr, *dmalloci(int);
  double *tptr;

  ind = dmalloci(p);
  dblcp(nt, t, p);
  isort(nt, ind, p);
  *np = 1; 
  intcp(nct, ct+ind[0]*n2, n2);
  iptr = ind + 1; 
  for (tptr = nt+1; --p > 0; tptr++, iptr++) {
    if (fabs(tptr[0] - tptr[-1]) < .5) 
      intadd(nct, ct+(*iptr)*n2, n2);
    else {
      nct += n2; 
      intcp(nct, ct+(*iptr)*n2, n2);
      nt[*np] = *tptr;
      (*np)++;
    }
  }
  free(ind);
}


int not_gap(char a)
{
  if (a == '-')
    return 0;
  else
    return 1;
}
      

/* 0,1,2,3 to base */
char lcun(int x)
{
  if (x == 0) return 'A';
  if (x == 1) return 'C';
  if (x == 2) return 'G';
  if (x == 3) return 'T';
  return 'Z';
}

/* bases to 0,1,2,3 */
int nucl(char b)
{
  if (b == 'a' || b == 'A') return 0;
  if (b == 'c' || b == 'C') return 1;
  if (b == 'g' || b == 'G') return 2;
  if (b == 't' || b == 'T') return 3;
  return -1;
}


/* aligned bases to 0,1,...,15 */
int aln2(char a, char b)
{
  int nucl(char);

  return 4 * nucl(a) + nucl(b);
}


/* aligned bases to 0,1,...,63 */
int aln3(char a, char b, char c)
{
  int nucl(char);

  return 16 * nucl(a) + 4 * nucl(b) + nucl(c);
}

/* amino acid to 0,...,19 */
int amino(char a)
{
  if (a == 'a' || a == 'A') return 0;
  if (a == 'r' || a == 'R') return 1;
  if (a == 'n' || a == 'N') return 2;
  if (a == 'd' || a == 'D') return 3;
  if (a == 'c' || a == 'C') return 4;
  if (a == 'q' || a == 'Q') return 5;
  if (a == 'e' || a == 'E') return 6;
  if (a == 'g' || a == 'G') return 7;
  if (a == 'h' || a == 'H') return 8;
  if (a == 'i' || a == 'I') return 9;
  if (a == 'l' || a == 'L') return 10;
  if (a == 'k' || a == 'K') return 11;
  if (a == 'm' || a == 'M') return 12;
  if (a == 'f' || a == 'F') return 13;
  if (a == 'p' || a == 'P') return 14;
  if (a == 's' || a == 'S') return 15;
  if (a == 't' || a == 'T') return 16;
  if (a == 'w' || a == 'W') return 17;
  if (a == 'y' || a == 'Y') return 18;
  if (a == 'v' || a == 'V') return 19;
  return -1;
}

/* 0,...,19 to amino acid */
char onima(int x)
{
  if (x == 0) return 'A';
  if (x == 1) return 'R';
  if (x == 2) return 'N';
  if (x == 3) return 'D';
  if (x == 4) return 'C';
  if (x == 5) return 'Q';
  if (x == 6) return 'E';
  if (x == 7) return 'G';
  if (x == 8) return 'H';
  if (x == 9) return 'I';
  if (x == 10) return 'L';
  if (x == 11) return 'K';
  if (x == 12) return 'M';
  if (x == 13) return 'F';
  if (x == 14) return 'P';
  if (x == 15) return 'S';
  if (x == 16) return 'T';
  if (x == 17) return 'W';
  if (x == 18) return 'Y';
  if (x == 19) return 'V';
  return 'Z';
}


/* aligned amino acids to 0,1,...,399 */
int aminoaln2(char a, char b)
{
  int amino(char);

  return 20 * amino(a) + amino(b);
}


/* collect data from pairwise gapped alignment */
void getdata(char *x, char *y, int length, int *subct, int *gapct, int *gapsm, int res)
{
  int i, type = 0;
  int nucl(char), amino(char), aln2(char, char), aminoaln2(char, char);
  int (*rtoc)(char), (*rrtoc)(char, char);
 
  if (res == 4) {
    rtoc = &nucl;
    rrtoc = &aln2;
  }
  else {
    assert( res == 20 );
    rtoc = &amino;
    rrtoc = &aminoaln2;
  }

  for (i = 0; i < res*res; i++) subct[i] = 0;
  for (i = 0; i < 2*res; i++) gapct[i] = 0;
  gapsm[0] = gapsm[1] = 1; for (i = 2; i < 7; i++) gapsm[i] = 0;
  for (i = 0; i < length; i++) {
    if (x[i] != '-' && y[i] != '-') {
      subct[(*rrtoc)(x[i],y[i])]++;
      gapsm[2]++;      
      gapsm[3]++;
      type = 1;
    }
    if (x[i] != '-' && y[i] == '-' ) {
      gapct[(*rtoc)(x[i])]++;
      gapsm[4]++;
      type = 3;
    }
    if (x[i] == '-' && y[i] != '-' ) {
      gapct[res+(*rtoc)(y[i])]++;
      if (type == 0) gapsm[1]++; 
      if (type == 1) {
        gapsm[3]++;
        type = 2;
      }
      else if (type == 3) {
        gapsm[4]--;
        gapsm[5]++;
        gapsm[6]++;
        type = 4;
      }
      else if (type == 4) gapsm[6]++;
      else if (type == 2) gapsm[3]++;
    }
  }
}   


/* read mfa file */
void read_mfa(FILE *fl, char **name, char **seq)
{
  char *p = NULL, w = '\0';
  int i_n, i_s, ct, c;

  i_n = i_s = -1;  
  ct = 0;
  while ((c = fgetc(fl)) != EOF) {
    if (c == '>') {
      if (i_s >= 0) {
        *p = '\0';
        ct = 0;
      }
      w = 'n';
      p = name[++i_n];
    }
    if (w == 'n' && c == '\n') {
      *p = '\0';
      ct = 0;
      w = 's';
      p = seq[++i_s];
    }
    if (c != '\n') {
      if (w == 's' || ct < 15) 
        *p++ = c;
      ct++;
    }
  }
  *p = '\0';
}



/* pick a subalignment by base coordinates */
void pick_sub(char *y[], char *x[], int start, int end, int which) 
{
  int i, base_pos, align_pos, left, right;

  base_pos = align_pos = left = right = -1;
  while (*x[which] != '\0') {
    align_pos++;
    if (*x[which] != '-')
      base_pos++;
    if (base_pos == start)
      left = align_pos;
    if (base_pos == end)
      right = align_pos;
    if (start <= base_pos && base_pos <= end)
      for (i = 0; i < 3; i++)
        *(y[i]++) = *x[i];
    for (i = 0; i < 3; i++)
      x[i]++;
  }
  for (i = 0; i < 3; i++)
    *y[i] = '\0';
  if (left < 0 || right < 0) {
    printf("Not found: start = %d. end = %d\n", start, end);
  }
}


/* Diagonalizing a reversible rate matrix */
/* Input: rate matrix: x, number of rows of x: n. */
/* Output: orthogonal matrix:y, eigenvalues:d, stationary distribution:sta */
/* ydy' = sqrt(sta) * x * sqrt(1/sta) */
void diarama(double *x, double *y, double *d, double *sta, long n)
{
  long i, j;
  double *e, *dmalloc(long);

  qtosta(sta, x, n);
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      y[j*n+i] = x[j*n+i] * sqrt(sta[i]/sta[j]);
  e = dmalloc(n);
  tred2(y, n, d, e);
  tqli(d, e, n, y);
  free(e);
}   



/* Joint distribution from diagonalized rate matrix */
void dqtof(double *f, double t, double *y, double *d, double *sta, long n)
{
  long i, j, k;
  double *e, *dmalloc(long);
 
  e = dmalloc(n);
  for (i = 0; i < n; i++) 
    e[i] = exp(d[i]*t); 
  for (i = 0; i < n; i++) { 
    for (j = 0; j < n; j++) {
      f[j*n+i] = 0;
      for (k = 0; k < n; k++)
        f[j*n+i] += e[k] * y[k*n+i] * y[k*n+j];
      f[j*n+i] *= sqrt(sta[i]*sta[j]);
    }
  }
  free(e);
}   
  

/* Transition matrix from diagonalized rate matrix */
void dqtop(double *pr, double t, double *x, double *d, double *sta, long n)
{
  long i, j, k;
  double *e, *dmalloc(long);
 
  e = dmalloc(n);
  for (i = 0; i < n; i++) 
    e[i] = exp(d[i]*t); 
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      pr[n*i+j] = 0;
      for (k = 0; k < n; k++)
	pr[n*i+j] += e[k] * x[k*n+i] * x[k*n+j];
      pr[n*i+j] *= sqrt(sta[i] / sta[j]); 
    }
  }
  free(e);
}   
  

/* Loglikelihood of a triplet: pointer version */
void dqtof3(double *F, double *t, double *y, double *d, double *sta, long n)
{
  long i, n2 = n*n, n3 = n2*n;
  double *P, *p[3], *f[3], *tr[3], *s, *dmalloc(long);
 
  P = dmalloc(3*n2);
  for (i = 0; i < 3; i++) { 
    p[i] = P + i*n2;
    dqtop(p[i], t[i], y, d, sta, n);  
  }
  for (p[0] = P, f[0] = F; f[0] < F+n3; f[0] += n*n, p[0] += n)
    for (p[1] = P+n2, f[1] = f[0]; f[1] < f[0]+n2; f[1] += n, p[1] += n)
      for (p[2] = P+2*n2, f[2] = f[1]; f[2] < f[1]+n; f[2]++, p[2] += n) {
        for (i = 0; i < 3; i++)
          tr[i] = p[i];
        for (*f[2] = 0, s = sta; s < sta+n; s++)
          *f[2] += *s * *tr[0]++ * *tr[1]++ * *tr[2]++;
      }
  free(P);
}                 


/* Chi-squares goodness-of-fit diagnostic */
void chkfit(double *chi, double *x, double *t, int *ct, int n, int p)
{
  int i, n2 = n*n, vec_isum(int *, int, int);
  double *A, *y, *d, *sta, *f, E, length; 
  double *dmalloc(long);

  A = dmalloc(2*(n2+n));
  y = A; 
  d = y + n2; 
  sta = d + n; 
  f = sta + n;
  diarama(x, y, d, sta, n);    
 
  for (; p-- > 0; chi++, t++, ct += n2) {
    *chi = 0;
    length = vec_isum(ct, n2, 1);
    dqtof(f, *t, y, d, sta, n);
    for (i = 0; i < n2; i++) {
      E = length * f[i];
      *chi += (E - ct[i]) * (E - ct[i]) / E;
    } 
  }
  free(A);
}


/* Transition matrix from rate matrix */
void qtop(double *p, double t, double *x, long n)
{
  double *A, *y, *d, *sta, *dmalloc(long);

  A = dmalloc((n+2)*n);
  y = A;
  d = y + n*n;
  sta = d + n;
  diarama(x, y, d, sta, n);
  dqtop(p, t, y, d, sta, n);
  free(A);
}


/* Joint distribution from rate matrix */
void qtof(double *f, double t, double *x, long n)
{
  double *A, *y, *d, *sta, *dmalloc(long);

  A = dmalloc((n+2)*n);
  y = A;
  d = y + n*n;
  sta = d + n;
  diarama(x, y, d, sta, n);
  dqtof(f, t, y, d, sta, n);
  free(A);
}


/* Joint distribution from rate matrix */
void qtof3(double *f, double *t, double *x, long n)
{
  double *A, *y, *d, *sta, *dmalloc(long);

  A = dmalloc((n+2)*n);
  y = A;
  d = y + n*n;
  sta = d + n;
  diarama(x, y, d, sta, n);
  dqtof3(f, t, y, d, sta, n);
  free(A);
}


/* Asymptotic variance of 3 distance estimates */
void asyvar3(double *V, double *x, double *t, long N, long n)
{
  long i, j, k, n3 = n*n*n;
  double *F, *Fh, *G, u[9], th[3], h = 1e-9;
  double *dmalloc(long);

  F = dmalloc(n3);
  Fh = dmalloc(n3);
  G = dmalloc(3*n3);
  
  qtof3(F, t, x, n);
  for (i = 0; i < 3; i++) {
    dblcp(th, t, 3);
    th[i] += h;
    qtof3(Fh, th, x, n);
    for (j = 0; j < n3; j++) G[n3*i+j] = (Fh[j] - F[j]) * 1e9;
    th[i] -= h;
  }
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      V[3*j+i] = 0;
      for (k = 0; k < n3; k++) V[3*j+i] += G[n3*i+k] * G[n3*j+k] / F[k];
      V[3*j+i] *= N;
    }
  }
  cholesky(V, 3);
  upper(V, 3, u);
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      V[3*j+i] = 0;
      for (k = 0; k < 3; k++) V[3*j+i] += u[3*k+i] * u[3*k+j];
    }
  }
  free(F); free(Fh); free(G);
} 
     
  
/* Computing score from TKF model parameters and distance */
/* if indel[0] >= indel[1], then gives usual substitution score */ 
void qtos(double *s, double *x, double *indel, double t, int n)
{
  int i, j;
  double *dmalloc(long), *sta, *temp, be, r, a, g;

  temp = dmalloc(n*n);
  sta = dmalloc(n);
  qtop(temp, t, x, n); 
  qtosta(sta, x, n); 
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      temp[j*n+i] = log(temp[j*n+i] / sta[j]);
    }
  } 
  if (indel[0] < indel[1]) {
    r = exp((indel[0] - indel[1]) * t);
    be = (1 - r) / (indel[1] - indel[0] * r); 
    a = -indel[1] * t + log(indel[1] * (1 - indel[0] * be) / indel[0]);
    g = log(indel[1] * be);
    for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++) {
        s[j*(n+1)+i] = temp[j*n+i] + a;
      }
      s[j*(n+1)+n] = g;
    }
    for (i = 0; i < n; i++) s[n*(n+1)+i] = g;
    s[n*(n+1)+n] = 0;
  }
  else dblcp(s, temp, n*n);
  // delete me!
  for( int i = 0; i < n; i++ ){
    printf( "%f\t", sta[i] );
  }
  printf( "\n\n" );
  free(sta); free(temp);
}


/* Scores from a fixed TKF model */
void sc_fix(double *s, double t)
{
  double x[16], indel[2];
  x[0] = -0.009363; x[1] = 0.002101; x[2] = 0.006764; x[3] = 0.001594;
  x[4] = 0.001877; x[5] = -0.010536; x[6] = 0.001933; x[7] = 0.005980;
  x[8] = 0.005929; x[9] = 0.001897; x[10] = -0.010791; x[11] = 0.001879;
  x[12] = 0.001557; x[13] = 0.006539; x[14] = 0.002094; x[15] = -0.009453;
  indel[0] = 0.000392; indel[1] = 0.000394;

  qtos(s, x, indel, t, 4);
}


void sc_fix_gapped(double *s,double t)
{
  double x[25],indel[2];
  x[0] = -0.008593; x[1] = 0.001520; x[2] = 0.004743; x[3] = 0.001267; x[4] = 0.001062;
  x[5] = 0.001656; x[6] = -0.009451; x[7] = 0.001535; x[8] = 0.005148; x[9] = 0.001111;
  x[10] = 0.005250; x[11] = 0.001560; x[12] = -0.009542; x[13] = 0.001630; x[14] = 0.001102;
  x[15] = 0.001288; x[16] = 0.004805; x[17] = 0.001497; x[18] = -0.008690; x[19] = 0.001100;
  x[20] = 0.019469; x[21] = 0.018694; x[22] = 0.018249; x[23] = 0.019825; x[24] = -0.076237;
  indel[0] = 0.000392; indel[1] = 0.000394;

  qtos(s,x,indel,t,5);
}


void score_matrix(double **&scores,double t)
{

  double *s = new double [25];

  double x[16], indel[2];
  x[0] = -0.009363; x[1] = 0.002101; x[2] = 0.006764; x[3] = 0.001594;
  x[4] = 0.001877; x[5] = -0.010536; x[6] = 0.001933; x[7] = 0.005980;
  x[8] = 0.005929; x[9] = 0.001897; x[10] = -0.010791; x[11] = 0.001879;
  x[12] = 0.001557; x[13] = 0.006539; x[14] = 0.002094; x[15] = -0.009453;
  indel[0] = 0.000392; indel[1] = 0.000394;

  qtos(s, x, indel, t, 4);

  // delete me!
  for( int i = 0; i < 5; i++ ){
    for( int j = 0; j < 5; j++ ){
      printf( "%f\t", s[5*i + j] );
    }
    printf( "\n" );
  }
  printf( "\n\n" );

  scores = new double* [6];
  for( int i = 0; i < 4; i++ ){
    scores[i] = new double [6];
    for( int j = 0; j < 4; j++ ){
      scores[i][j] = s[5*i + j];
    }
    scores[i][5] = s[5*i + 4];
  }
  scores[4] = new double [6];
  scores[5] = new double [6];
  for( int j = 0; j < 4; j++ ){
    scores[5][j] = s[5*4 + j];
  }
  scores[5][5] = s[5*4 + 4];

  for( int i = 0; i < 4; i++ ){
    scores[i][4] = scores[i][0];
    for( int j = 1; j < 4; j++ ){
      if( scores[i][j] < scores[i][4] )
	scores[i][4] = scores[i][j];
    }
    scores[4][i] = scores[i][4];
  }
  scores[4][4] = scores[0][4];
  for( int i = 1; i < 4; i++ ){
    if( scores[i][4] < scores[4][4] )
      scores[4][4] = scores[i][4];
  }

}



void gapped_score_matrix(double **&scores,double t)
{

  double *s = new double [36];
  double x[25],indel[2];

  x[0] = -0.008593; x[1] = 0.001520; x[2] = 0.004743; x[3] = 0.001267; x[4] = 0.001062;
  x[5] = 0.001656; x[6] = -0.009451; x[7] = 0.001535; x[8] = 0.005148; x[9] = 0.001111;
  x[10] = 0.005250; x[11] = 0.001560; x[12] = -0.009542; x[13] = 0.001630; x[14] = 0.001102;
  x[15] = 0.001288; x[16] = 0.004805; x[17] = 0.001497; x[18] = -0.008690; x[19] = 0.001100;
  x[20] = 0.019469; x[21] = 0.018694; x[22] = 0.018249; x[23] = 0.019825; x[24] = -0.076237;
  indel[0] = 0.000392; indel[1] = 0.000394;

  qtos(s,x,indel,t,5);

  // delete me!
  for( int i = 0; i < 6; i++ ){
    for( int j = 0; j < 6; j++ ){
      printf( "%f\t", s[6*i + j] );
    }
    printf( "\n" );
  }
  printf( "\n\n" );

  scores = new double* [6];
  for( int i = 0; i < 4; i++ ){
    scores[i] = new double [6];
    for( int j = 0; j < 4; j++ ){
      scores[i][j] = s[6*i + j];
    }
    scores[i][5] = s[6*i + 4];
  }
  scores[4] = new double [6];
  scores[5] = new double [6];
  for( int j = 0; j < 4; j++ ){
    scores[5][j] = s[6*4 + j];
  }
  scores[5][5] = s[6*4 + 4];

  for( int i = 0; i < 4; i++ ){
    scores[4][i] = scores[i][4] = 0.0;
    for( int j = 0; j < 4; j++ ){
      if( scores[i][j] < scores[i][4] )
	scores[i][4] = scores[i][j];
      if( scores[j][i] < scores[4][i] )
	scores[4][i] = scores[j][i];
    }
  }


}



// Calculate joint probability matrix
void f_fix(double *P_t,double t)
{

  double Q[16];

  Q[0] = -0.009363; Q[1] = 0.002101; Q[2] = 0.006764; Q[3] = 0.001594;
  Q[4] = 0.001877; Q[5] = -0.010536; Q[6] = 0.001933; Q[7] = 0.005980;
  Q[8] = 0.005929; Q[9] = 0.001897; Q[10] = -0.010791; Q[11] = 0.001879;
  Q[12] = 0.001557; Q[13] = 0.006539; Q[14] = 0.002094; Q[15] = -0.009453;

  /*Q[0] = -0.0093630; Q[1] = 0.0019890; Q[2] = 0.0063465; Q[3] = 0.0015755;
  Q[4] = 0.0019890; Q[5] = -0.0105360; Q[6] = 0.0019150; Q[7] = 0.0062595;
  Q[8] = 0.0063465; Q[9] = 0.0019150; Q[10] = -0.0107910; Q[11] = 0.0019865;
  Q[12] = 0.0015755; Q[13] = 0.0062595; Q[14] = 0.0019865; Q[15] = -0.0094530;*/

  qtof(P_t,t,Q,4);

}


// Calculate transition matrix P_t = e^(Q*t) where Q is a fixed rate matrix.
void p_fix(double *P_t,double t)
{

  double Q[16];

  Q[0] = -0.009363; Q[1] = 0.002101; Q[2] = 0.006764; Q[3] = 0.001594;
  Q[4] = 0.001877; Q[5] = -0.010536; Q[6] = 0.001933; Q[7] = 0.005980;
  Q[8] = 0.005929; Q[9] = 0.001897; Q[10] = -0.010791; Q[11] = 0.001879;
  Q[12] = 0.001557; Q[13] = 0.006539; Q[14] = 0.002094; Q[15] = -0.009453;

  /*Q[0] = -0.0093630; Q[1] = 0.0019890; Q[2] = 0.0063465; Q[3] = 0.0015755;
  Q[4] = 0.0019890; Q[5] = -0.0105360; Q[6] = 0.0019150; Q[7] = 0.0062595;
  Q[8] = 0.0063465; Q[9] = 0.0019150; Q[10] = -0.0107910; Q[11] = 0.0019865;
  Q[12] = 0.0015755; Q[13] = 0.0062595; Q[14] = 0.0019865; Q[15] = -0.0094530;*/

  qtop(P_t,t,Q,4);

}




void f_fix_gapped(double *s,double t)
{
  double x[25],indel[2];
  x[0] = -0.008593; x[1] = 0.001520; x[2] = 0.004743; x[3] = 0.001267; x[4] = 0.001062;
  x[5] = 0.001656; x[6] = -0.009451; x[7] = 0.001535; x[8] = 0.005148; x[9] = 0.001111;
  x[10] = 0.005250; x[11] = 0.001560; x[12] = -0.009542; x[13] = 0.001630; x[14] = 0.001102;
  x[15] = 0.001288; x[16] = 0.004805; x[17] = 0.001497; x[18] = -0.008690; x[19] = 0.001100;
  x[20] = 0.019469; x[21] = 0.018694; x[22] = 0.018249; x[23] = 0.019825; x[24] = -0.076237;
  indel[0] = 0.000392; indel[1] = 0.000394;

  qtof(s,t,x,5);
}


void scoreMatrix( double **&mat, double t )
{

  long i,j,k;
  double *dmalloc(long);
  double *d,*sta,*y,*x,*z;
  d = dmalloc(4);
  sta = dmalloc(4);
  y = dmalloc(4*4);
  z = dmalloc(4*4);
  x = dmalloc(4*4);

  x[0] = -0.009363; x[1] = 0.002101; x[2] = 0.006764; x[3] = 0.001594;
  x[4] = 0.001877; x[5] = -0.010536; x[6] = 0.001933; x[7] = 0.005980;
  x[8] = 0.005929; x[9] = 0.001897; x[10] = -0.010791; x[11] = 0.001879;
  x[12] = 0.001557; x[13] = 0.006539; x[14] = 0.002094; x[15] = -0.009453;
  /*  x[0] = -0.009774; x[1] = 0.002835; x[2] = 0.004681; x[3] = 0.002768;
  x[4] = 0.002677; x[5] = -0.010113; x[6] = 0.002675; x[7] = 0.004426;
  x[8] = 0.004355; x[9] = 0.002635; x[10] = -0.010235; x[11] = 0.002704;
  x[12] = 0.002743; x[13] = 0.004644; x[14] = 0.002880; x[15] = -0.009899;*/
  diarama(x,z,d,sta,4);
  for (i=0;i<4;i++) {
    d[i] = exp(d[i]*t);
  }
  for (j=0;j<4;j++) {
    for (i=0;i<4;i++) {
      y[j*4+i] = 0.0;
      for (k=0;k<4;k++) {
        y[j*4+i] += d[k]*z[k*4+i]*z[k*4+j];
      }
    }
  }

  mat = new double* [4];
  for( i = 0; i < 4; i++ )
    mat[i] = new double[4];

  for (j=0;j<4;j++) {
    for (i=0;i<4;i++) {
      //mat[j][i] = log(y[j*4+i]/sqrt(sta[i]*sta[j]));
      mat[j][i] = y[j*4+i];
    }
  }
}


void transition_matrix_gapped( double **&mat, double t )
{

  long i,j,k;
  double *dmalloc(long);
  double *d,*sta,*y,*x,*z;
  d = dmalloc(5);
  sta = dmalloc(5);
  y = dmalloc(5*5);
  z = dmalloc(5*5);
  x = dmalloc(5*5);

  x[0] = -0.008593; x[1] = 0.001520; x[2] = 0.004743; x[3] = 0.001267; x[4] = 0.001062;
  x[5] = 0.001656; x[6] = -0.009451; x[7] = 0.001535; x[8] = 0.005148; x[9] = 0.001111;
  x[10] = 0.005250; x[11] = 0.001560; x[12] = -0.009542; x[13] = 0.001630; x[14] = 0.001102;
  x[15] = 0.001288; x[16] = 0.004805; x[17] = 0.001497; x[18] = -0.008690; x[19] = 0.001100;
  x[20] = 0.019469; x[21] = 0.018694; x[22] = 0.018249; x[23] = 0.019825; x[24] = -0.076237;
  diarama(x,z,d,sta,5);
  for (i=0;i<5;i++) {
    d[i] = exp(d[i]*t);
  }
  for (j=0;j<5;j++) {
    for (i=0;i<5;i++) {
      y[j*5+i] = 0.0;
      for (k=0;k<5;k++) {
        y[j*5+i] += d[k]*z[k*5+i]*z[k*5+j];
      }
    }
  }

  mat = new double* [5];
  for( i = 0; i < 5; i++ )
    mat[i] = new double[5];

  for (j=0;j<5;j++) {
    for (i=0;i<5;i++) {
      //mat[j][i] = log(y[j*4+i]/sqrt(sta[i]*sta[j]));
      mat[j][i] = y[j*5+i];
    }
  }

  free(d);
  free(sta);
  free(y);
  free(z);
  free(x);

}



/* Rate matrix from symmetric joint distribution */
void ftoq(double *x, double *f, double *t, long n)
{
  long i, j, k, n2 = n*n;
  double tot, *sta, *d, *e, *y;
  double *dmalloc(long), vec_sum(double *, long, long);
 
  sta = dmalloc(n);
  d = dmalloc(n);
  e = dmalloc(n);
  y = dmalloc(n*n);
  tot = vec_sum(f, n2, 1);
  if (fabs(tot - 1) > 0.01) {
    printf("ftoq: tot = %lf, normalizing input f\n", tot);
    for (i = 0; i < n2; i++)
      f[i] /= tot;
  }
  row_sum(sta, f, n, n); 
  for (i = 0; i < n; i++) {
    if (sta[i] < 10e-6) {
      printf("ftoq: Abnormal background frequencies\n"); 
      printm(sta, 1, n);
      *t = -1;
      free(d); free(e); free(sta); 
      return;
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      y[i*n+j] = f[i*n+j] / sqrt(sta[i] * sta[j]);
  }  
  tred2(y, n, d, e);
  tqli(d, e, n, y);
  for (i = 0; i < n; i++) 
    if (d[i] < 1e-6) {
      printf("ftoq: Negative eigenvalues\n"); 
/*      printm(d, 1, n); */
      *t = -1;
      free(d); free(e); free(y); free(sta); 
      return;
    }
  for (i = 0; i < n; i++)  
    d[i] = log(d[i]);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      x[i*n+j] = 0;
      for (k = 0; k < n; k++)
        x[i*n+j] += d[k] * y[k*n+i] * y[k*n+j];
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      x[i*n+j] *= sqrt(sta[i] / sta[j]);
  }
  k = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j != i && x[i*n+j] < 1e-6)
        k++;
    }
  }
  if (k > 0) {
    printf("ftoq: Negative rates\n"); 
/*    printm(x, n, n); */
    *t = -1;
    free(d); free(e); free(y); free(sta); 
    return;
  }  
  *t = 0;
  for (i = 0; i < n; i++) 
    *t -= sta[i] * x[i*n+i]; 
  *t *= 100;
  for (i = 0; i < n2; i++) 
    x[i] /= *t; 
  free(sta); free(d); free(e); free(y);
}


/* Loglikelihood of pair data */
double llhp(double *L, int *ct, double *x, double *t, int n, int p)
{
  int k, n2 = n*n;
  double *A, *y, *d, *sta, *f, *e, a, L0 = 0;
  double *dmalloc(long), gmin(double *, int, int *);
  double multi_llh(int *, double *, int);

  A = dmalloc(2*n2+3*n);
  y = A;
  d = A + n2;
  sta = d + n;
  f = sta + n;
  e = f + n2;
  diarama(x, y, d, sta, n);  
  for (; p-- > 0; L++, t++, ct += n2) {
    if (*t == 0) 
      for (*L = 0, k = 0; k < n; k++)
        *L += log(sta[k]) * ct[k*(n+1)];
    else {
      dqtof(f, *t, y, d, sta, n);
      a = gmin(f, n2, &k); 
      if (a <= 0) {
/*        printf("llhp: negative probabilities.\n"); */
        return 1;
      }
      else
        *L = multi_llh(ct, f, n2);
    }
    L0 += *L; 
  } 
  free(A);
  return L0;
}


/* Search routine */
double look(double *t, double *x, double *d, double *sta, int *ct, double step, int n) 
{
  int n2 = n*n, cont = 1;
  double L, newL, newt, *f, *newf;
  double *dmalloc(long), multi_llh(int *, double *, int);

  f = dmalloc(n2);
  newf = dmalloc(n2);
  dqtof(f, *t, x, d, sta, n);
  L = multi_llh(ct, f, n2);
  newt = *t + step;
  dqtof(newf, newt, x, d, sta, n);
  newL = multi_llh(ct, newf, n2);
  if (newL < L) 
    step *= -1.0;
  while (cont == 1) {
    newt = *t + step;
    if (newt <= 0) 
      cont = 0;
    else {
      dqtof(newf, newt, x, d, sta, n);
      newL = multi_llh(ct, newf, n2);
      if (newL > L) {
        *t = newt;
	L = newL;
      }
      else 
        cont = 0;
    }
  }
  free(f); free(newf);
  return L;
}  	 


/* Distance estimation */
double timest(double *t, double *L, double *x, int *ct, int n, int p)
{
  int n2 = n*n, i, vec_isum(int *, int, int);
  double *A, *y, *d, *sta, per, L0 = 0, *dmalloc(long);
  double look(double *, double *, double *, double *, int *, double, int);

  A = dmalloc(n2+2*n);
  y = A;
  d = A + n2;
  sta = d + n;
  diarama(x, y, d, sta, n);
  for (; p-- > 0; ct += n2, t++, L++) {
    *L = 0;
    per = (double)vec_isum(ct, n2, n+1) / vec_isum(ct, n2, 1);
    if (per == 1.0) {
      *t = 0; 
      for (i = 0; i < n; i++)
        *L += log(sta[i]) * ct[i*(n+1)]; 
    }
    else {
      *t = (double)floor(-100 * log(per));
      if (*t == 0)  
        *t = 1; 
      *L = look(t, y, d, sta, ct, 10., n);
      *L = look(t, y, d, sta, ct, 1., n);
    }  
    L0 += *L;
  }
  free(A);
  return L0;
}


/* Initial estimate of rate matrix */
void ini(double *x, int *ct, int n, int p)
{
  int i, n2 = n*n, *y, *dmalloci(int);
  double t, *tr, *z, *sta;
  double *dmalloc(long), vec_sum(double *, long, long);

  y = dmalloci(n2);
  tr = dmalloc(n2);
  z = dmalloc(n2);
  sta = dmalloc(n);
  row_isum(y, ct, n2, p);
  for (i = 0; i < n2; i++)
    z[i] = (double)y[i];
  transpose(tr, z, n, n);
  dbladd(z, tr, n2);
  t = vec_sum(z, n2, 1);
  for (i = 0; i < n2; i++) 
    z[i] /= t;
  ftoq(x, z, &t, n);
  if (t < 0) {
    printf("ftoq failed: t = %lf\n", t);
    row_sum(sta, z, n, n);
    flat(x, sta, n);
  }
  free(y); free(tr); free(z); free(sta);
}


/* Initial estimate of rate matrix */
/*double ini1(double *x, int *ct, int ns)
{
  int i, n2 = ns*ns, *y, *dmalloci(int);
  double t, *tr, *sta;
  double *dmalloc(long), vec_sum(double *, long, long);

  y = dmalloci(n2);
  tr = dmalloc(n2);
  sta = dmalloc(ns);
  transpose(tr, ct, ns, ns);
  dbladd(ct, tr, n2);
  t = vec_sum(ct, n2, 1);
  dbldiv(ct, t, (long)n2);
  ftoq(x, ct, &t, ns);
  if (t < 0) {
    printf("ftoq failed: t = %lf\n", t);
    row_sum(sta, ct, ns, ns);
    flat(x, sta, ns);
  }
  free(tr); free(sta);
  return t;
}
*/

/* Parameters from rate matrix */
void qtopar(double *par, double *x, long n)
{
  long i, j, k;
  double *sta, *dmalloc(long);

  sta = dmalloc(n);
  qtosta(sta, x, n);
  for (i = 0; i < n-1; i++) 
    par[i] = sta[i];
  k = n-2;
  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      k++;
      par[k] = x[j*n+i] / sta[j];
    }
  }
  free(sta);
}  


/* Rate matrix from parameters */
void partoq(double *x, double *par, long n)
{
  long i, j, k;
  double *sta, a, b = 0;
  double *dmalloc(long);

  sta = dmalloc(n); 
  sta[n-1] = 1;
  for (i = 0; i < n-1; i++) {
    sta[i] = par[i];
    sta[n-1] -= par[i];
  }
  k = n-2;
  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      k++; 
      x[i*n+j] = par[k] * sta[i];
      x[j*n+i] = par[k] * sta[j];
    }
  }
  for (i = 0; i < n; i++) {
    x[i*n+i] = 0;
    a = 0;
    for (j = 0; j < n; j++) 
      a -= x[j*n+i];
    x[i*n+i] = a;
    b -= a * sta[i];
  }
  b *= 100; 
  for (i = 0; i < n*n; i++) 
    x[i] /= b;
  free(sta);
}

       
/* Gradient */
void grad(double *g, double *x, double *t, int *ct, int n, int p)
{
  int i, npar;
  double *dmalloc(long);
  double *y, *par, *LL, L0, L1;
  double llhp(double *, int *, double *, double *, int, int);

  npar = (n-1)*(n+2)/2;
  y = dmalloc(n*n);
  par = dmalloc(npar);
  LL = dmalloc(p);
  L0 = llhp(LL,ct,x,t,n,p);
  qtopar(par,x,n);
  for (i=0;i<npar;i++) {
    par[i] += 1e-11;
    if (i>0) par[i-1] -= 1e-11;
    partoq(y,par,n);
    L1 = llhp(LL,ct,y,t,n,p);
    g[i] = (L1-L0)*1e11;
  }
  free(y); free(par); free(LL);
}


/* Rate estimation */
double ratest(double L0, double *x, double *L, double *t, int *ct, int n, int p, double s)
{
  int i, n2 = n*n, npar, cont;
  double *y, *par0, *par1, *g, *mag, *LL, a, L1;
  double *dmalloc(long);
  double gmax(double *, int, int *);
  double gmin(double *, int, int *);
  double llhp(double *, int *, double *, double *, int, int);

  npar = (n-1)*(n+2)/2;
  y = dmalloc(n2);
  par0 = dmalloc(npar);
  par1 = dmalloc(npar);
  g = dmalloc(npar);
  mag = dmalloc(npar);
  LL = dmalloc(p);

  qtopar(par0, x, n);
  grad(g, x, t, ct, n, p);
  for (i = 0; i < npar; i++) mag[i] = fabs(g[i] / par0[i]);
  a = gmax(mag, npar, &i);
  for (i = 0; i < npar; i++) g[i] /= s*a;
  dblcp(par1, par0, npar);
  dbladd(par1, g, npar);
  a = gmin(par1, npar, &i);
  if (a > 0) {
    cont = 1;
    while (cont == 1) {
      partoq(y, par1, n);
      L1 = llhp(LL, ct, y, t, n, p);
      if (L1 < L0) cont = 0;
      else {
        L0 = L1; 
        dblcp(x, y, n2);
        dbladd(par1, g, npar);
        a = gmin(par1, npar, &i);
        if (a < 0) cont = 0;
      }
    }
  }
  free(y); free(par0); free(par1); free(g); free(mag); free(LL);
  return(L0);
}


/* MLE of rate matrix and distance */
double mle2(double *x, double *t, double *L, int *ct, int n, int p, double s)
{
  int *nct, *dmalloci(int), n2 = n*n, cont = 1, climb = 0, np;
  double L0, L1, *LL, *y, *nt, *newt, diff;
  double *dmalloc(long);
  double timest(double *, double *, double *, int *, int, int);
  double ratest(double, double *, double *, double *, int *, int, int, double);
  double dist_l1(double *, double *, long);

  LL = dmalloc(p);
  nct = dmalloci(n2*p);
  y = dmalloc(n2);
  nt = dmalloc(p);
  newt = dmalloc(p);
  L0 = timest(t, L, x, ct, n, p); 
  coll(nct, nt, &np, ct, t, n, p);
  dblcp(y, x, n2);
  while (cont == 1) {
    L1 = ratest(L0, y, LL, nt, nct, n, np, s);
    if (L1 - L0 < p*0.0001) {
      if (climb == 0) cont = 0; 
      else {
        L0 = timest(newt, L, x, ct, n, p);     
        climb = 0;
        diff = dist_l1(t, newt, p);
        if (diff < 0.5) cont = 0; 
        else {
          dblcp(t, newt, p);
          coll(nct, nt, &np, ct, t, n, p);
        }
      }
    }
    else {
      climb++; 
      L0 = L1;
      dblcp(x, y, n2);
    }
  }
  free(nct); free(y); free(nt); free(newt);
  return(L0);
} 

     
/* Functions for three sequences */
/* Loglikelihood of triples */
double llh3(double *L, int *ct, double *x, double *t, int n, int p)
{
  int k, n2 = n*n, n3 = n2*n;
  double *A, a, *y, *d, *sta, *f, L0 = 0; 
  double *dmalloc(long), gmin(double *, int, int *);
  double multi_llh(int *, double *, int);
 
  A = dmalloc(n3+n2+2*n);
  y = A;
  d = y + n2;
  sta = d + n;
  f = sta + n;
  diarama(x, y, d, sta, n);
  for (; p-- > 0; ct += n3, L++, t += 3) {
    dqtof3(f, t, y, d, sta, n);
    a = gmin(f, n3, &k);
    if (a <= 0) 
      return 1;
    *L = multi_llh(ct, f, n3);
    L0 += *L;
  }
  free(A);
  return L0;
}                 


/* Output for ungapped pair analysis */
void result2(double *x, int *ct, char *name, long n, long p)
{
  char out_name[100];
  long i;
  double *t, *L, *chi, L0, mean, var;
  double *dmalloc(long);
  double mle2(double *, double *, double *, int *, int, int, double);
  FILE *out;

  t = dmalloc(p);
  L = dmalloc(p);
  chi = dmalloc(p);

  L0 = mle2(x, t, L, ct, n, p, 50.0);
  chkfit(chi, x, t, ct, n, p);
  desk(chi, p, &mean, &var);
  printf("The mean and SD of the chisquare statistics are %.1lf and %.1lf\n", mean, sqrt(var));

  strcpy(out_name, name);
  strcat(out_name, ".q");
  out = fopen(out_name, "w"); 
  fprintm(out, x, 4, 4);
  fclose(out);
 
  strcpy(out_name, name);
  strcat(out_name, ".out");
  out = fopen(out_name, "w");
  for (i = 0; i < p; i++) {
    fprintf(out, "%.1lf %.1lf %.1lf\n", t[i], L[i], chi[i]);
  }
  fclose(out);
}



/* Three-way counts to two-way counts */
void thtotw(int *ct2, int *ct3, int marg, int n, int p)
{
  int i, j, k, l, n2 = n*n, n3 = n2*n;

  if (marg == 3) {
    for (i = 0; i < p; i++) {
      for (j = 0; j < n2; j++) {
        ct2[i*n2+j] = 0;
        for (k = 0; k < n; k++) ct2[i*n2+j] += ct3[i*n3+j*n+k];
      }
    }
  }
  if (marg == 1) {
    for (i = 0; i < p; i++) {
      for (j = 0; j < n2; j++) {
        ct2[i*n2+j] = 0;
        for (k = 0; k < n; k++) ct2[i*n2+j] += ct3[i*n3+k*n2+j];
      }
    }
  }      
  if (marg == 2) {
    for (i = 0; i < p; i++) {
      for (j = 0; j < n; j++) {
        for (k = 0; k < n; k++) {
          ct2[i*n2+j*n+k] = 0; 
          for (l = 0;l < n; l++) ct2[i*n2+j*n+k] += ct3[i*n3+j*n2+l*n+k];
        }
      }
    }
  }
} 


/* Branch length estimation for triple by pairwise comparisons */
double ble(double *t, double *L, double *x, int *ct, int n, int p)
{
  int i, *ct0, *dmalloci(int);
  double *t1, *t2, *t3, L0;
  double *dmalloc(long);
  double timest(double *, double *, double *, int *, int, int);
  double llh3(double *, int *, double *, double *, int, int);

  ct0 = dmalloci(p*n*n);
  t1 = dmalloc(p);
  t2 = dmalloc(p);
  t3 = dmalloc(p);
  thtotw(ct0, ct, 1, n, p);
  L0 = timest(t1, L, x, ct0, n, p);
  thtotw(ct0, ct, 2, n, p);
  L0 = timest(t2, L, x, ct0, n, p);
  thtotw(ct0, ct, 3, n, p);
  L0 = timest(t3, L, x, ct0, n, p);

  for (i = 0; i < p; i++) {
    t[3*i+0] = (-t1[i] + t2[i] + t3[i]) / 2;
    t[3*i+1] = (-t2[i] + t3[i] + t1[i]) / 2;
    t[3*i+2] = (-t3[i] + t1[i] + t2[i]) / 2;
  }
  for (i = 0; i < 3*p; i++) {
    if (t[i] <= 0.1) t[i] = 0.1;
  }
  L0 = llh3(L, ct, x, t, n, p);
  free(ct0); free(t1); free(t2); free(t3);
  return L0;
}


/* Branch length gradient */
void gradt3(double *g, double *t, double *x, double *d, double *sta, int *ct, int n)
{
  int j, n3 = n*n*n;
  double L0, L1, *f;
  double *dmalloc(long);
  double multi_llh(int *, double *, int);

  f = dmalloc(n3);
  dqtof3(f, t, x, d, sta, n);
  L0 = multi_llh(ct, f, n3);
  for (j = 0; j < 3; j++) {
    t[j] += 1e-11;
    dqtof3(f, t, x, d, sta, n);
    L1 = multi_llh(ct, f, n3);
    g[j] = (L1 - L0) * 1e11;
    t[j] -= 1e-11;
  }
  free(f);
}

  
/* Branch length estimation for triple */
double lenest(double *t, double *L, double *x, int *ct, int n, int p)
{
  int j, cont, n2 = n*n, n3 = n2*n;
  double a, *A, g[3], mag[3], t1[3], *y, *d, *sta, *f, L0 = 0, L1;
  double *dmalloc(long), multi_llh(int *, double *, int);
  double gmax(double *, int, int *), gmin(double *, int, int *);
 
  A = dmalloc(n3+n2+2*n);
  y = A;
  d = A + n2;
  sta = d + n;
  f = sta + n;
  diarama(x, y, d, sta, n);
  for (; p-- > 0; L++, t += 3, ct += n3) {
    gradt3(g, t, y, d, sta, ct, n);
    for (j = 0; j < 3; j++) 
      mag[j] = 1 + fabs(g[j] / t[j]);
    a = gmax(mag, 3, &j);
    dblmul(g, 1/(40*a), 3); 
    dblcp(t1, t, 3);
    dbladd(t1, g, 3);
    a = gmin(t1, 3, &j);
    if (a > 0.01) {
      cont = 1;
      while (cont == 1) {
        dqtof3(f, t1, y, d, sta, n);
        L1 = multi_llh(ct, f, n3); 
        if (L1 < *L+1e-6) 
          cont = 0; 
        else {
          *L = L1;
          dblcp(t, t1, 3);
          dbladd(t1, g, 3);
          a = gmin(t1, 3, &j);
          if (a < 0.01) 
            cont = 0; 
        }
      }
    }
    L0 += *L;
  }
  free(A);
  return L0;
}  
                       

/* Rate matrix gradient */
void gradr3(double *g, double *x, double *t, int *ct, int n, int p)
{
  int i, npar;
  double *A, *y, *par, *L, L0, L1, *dmalloc(long);
  double llh3(double *, int *, double *, double *, int, int);

  npar = (n-1)*(n+2)/2;
  A = dmalloc(n*n+npar+p);
  y = A;
  par = A + n*n;
  L = par + npar;
  L0 = llh3(L, ct, x, t, n, p);
  qtopar(par, x, n);
  for (i = 0; i < npar; i++) {
    par[i] += 1e-9;
    if (i > 0) par[i-1] -= 1e-9;
    partoq(y, par, n);
    L1 = llh3(L, ct, y, t, n, p);
    g[i] = (L1 - L0) * 1e9;
  }
  free(A);
}    
  
/* Rate estimation for triple */
double ratest3(double *x, double *L, double *t, int *ct, int n, int p)
{
  int i, npar, cont;
  double *A, *y, *par0, *par1, *g, *mag, *LL, a, L0 = 0, L1, dir = 1;
  double gmax(double *, int, int *), gmin(double *, int, int *);
  double llh3(double *, int *, double *, double *, int, int); 
  double *dmalloc(long);

  npar = (n-1)*(n+2)/2;
  A = dmalloc(n*n+4*npar+p);
  y = A;
  par0 = y + n*n;
  par1 = par0 + npar;
  g = par1 + npar;
  mag = g + npar;
  LL = mag + npar;
  for (i = 0; i < p; i++) 
    L0 += L[i];
  qtopar(par0, x, n);
  gradr3(g, x, t, ct, n, p);
  for (i = 0; i < npar; i++) 
    mag[i] = fabs(g[i] / par0[i]);
  a = gmax(mag, npar, &i);
  for (i = 0; i < npar; i++) {
    g[i] /= (50.0*a); 
    par1[i] = par0[i] + g[i];
  }
  a = gmin(par1, npar, &i);
  if (a > 0.001) {
    partoq(y, par1, n);
    L1 = llh3(LL, ct, y, t, n, p);
    if (L1 < L0) { 
      dir = -1.0;    
      for (i = 0; i < npar; i++) 
        par1[i] = par0[i] + dir * g[i];
    }
  }  
  a = gmin(par1, npar, &i);
  if (a > 0.001) {
    cont = 1;
    while (cont == 1) {
      partoq(y, par1, n);
      L1 = llh3(LL, ct, y, t, n, p);
      if (L1 < L0) 
        cont = 0;
      else {
        L0 = L1;
        dblcp(x, y, n*n); 
        dblcp(L, LL, p);
        for (i = 0; i < npar; i++) 
          par1[i] += dir * g[i];
        a = gmin(par1, npar, &i);
        if (a < 0.001) 
          cont = 0;
      }
    }
  }
  free(A);
  return L0;
}      


/* MLE of rate matrix and branch lengths for triples */
double mle3(double *x, double *t, double *L, int *ct, int n, int p)
{
  long cont = 1;
  double L0, L1, L2, *newt, rgain, lgain;
  double *dmalloc(long);
  double ble(double *, double *, double *, int *, int, int);
  double ratest3(double *, double *, double *, int *, int, int);
  double lenest(double *, double *, double *, int *, int, int);

  newt = dmalloc(p);
  L0 = ble(t, L, x, ct, n, p);  
  while (cont == 1) {
    L1 = ratest3(x, L, t, ct, n, p);
    rgain = L1 - L0;
    L2 = lenest(t, L, x, ct, n, p);
    lgain = L2 - L1;
    if (rgain + lgain < 0.001*p) cont = 0;
    else L0 = L2;
  }
  free(newt);
  return(L0);
}  



/* three-way estimation */
double three_2(double *Q, double *T, double *L, int *ct, int n, int p) 
{
  int n2 = n*n, *ct0, *dmalloci(int);
  double *q, *t, *l, l0, L0;
  double *dmalloc(long);
  double llh3(double *, int *, double *, double *, int, int);
  double ble(double *, double *, double *, int *, int, int);
  double mle3(double *, double *, double *, int *, int, int);

  ct0 = dmalloci(n2*p);
  q = dmalloc(n2);
  t = dmalloc(3*p);
  l = dmalloc(p);

  thtotw(ct0, ct, 1, n, p);
  ini(q, ct0, n, p);
  l0 = ble(t, l, q, ct, n, p);
  dblcp(Q, q, n2); dblcp(T, t, 3*p); dblcp(L, l, p); L0 = l0;

  thtotw(ct0, ct, 2, n, p);
  ini(q, ct0, n, p);
  l0 = ble(t, l, q, ct, n, p);
  if (l0 > L0 || L0 > 0) {
    dblcp(Q, q, n2); dblcp(T, t, 3*p); dblcp(L, l, p); L0 = l0;
  }

  thtotw(ct0, ct, 3, n, p);
  ini(q, ct0, n, p); 
  l0 = ble(t, l, q, ct, n, p);
  if (l0 > L0 || L0 > 0) {
    dblcp(Q, q, n2); dblcp(T, t, 3*p); dblcp(L, l, p); L0 = l0;
  }

  L0 = mle3(Q, T, L, ct, n, p); 
/*  printf("max LLH = %lf\n",L0); */
  free(ct0); free(q); free(t); free(l);
  return L0;
}


void three_print(char *name, double *Q, long *apos, double *in_gene, double *in_exon, double *in_utr, double *T, double *L, long n, long p) {
  char out_name[80];
  long i, j;
  FILE *out;

  sprintf(out_name, "%s_q", name);
  out = fopen(out_name, "w");
  fprintm(out, Q, n, n);
  fclose(out);

  sprintf(out_name, "%s_cns", name);
  out = fopen(out_name, "w");
  for (i = 0; i < p; i++) {
    for (j = 6*i; j < 6*i+3; j++)
      fprintf(out, "%ld %ld ", apos[j], apos[j+3]);
    fprintf(out, "%.2lf %.2lf %.2lf ", in_gene[i], in_exon[i], in_utr[i]);
    fprintf(out, "%.2lf %.2lf %.2lf %.2lf\n", T[3*i], T[3*i+1], T[3*i+2], L[i]);
  }  
  fclose(out);
}


void three_print_noanno(char *name, double *Q, long *apos, long *ali_n, long *gap_n, double *T, double *L, long n, long p) {
  char out_name[80];
  long i, j;
  FILE *out;

  sprintf(out_name, "%s_q", name);
  out = fopen(out_name, "w");
  fprintm(out, Q, n, n);
  fclose(out);

  sprintf(out_name, "%s_cns", name);
  out = fopen(out_name, "w");
  for (i = 0; i < p; i++) {
    for (j = 6*i; j < 6*i+3; j++)
      fprintf(out, "%ld %ld ", apos[j], apos[j+3]);
    fprintf(out, "%.2lf %.2lf %.2lf ", T[3*i], T[3*i+1], T[3*i+2]);
    fprintf(out, "%.1lf %ld %ld\n", L[i]/ali_n[i], ali_n[i], gap_n[i]);
  }  
  fclose(out);
}


/* Strand-symmetric rate matrix */
void SS_rate(double *Q, double *sta, double *k) 
{
  int i;
  double tot;

  /* Q[0] */        Q[4]=sta[1];      Q[8]=sta[2]*k[0]; Q[12] = sta[3]*k[1];
  Q[1]=sta[0];      /* Q[5] */        Q[9]=sta[2]*k[2]; Q[13] = sta[3]*k[0];
  Q[2]=sta[0]*k[0]; Q[6]=sta[1]*k[2]; /* Q[10] */       Q[14] = sta[3];
  Q[3]=sta[0]*k[1]; Q[7]=sta[1]*k[0]; Q[11] = sta[2];   /* Q[15] */
  Q[0] = - Q[4] - Q[8] - Q[12];
  Q[5] = - Q[1] - Q[9] - Q[13];
  Q[10] = - Q[2] - Q[6] - Q[14];
  Q[15] = - Q[3] - Q[7] - Q[11];
  tot = - sta[0]*Q[0] - sta[1]*Q[5] - sta[2]*Q[10] - sta[3]*Q[15];
  for (i = 0; i < 16; i++)
    Q[i] /= (tot*100);
}


/* SS branch estimate */
void SS_br(double *Res, int *ct, int *base) 
{
  int i, tot, vec_isum(int *, int, int);
  double sta[4], Q[16], k[3], L, L0, T[3];
  double lenest(double *, double *, double *, int *, int, int);
  double ble(double *, double *, double *, int *, int, int);

  tot = vec_isum(base, 4, 1);
  for (i = 0; i < 4; i++) 
    sta[i] = (double)base[i]/tot;
  sta[0] = (sta[0] + sta[3]) / 2; 
  sta[3] = sta[0];
  sta[1] = (sta[1] + sta[2]) / 2; 
  sta[2] = sta[1];
  k[0] = 1; k[1] = 1; k[2] = 1;
  SS_rate(Q, sta, k);
  L0 = ble(T, &L, Q, ct, 4, 1);
  L0 = lenest(T, &L, Q, ct, 4, 1);  
  Res[12] = L0;
  dblcp(Res, T, 3);
  for (k[0] = 0.5; k[0] < 4.1; k[0] += 0.5) { 
    for (k[1] = 0.5; k[1] < 1.6; k[1] += 0.5) { 
      for (k[2] = 0.5; k[2] < 1.6; k[2] += 0.5) { 
        SS_rate(Q, sta, k);
        L0 = ble(T, &L, Q, ct, 4, 1);
        L0 = lenest(T, &L, Q, ct, 4, 1);  
        if (L0 > Res[12]) {
          Res[12] = L0;
          dblcp(Res, T, 3);
        } 
      }
    }
  }
  tot = vec_isum(ct, 64, 1); 
  asyvar3(Res+3, Q, Res, (long)tot, 4);
}



/* HKY rate matrix */
void HKY_rate(double *Q, double *sta, double k) 
{
  int i;
  double tot;

  /* Q[0] */       Q[4] = sta[1];   Q[8] = sta[2]*k; Q[12] = sta[3];
  Q[1] = sta[0];   /* Q[5] */       Q[9] = sta[2];   Q[13] = sta[3]*k;
  Q[2] = sta[0]*k; Q[6] = sta[1];   /* Q[10] */      Q[14] = sta[3];
  Q[3] = sta[0];   Q[7] = sta[1]*k; Q[11] = sta[2];  /* Q[15] */
  Q[0] = - Q[4] - Q[8] - Q[12];
  Q[5] = - Q[1] - Q[9] - Q[13];
  Q[10] = - Q[2] - Q[6] - Q[14];
  Q[15] = - Q[3] - Q[7] - Q[11];
  tot = - sta[0]*Q[0] - sta[1]*Q[5] - sta[2]*Q[10] - sta[3]*Q[15];
  for (i = 0; i < 16; i++)
    Q[i] /= (tot*100);
}


/* HKY branch estimate, fixed k */
void jc_br(double *Res, int *ct, int *base) 
{
  int i, tot, vec_isum(int *, int, int);
  double sta[4], Q[16], L;
  double lenest(double *, double *, double *, int *, int, int);
  double ble(double *, double *, double *, int *, int, int);

  for (i = 0; i < 4; i++) 
    sta[i] = 0.25;
  flat(Q, sta, 4);
  Res[12] = ble(Res, &L, Q, ct, 4, 1);  
  Res[12] = lenest(Res, &L, Q, ct, 4, 1);  
  tot = vec_isum(ct, 64, 1); 
  asyvar3(Res+3, Q, Res, (long)tot, 4);
}


/* HKY branch estimate, fixed k */
void hky_br(double *Res, int *ct, int *base) 
{
  int i, tot, vec_isum(int *, int, int);
  double sta[4], Q[16], L;
  double lenest(double *, double *, double *, int *, int, int);
  double ble(double *, double *, double *, int *, int, int);

  tot = vec_isum(base, 4, 1);
  for (i = 0; i < 4; i++) 
    sta[i] = base[i]/tot;
  sta[0] = (sta[0] + sta[3]) / 2; 
  sta[3] = sta[0];
  sta[1] = (sta[1] + sta[2]) / 2; 
  sta[2] = sta[1];
  HKY_rate(Q, sta, 2.);
  Res[12] = ble(Res, &L, Q, ct, 4, 1);  
  Res[12] = lenest(Res, &L, Q, ct, 4, 1);  
  tot = vec_isum(ct, 64, 1); 
  asyvar3(Res+3, Q, Res, (long)tot, 4);
}




/* HKY branch estimate */
void hky_br1(double *Res, int *ct, int *base) 
{
  int tot, vec_isum(int *, int, int);
  double sta[4], Q[16], L, k0, k1, k2, l0, l1, T[3];
  double llh3(double *, int *, double *, double *, int, int);
  double lenest(double *, double *, double *, int *, int, int);
  double ble(double *, double *, double *, int *, int, int);

  tot = vec_isum(base, 4, 1);
  for( int i = 0; i < 4; i++ )
    sta[i] = base[i]/tot;
  k0 = 2;
  k1 = 5;
    
  for( int i = 0; i < 5; i++ ){
    HKY_rate(Q, sta, k0);
    l0 = ble(T, &L, Q, ct, 4, 1);
    l0 = lenest(T, &L, Q, ct, 4, 1);  
    HKY_rate(Q, sta, k1);
    l1 = ble(T, &L, Q, ct, 4, 1);
    l1 = lenest(T, &L, Q, ct, 4, 1);  

    k2 = (k0 + k1) / 2;
    if (l0 > l1) 
      k1 = k2;
    else
      k0 = k2;
  }
  k2 = (k0 + k1) / 2;
  printf("k = %lf\n", k2);
  HKY_rate(Q, sta, k2);
  Res[12] = ble(T, &L, Q, ct, 4, 1);
  Res[12] = lenest(T, &L, Q, ct, 4, 1);  
  dblcp(Res, T, 3);
  tot = vec_isum(ct, 64, 1); 
  asyvar3(Res+3, Q, Res, (long)tot, 4);
}

/* HKY branch estimate */
void hky_br2(double *Res, int *ct, int *base) 
{
  int i, tot, vec_isum(int *, int, int);
  double sta[4], Q[16], L, k0, k1, k2, l0, l1, T[3];
  double llh3(double *, int *, double *, double *, int, int);
  double lenest(double *, double *, double *, int *, int, int);
  double ble(double *, double *, double *, int *, int, int);

  tot = vec_isum(base, 4, 1);
  for (i = 0; i < 4; i++) 
    sta[i] = base[i]/tot;
  sta[0] = (sta[0] + sta[3]) / 2; 
  sta[3] = sta[0];
  sta[1] = (sta[1] + sta[2]) / 2; 
  sta[2] = sta[1]; 
  k0 = 2;
  k1 = 5;
    
  for (i = 0; i < 3; i++) {
    HKY_rate(Q, sta, k0);
    l0 = ble(T, &L, Q, ct, 4, 1);
    l0 = lenest(T, &L, Q, ct, 4, 1);  
    HKY_rate(Q, sta, k1);
    l1 = ble(T, &L, Q, ct, 4, 1);
    l1 = lenest(T, &L, Q, ct, 4, 1);  

    k2 = (k0 + k1) / 2;
    if (l0 > l1) 
      k1 = k2;
    else
      k0 = k2;
  }
  k2 = (k0 + k1) / 2;
  printf("k = %lf\n", k2);
  HKY_rate(Q, sta, k2);
  Res[12] = ble(T, &L, Q, ct, 4, 1);
  Res[12] = lenest(T, &L, Q, ct, 4, 1);  
  dblcp(Res, T, 3);
  tot = vec_isum(ct, 64, 1); 
  asyvar3(Res+3, Q, Res, (long)tot, 4);
}


/* HKY branch estimate */
void HKY_br(double *Res, int *ct, int *base) 
{
  int i, tot, vec_isum(int *, int, int);
  double sta[4], Q[16], k, L, L0, T[3];
  double lenest(double *, double *, double *, int *, int, int);
  double ble(double *, double *, double *, int *, int, int);

  tot = vec_isum(base, 4, 1);
  for (i = 0; i < 4; i++) 
    sta[i] = base[i]/tot;
  sta[0] = (sta[0] + sta[3]) / 2; 
  sta[3] = sta[0];
  sta[1] = (sta[1] + sta[2]) / 2; 
  sta[2] = sta[1];
  HKY_rate(Q, sta, 1.);
  L0 = ble(T, &L, Q, ct, 4, 1);
  L0 = lenest(T, &L, Q, ct, 4, 1);  
  Res[12] = L0;
  dblcp(Res, T, 3);
  for (k = 1.5; k < 4.1; k += 0.5) { 
    HKY_rate(Q, sta, k);
    L0 = ble(T, &L, Q, ct, 4, 1);
    L0 = lenest(T, &L, Q, ct, 4, 1);  
    if (L0 > Res[12]) {
      Res[12] = L0;
      dblcp(Res, T, 3);
    }
  }
  tot = vec_isum(ct, 64, 1); 
  asyvar3(Res+3, Q, Res, (long)tot, 4);
}


/* flip coordinates */
void flip(int *x)
{
  int temp;
  if (x[0] > x[1]) {
    temp = x[1];
    x[1] = x[0];
    x[0] = temp;
  }
}

/* reorient strands */
void reorient(char *x, char *y) 
{
  if ((*x) == '-') {
    if ((*y) == '-')
      (*y) = '+';
    else
      (*y) = '-';
  }
}


/* print summary */
void pr_sm(FILE *fl, int block, int ali_n, char chr[][6], int *chrpos, char *str, double *Res, double se)
{
  int i;

  fprintf(fl, "%d\t%d\t", block, ali_n); 
  for (i = 0; i < 3; i++)
    fprintf(fl, "%s\t%d\t%d\t%c\t%.2lf\t", chr[i], chrpos[2*i], 
    chrpos[2*i+1], str[i], Res[i]); 
  fprintf(fl, "%.2lf\t%.2lf\n", se, Res[12]);
}


/* get chromosomal positions */
void get_chrpos(int *chrpos, int *coor, int *apos, char str)
{
  if (str == '+') {
    chrpos[0] = coor[0] + apos[0];
    chrpos[1] = coor[0] + apos[1];
  }
  else {
    chrpos[0] = coor[1] - apos[1];
    chrpos[1] = coor[1] - apos[0];
  }
}





/* make frequency tables */
void gather(int *ct, int *base, int *id, char *seq[], int length)
{
  int i, j, nucl(char), aln3(char, char, char);

  for (i = 0; i < 64; i++)
    ct[i] = 0;
  for (i = 0; i < 4; i++)
    base[i] = 0;
  for (i = 0; i < 3; i++)
    id[i] = 0;
  for (i = 0; i < length; i++) {
    if (seq[0][i] != '-' && seq[0][i] != 'N' && seq[1][i] != '-' &&
    seq[1][i] != 'N' && seq[2][i] != '-' && seq[2][i] != 'N') {
      ct[aln3(seq[0][i], seq[1][i], seq[2][i])]++;
      if (seq[0][i] == seq[1][i])
        id[0]++;
      if (seq[0][i] == seq[2][i])
        id[1]++;
      if (seq[1][i] == seq[2][i])
        id[2]++;
      for (j = 0; j < 3; j++)
        if (seq[j][i] != 'N' && seq[j][i] != '-')
          base[nucl(seq[j][i])]++;
    }
  }
}


void explore(char *x[], int *par)
{
  char *p[3], *seq[3];
  long left = 0, right = 0, counter = -1, hm = 0, hr = 0, mr = 0, num_cs = 0;
  long i, j, ali_n, gap_n, is_align, pos[3] = {0, 0, 0}, apos[6];
  int nucl(char), aln3(char, char, char);
  double ct[64], base[4], zero[64], tot, hmper, hrper, mrper;

  for (i = 0; i < 64; i++) 
    zero[i] = 0;
  for (i = 0; i < 4; i++) 
    base[i] = 0;
  ali_n = 0; 
  gap_n = 0;
  is_align = 0;
  for (i = 0; i < 3; i++)
    p[i] = x[i]; 
  dblcp(ct, zero, 64); 
  dblcp(base, zero, 4); 
  while (*p[0] != '\0') {
    counter++;
    for (i = 0; i < 3; i++)
      if (*p[i] != '-') 
        pos[i]++; 
    if (*p[0] != '-' && *p[1] != '-' && *p[2] != '-') {
      if (is_align == 0) {
        if (gap_n > par[2]) {
          hmper = (double)hm / ali_n;
          hrper = (double)hr / ali_n;
          mrper = (double)mr / ali_n;
          tot = vec_sum(ct, 64, 1);
          if (tot >= par[0] && hmper*100 > par[1] && hrper*100 > par[1] && 
          mrper*100 > par[1]) {
            num_cs++;
            printf("%ld ", ali_n); 
            printf("hm = %.2lf hr = %.2lf mr = %.2lf\n", hmper, hrper, mrper); 
            for (j = 0; j < 3; j++)
              seq[j] = x[j];
            prmual(stdout, left, right, seq, 3); 
          }
          dblcp(ct, zero, 64); 
          dblcp(base, zero, 4); 
          gap_n = 0;
          ali_n = 0;
          hm = 0; hr = 0; mr = 0;
        }
        else {
          gap_n = 0;
          right = counter;
        }
      }
      if (*p[0] != 'N' && *p[1] != 'N' && *p[2] != 'N') {
        ct[aln3(*p[0],*p[1],*p[2])]++;
        if (*p[0] == *p[1])
          hm++;
        if (*p[0] == *p[2])
          hr++;
        if (*p[1] == *p[2])
          mr++;
      }
      for (i = 0; i < 3; i++)  
        if (*p[i] != 'N')
          base[nucl(*p[i])]++; 
      ali_n++;
      if (ali_n == 1) {
        for (i = 0; i < 3; i++)  
          apos[2*i] = pos[i];
        left = counter;
      }
      for (i = 0; i < 3; i++)  
        apos[2*i+1] = pos[i];
      right = counter;
      is_align = 1;
    }
    else {
      gap_n++;
      is_align = 0;
    }   
    for (i = 0; i < 3; i++)
      p[i]++;
  }
  hmper = (double)hm / ali_n;
  hrper = (double)hr / ali_n;
  mrper = (double)mr / ali_n;
  if (ali_n >= par[0] && hmper*100 > par[1] && hrper*100 > par[1] && 
  mrper*100 > par[1]) {
    num_cs++;
    printf("%ld", ali_n); 
    printf("hm = %.2lf hr = %.2lf mr = %.2lf\n", hmper, hrper, mrper);
    for (j = 0; j < 3; j++)
      seq[j] = x[j];
    prmual(stdout, left, right, seq, 3); 
  }
}
 

/****** TKF functions ******/

#define N 100000000



/* collect data from pairwise gapped alignment */
void tkfdata(char *x, char *y, long *subct, long *xlen, int *gapct, int *gapsm)
{
  int i, n, type = 0;
  int nucl(char), aln2(char, char);
 
  *xlen = 0;
  for (i = 0; i < 16; i++) 
    subct[i] = 0;
  for (i = 0; i < 4; i++) 
    gapct[i] = 0;
  gapsm[0] = gapsm[1] = 1; 
  for (i = 2; i < 7; i++) 
    gapsm[i] = 0;
  n = strlen(x);
  for (i = 0; i < n; i++) {
    if (x[i] != '-' && y[i] != '-') {
      subct[aln2(x[i],y[i])]++;
      (*xlen)++;
      gapsm[2]++;      
      gapsm[3]++;
      type = 1;
    }
    if (x[i] != '-' && y[i] == '-' ) {
      gapct[nucl(x[i])]++;
      (*xlen)++;
      gapsm[4]++;
      type = 3;
    }
    if (x[i] == '-' && y[i] != '-' ) {
      gapct[nucl(y[i])]++;
      if (type == 0) 
        gapsm[1]++; 
      if (type == 1) {
        gapsm[3]++;
        type = 2;
      }
      else if (type == 3) {
        gapsm[4]--;
        gapsm[5]++;
        gapsm[6]++;
        type = 4;
      }
      else if (type == 4) 
        gapsm[6]++;
      else if (type == 2) 
        gapsm[3]++;
    }
  }
}   


/* Simulate link structure */
void simlink(int *gapsm, int n, double *indel, double t) 
{
  int *geo, i, *geometric(int, double);
  double *x, la, mu, be, r, labe, s1, s2;
  double *rdm(long);

  la = indel[0]; 
  mu = indel[1];
  r = exp((la-mu) * t);
  be = (1 - r)/(mu - la*r);
  labe = la*be; 
  s1 = exp(-mu*t); 
  s2 = s1 + mu*be;

  geo = geometric(n+1, 1-labe);
  x = rdm(n);
  
  gapsm[0] = 1; gapsm[1] = geo[n];
  for (i = 2; i < 7; i++) gapsm[i] = 0;
  for (i = 0; i < n; i++) {
    if (x[i] < s1) {
      gapsm[2]++;
      gapsm[3] += geo[i];
    }
    else if (x[i] > s2) {
      gapsm[5]++;
      gapsm[6] += geo[i];
    }
    else gapsm[4]++; 
  }
}
  
  
  
/* Loglikelihood of a pairwise alignment structure. 
Adding this to the loglikelihood of the unaligned and aligned bases gives the 
loglikelihood of the aligned sequences. */
double llhlink(double *L, int *gsm, double *indel, double *t, int p)
{

  double la, mu, be, r, labe, mube, lglabe, lgmlabe, LL = 0;
  
  la = indel[0]; 
  mu = indel[1];
  for (; p-- > 0; L++, t++, gsm += 7) {
    r = exp((la-mu)*(*t));
    be = (1 - r)/(mu - la*r);
    labe = la*be; 
    mube = mu*be; 
    lglabe = log(labe); 
    lgmlabe = log(1-labe); 
   
    *L = 0;

    /* immortal link */
    *L += gsm[0] * lgmlabe + (gsm[1] - gsm[0]) * lglabe;
  
    /* mortal links that survived */
    *L += gsm[2] * (-mu*(*t) + lgmlabe) + (gsm[3] - gsm[2]) * lglabe;
  
    /* dead mortal links with 0 descendants */
    *L += gsm[4] * log(mube);
  
    /* dead mortal links with >0 descendants */
    *L += gsm[5] * (log(1-exp(-mu*(*t))-mube) + lgmlabe) + 
    (gsm[6] - gsm[5]) * lglabe;
    LL += *L;
  } 
  return LL;
}  
  

/* Simulate indel process */
void slink(char *x1, char *x2, int *gapsm, int n, double *indel, double t) 
{
  int *geo, i, j, *geometric(int, double);
  double *x, la, mu, be, r, labe, s1, s2;
  double *rdm(long);
 
  la = indel[0]; 
  mu = indel[1];
  r = exp((la-mu) * t);
  be = (1 - r)/(mu - la*r);
  labe = la*be; 
  s1 = exp(-mu*t); 
  s2 = s1 + mu*be;

  geo = geometric(n+1, 1-labe);
  x = rdm(n);
  
  gapsm[0] = 1; 
  gapsm[1] = geo[n];
  for (j = 0; j < geo[n]-1; j++) {
    *x1++ = '-';
    *x2++ = 'A';
  }
  for (i = 2; i < 7; i++) 
    gapsm[i] = 0;
  for (i = 0; i < n; i++) {
    if (x[i] < s1) {
      *x1++ = 'A';
      *x2++ = 'A';
      gapsm[2]++;
      gapsm[3] += geo[i];
      for (j = 0; j < geo[i]-1; j++) {
        *x1++ = '-';
        *x2++ = 'A';
      }
    }
    else {
      *x1++ = 'A';
      *x2++ = '-';
      if (x[i] > s2) {
        gapsm[5]++;
        gapsm[6] += geo[i];
        for (j = 0; j < geo[i]; j++) {
          *x1++ = '-';
          *x2++ = 'A';
        }
      }
      else
        gapsm[4]++; 
    }
  }
  *x1 = '\0';
  *x2 = '\0';
}
  


/* Maximize the llh of indel, given alignment structures and distances*/
double maxlk(double *indel, double *L, int *gsm, double *t, int p)
{
  long i, j, cont = 1;
  double *L1, a[2], g[2], Nindel[2], mag, LL0, LL1;
  double llhlink(double *, int *, double *, double *, int);
  double *dmalloc(long);

  L1 = dmalloc(p);
  LL0 = llhlink(L, gsm, indel, t, p);

  for (i = 0; i < 2; i++)
    a[i] = log(indel[i]); 
  while (cont == 1) {
    for (i = 0; i < 2; i++) {
      a[i] += 1e-11;
      for (j = 0; j < 2; j++)
        Nindel[j] = exp(a[j]); 
      LL1 = llhlink(L1, gsm, Nindel, t, p);
      g[i] = (LL1-LL0) * 1e11;
      a[i] -= 1e-11;
    }
    mag = fabs(g[0]);
    if (fabs(g[1]) > mag) 
      mag = fabs(g[1]);
    for (i = 0; i < 2; i++) {
      g[i] /= 20*mag;
      a[i] += g[i];
      Nindel[i] = exp(a[i]);
    }
    if (fabs(Nindel[0]-Nindel[1]) < 1e-7 || Nindel[0] < 1e-7 || Nindel[1] < 1e-7) 
      return 1;
    LL1 = llhlink(L1, gsm, Nindel, t, p);
    while (LL1 - LL0 > 0.0001*p) {
      for (i = 0; i < p; i++) 
        L[i] = L1[i];
      LL0 = LL1;
      for (i = 0; i < 2; i++) {
        indel[i] = Nindel[i];
        a[i] += g[i];
        Nindel[i] = exp(a[i]);
      }
      if (fabs(Nindel[0]-Nindel[1]) < 1e-7 || Nindel[0] < 1e-7 || Nindel[1] < 1e-7) 
        return 1;
      LL1 = llhlink(L1, gsm, Nindel, t, p);
    }
    if (LL1 - LL0 <= 0.0001*p) cont = 0;
  }
  free(L1);
  return LL0;
}

       
double MAXlk(double *indel, double *L, int *gapsm, double *t, int p)
{
  int len;
  double L0, temp, maxlk(double *, double *, int *, double *, int);

  len = gapsm[2]+gapsm[4]+gapsm[5];
  temp = (double) (gapsm[2]+1) / (len+2);
  indel[1] = - log(temp) / *t;
  indel[0] = indel[1] / 2;
  if (indel[0] < 1e-6 || indel[1] < 1e-6 || fabs(indel[0]-indel[1]) < 1e-6) {
    printf("MAXlk: BANG!\n");
    indel[0] = 0.0001;
    indel[1] = 0.0002;
  }
  L0 = maxlk(indel, L, gapsm, t, p);
  L0 = maxlk(indel, L, gapsm, t, p);
  return L0; 
}


/* check goodness of fit for links */
void chkfit_link(double *chi, double *indel, double *t, int *gapsm, int p)
{
  int i, j, len;
  double la, mu, emu, r, be, labe, mube, Egsm[7];
  double dist_l1(double *, double *, long);

  la = indel[0]; 
  mu = indel[1];

  for (i = 0; i < p; i++) {
    r = exp((la-mu) * t[i]);
    be = (1 - r)/(mu - la*r);
    labe = la * be; 
    mube = mu * be;
    emu = exp(-mu*t[i]);
    len = gapsm[7*i+2] + gapsm[7*i+4] + gapsm[7*i+5];
    Egsm[0] = 1;
    Egsm[1] = 1/(1-labe);
    Egsm[2] = len * emu;
    Egsm[3] = Egsm[2] * Egsm[1];
    Egsm[4] = len * mube;
    Egsm[5] = len * (1 - emu - mube);
    Egsm[6] = Egsm[5] * Egsm[1];
/*    for (j = 0; j < 7; j++) gsm[j] = (double)gapsm[7*i+j]; 
    chi[i] = dist_l1(gsm, Egsm, 7);*/
    chi[i] = 0;
    for (j = 0; j < 7; j++)
      chi[i] += (gapsm[7*i+j] - Egsm[j]) * (gapsm[7*i+j] - Egsm[j]) / Egsm[j];
  }
} 

/* SE via simulation */
void indel_sim(double *stat, int sim, int *len, double *indel, double *t, int p)
{
  int i, j, *gapsm, *dmalloci(int);
  double *sim_indel, *chi, *L, *x, L0, *m, *v, *dmalloc(long);
  double MAXlk(double *, double *, int *, double *, int);

  gapsm = dmalloci(7*p);
  sim_indel = dmalloc(2*sim);
  chi = dmalloc(p);
  L = dmalloc(p);
  x = dmalloc(sim);
  m = dmalloc(sim);
  v = dmalloc(sim);
  for (i = 0; i < sim; i++) {
    for (j = 0; j < p; j++)
      simlink(gapsm+7*j, len[j], indel, t[j]);
    L0 = MAXlk(sim_indel+2*i, L, gapsm, t, p);
    chkfit_link(chi, sim_indel+2*i, t, gapsm, p);
    desk(chi, p, m+i, &L0);
    v[i] = sqrt(L0);
  }
  for (i = 0; i < sim; i++) 
    x[i] = sim_indel[2*i];
  desk(x, sim, stat, &L0);
  stat[1] = sqrt(L0);
  for (i = 0; i < sim; i++) 
    x[i] = sim_indel[2*i+1];
  desk(x, sim, stat+2, &L0);
  stat[3] = sqrt(L0);
  desk(m, sim, stat+4, &L0);
  stat[5] = sqrt(L0);
  desk(v, sim, stat+6, &L0);
  stat[7] = sqrt(L0);
  free(gapsm); free(sim_indel); free(chi); free(x); free(L); free(m); free(v);
}        


/* Likelihood of gap bases */
double lgap(int *gapct, double *x, int n, int p)
{
  int i,j,k;
  double L0 = 0.0, *sta;
  double *dmalloc(long);

  sta = dmalloc(n);
  qtosta(sta,x,n);
  for (i=0;i<n;i++) {
    k = 0;
    for (j=0;j<p;j++) {
      k += gapct[i+j*n];
    }
    L0 += (double)k * log(sta[i]);
  }
  free(sta);
  return(L0);
}      



/* Gradient */
void tkfgrad(double *g, double *x, double *t, int *subct, int *gapct, int n, int p)
{
  int i, npar;
  double *dmalloc(long);
  double *y, *par, *LL, L0, L1;
  double lgap(int *, double *, int, int);
  double llhp(double *, int *, double *, double *, int, int);

  npar = (n-1)*(n+2)/2;
  y = dmalloc(n*n);
  par = dmalloc(npar);
  LL = dmalloc(p);
  L0 = llhp(LL, subct, x, t, n, p) + lgap(gapct, x, n, 1); 
  qtopar(par,x,n);
  for (i=0;i<npar;i++) {
    par[i] += 1e-11;
    if (i>0) par[i-1] -= 1e-11;
    partoq(y,par,n);
    L1 = llhp(LL, subct, y, t, n, p) + lgap(gapct, y, n, 1); 
    g[i] = (L1-L0)*1e11;
  }
  free(y); free(par); free(LL);
}

             
/* Rate estimation in TKF*/
void tkfratest(double *x, double *L, double *t, int *subct, int *gapct, int n, int p, double s)
{
  int i, npar, cont;
  double *y, *par0, *par1, *g, *mag, *LL, a, L0;
  double *dmalloc(long);
  double gmax(double *, int, int *), gmin(double *, int, int *);
  double lgap(int *, double *, int, int);
  double llhp(double *, int *, double *, double *, int, int);

  npar = (n-1)*(n+2)/2;
  y = dmalloc(n*n);
  par0 = dmalloc(npar);
  par1 = dmalloc(npar);
  g = dmalloc(npar);
  mag = dmalloc(npar);
  LL = dmalloc(p);
  qtopar(par0,x,n);
  grad(g,x,t,subct,n,p);
  for (i=0;i<npar;i++) mag[i] = fabs(g[i]/par0[i]);
  a = gmax(mag,npar,&i);
  for (i=0;i<npar;i++) g[i] /= (s*a);
  for (i=0;i<npar;i++) par1[i] = par0[i] + g[i];
  a = gmin(par1,npar,&i);
  if (a > 0) {
    cont = 1;
    while (cont == 1) {
      partoq(y,par1,n);
      L0 = llhp(LL,subct,y,t,n,p) + lgap(gapct,y,n,1); 
      if (L0 < *L) cont = 0;
      else {
        *L = L0; 
        for (i=0;i<n*n;i++) x[i] = y[i]; 
        for (i=0;i<npar;i++) par1[i] += g[i];
        a = gmin(par1,npar,&i);
        if (a < 0) cont = 0;
      }
    }
  }
  free(y); free(par0); free(par1); free(g); free(mag); free(LL);
}


/* MLE of rate matrix and distance in TKF*/
double tkfmle2(double *x, double *t, double *L, int *subct, int *gapct, int n, int p, double s)
{
  int i, cont=1, climb=0, np, *nct, *dmalloci(int);
  double L0, L1, *y, *nt, *newt, diff;
  double *dmalloc(long);
  double timest(double *, double *, double *, int *, int, int);
  double lgap(int *, double *, int, int);

  nct = dmalloci(n*n*p);
  y = dmalloc(n*n);
  nt = dmalloc(p);
  newt = dmalloc(p);
  L0 = timest(t,L,x,subct,n,p) + lgap(gapct,x,n,1); 
  printf("Loglikelihood = %lf\n",L0);
  coll(nct,nt,&np,subct,t,n,p);
  for (i=0;i<n*n;i++) y[i] = x[i];
  printf("rate estimation...\n");
  while (cont == 1) {
    L1 = L0;
    tkfratest(y,&L1,nt,nct,gapct,n,np,s);
    printf("."); 
    if (L1 - L0 < p*0.0001) {
      printf("done.\n");
      if (climb == 0) cont = 0; 
      else {
        L0 = timest(newt,L,x,subct,n,p) + lgap(gapct,x,n,1);     
        climb = 0;
        diff = 0.0;
        for (i=0;i<p;i++) diff += fabs(t[i]-newt[i]);
        if (diff < 0.5) cont = 0; 
        else {
          printf("Loglikelihood = %lf\n",L0);
          printf("rate estimation...\n");
          for (i=0;i<p;i++) t[i] = newt[i];
          coll(nct,nt,&np,subct,t,n,p);
        }
      }
    }
    else {
      climb++; 
      L0 = L1;
      for (i=0;i<n*n;i++) x[i] = y[i];
    }
  }
  free(nct); free(y); free(nt); free(newt);
  return(L0);
} 


/* ****** */
/* BLOSUM */
/* ****** */

int blosum(double *blosum_d, int *ct, int n, int p, int per)
{
  int j, n2 = n*n, *z, *len, *dmalloci(int);
  double total, *temp, *per_i; 
  double *dmalloc(long);
  double vec_sum(double *, long, long);

  z = dmalloci(n2);
  len = dmalloci(p);
  temp = dmalloc(n2);
  per_i = dmalloc(p);
  summ(per_i, len, ct, n, p);
  for (j = 0; j < n2; j++) {
    z[j] = 0;
    blosum_d[j] = 0;
  }
  total = 0;
  for (j = 0; j < p; j++) {
    if (per_i[j] <= (double)per/100) {
      intadd(z, ct+j*n2, n2); 
      total += 1;
    }
  }
  for (j = 0; j < n2; j++)
    blosum_d[j] = (double)z[j];
  free(z); free(per_i); free(len);

  if (total == 0) {
    printf("BLOSUM failed at %d percent\n", per);
    return(0);
  }
  transpose(temp, blosum_d, n, n);
  dbladd(blosum_d, temp, n2);
  total = vec_sum(blosum_d, n2, 1);
  if (total == 0) {
    printf("BLOSUM failed at %d percent\n", per);
    free(temp);
    return(0);
  }
  dblmul(blosum_d, 1/total, n2);
  free(temp); 
  return(1);
}
    

/* ********************************** */
/* Functions for the resolvent method */
/* ********************************** */


/* Symmetrizing a rate matrix */
void sym_rate(double *x, long n)
{
  long i, j;
  double s = 0, *sta, *r;
  double *dmalloc(long);

  sta = dmalloc(n);
  r = dmalloc(n*n);
  stadis(sta,x,n);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      r[i*n+j] = x[i*n+j] / sta[i];
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      x[i*n+j] = (r[i*n+j] + r[j*n+i]) * sta[i];
    }
  }
  for (i = 0; i < n; i++) {
    x[i*n+i] = 0;
    for (j = 0; j < n; j++) {
      if (j != i) x[i*n+i] -= x[j*n+i]; 
    }
    s -= x[i*n+i] * sta[i];
  }
  s *= 100;
  for (i = 0; i < n*n; i++) {
    x[i] /= s;
  }
  free(sta); free(r);
}


/* Evaluating the interpolating integral */
/* int_t1^t2 exp(-at)( x1 + (x2-x1)/(t2-t1) * (t-t1) dt */
void resol_int(double *I, double a, double *t, double *x1, double *x2, long n)
{
  double t1, t2, e1, e2, *I_p, *x1_p, *x2_p, *z_p;
  
  t1 = *t++;
  t2 = *t;
  e1 = exp(- a * t1);
  e2 = exp(- a * t2);
  I_p = I;          /* I_p starts at first element of I */
  z_p = I + n * n;  /* will end at last element of I */
  x1_p = x1;        /* x1_p starts at first element of x1 */
  x2_p = x2;        /* x2_p starts at first element of x2 */
  while (I_p < z_p) {
    *I_p = (t2 * *x1_p - t1 * *x2_p + (*x2_p - *x1_p) / a) * (e1 - e2);
    *I_p += (*x2_p - *x1_p) * (t1 * e1 - t2 * e2);
    *I_p /= a * (t2 - t1);
    I_p++; x1_p++; x2_p++;
  }
} 


/* Converting counts to reversible transition matrix */
void cttop(double *p, int *ct, int n)
{
  int i, j, n2 = n*n;
  double *sta, *tran;
  double *dmalloc(long);

  sta = dmalloc(n);
  tran = dmalloc(n2);
  for (i = 0; i < n2; i++)
    p[i] = (double)ct[i];
  transpose(tran, p, n, n);
  dbladd(p, tran, n2);
  for (i = 0; i < n; i++) {
    sta[i] = 0;
    for (j = 0; j < n; j++) {
      sta[i] += p[j*n+i];
    }
  } 
  for (i = 0; i < n; i++) {
    if (sta[i] < 1) {
      for (j = 0; j < n; j++) p[j*n+i] = 1./n;
    }
    else {
      for (j = 0; j < n; j++) p[j*n+i] /= sta[i];
    }
  } 
  free(sta); free(tran);
}   


/* Rate estimation: fixed a */
void resol_fix(double *R, double a, double *P, double *T, long n, long p)
{
  long n2 = n*n, i, j;
  double *r;
  double *dmalloc(long);

  r = dmalloc(n2);
  for (i = 0; i < n2; i++) R[i] = 0;
  for (i = 1; i < p+1; i++) {
    resol_int(r, a, T+i-1, P+(i-1)*n2, P+i*n2, n);
    dbladd(R, r, n2);
  }
  bs_all(R, n);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      R[i*n+j] = -R[i*n+j];
    }
    R[i*n+i] += a;
  }
  sym_rate(R, n);
  free(r);
}


/* Rate estimation: optimized a */
double resol_rate(double *R, int *ct, double *t, int n, int p)
{
  int n2 = n*n, i, j; 
  double a, *P, *T, *L, LL[5], L0; 
  double *dmalloc(long);
  double llhp(double *, int *, double *, double *, int, int);
  double gmax(double *, int, int *), gmin(double *, int, int *);

  P = dmalloc(n2*(p+1));
  T = dmalloc(p+1);
  L = dmalloc(p);

  T[0] = 0; /* T = [0 t] */
  dblcp(T+1, t, p);

  /* first block of P is the identity matrix */
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      P[i*n+j] = 0;
    }
    P[i*n+i] = 1;
  }
  
  for (i = 1; i < p+1; i++)  cttop(P+i*n2, ct+(i-1)*n2, n);

  for (i = 0; i < 5; i++) {
    a = (double)(i+1)/100;
    resol_fix(R, a, P, T, n, p);
    LL[i] = llhp(L, ct, R, t, n, p); 
  }
  L0 = gmin(LL, 5, &i); 
  if (L0 > 0) a = 0;
  else {
    for (i = 0; i < 5; i++) {
      if (LL[i] > 0) LL[i] = L0 - 1;
    }
    L0 = gmax(LL, 5, &j);
    a = (double)(j+1)/100;
    resol_fix(R, a, P, T, n, p);
  }
  free(P); free(T); free(L);
  return(a);
} 


/* Resolvent estimation */
double resol_est(double *x, int *ct, int n, int p)
{
  int n2 = n*n, np, cont = 1, *nct, *dmalloci(int);
  double a1, a2, dis1, dis2, *y, *L, *t, *nt;
  double *dmalloc(long);
  double resol_rate(double *, int *, double *, int, int);
  double dist_l1(double *, double *, long);

  y = dmalloc(n2);
  L = dmalloc(p);
  t = dmalloc(p);
  nct = dmalloci(n2*p);
  nt = dmalloc(p); 
  
  ini(x, ct, n, p);
  timest(t, L, x, ct, n, p);
  coll(nct, nt, &np, ct, t, n, p); 
  a1 = resol_rate(x, nct, nt, n, np);
  if (a1 == 0) {
    free(L); free(t); free(nct); free(nt); 
    return(a1);
  }
  timest(nt, L, x, ct, n, p);
  dis1 = dist_l1(t, nt, p);
  dblcp(t, nt, p);

  do {
    coll(nct, nt, &np, ct, t, n, p);
    a2 = resol_rate(y, nct, nt, n, np);
    if (a2 == 0) {
      free(L); free(t); free(nct); free(nt);
      return(a1);
    }
    timest(nt, L, y, ct, n, p);
    dis2 = dist_l1(t, nt, p);
    if (dis2 >= dis1 ) cont = 0;
    else {
      dis1 = dis2; 
      a1 = a2;
      dblcp(x, y, n2); 
      dblcp(t, nt, p);
    }
    if (dis1 == 0) cont = 0; 
  } while (cont == 1); 
  free(y); free(L); free(t); free(nct); free(nt); 
  return(a1);
}


/* One simulation comparing BLOSUM, MLE and Resolvent */
double sim_BRM(double *result, double *F, int n_F, int *ct, int n, int p, int *per, int n_per)
{
  int BL, i, n2 = n*n;
  double t, edt, a;
  double *blosum_d, *q_mle, *t_mle, *l_mle, *q_res, *q_blo, *f, *ff;
  double *dmalloc(long);
  double resol_est(double *, int *, int, int);
  double dist_l1(double *, double *, long);
  int blosum(double *, int *, int, int, int);

  /* Resolvent estimate q_res */
  q_res = dmalloc(n2);
  printf("Resolvent estimation...\n");
  a = resol_est(q_res, ct, n, p);

  if (a == 0) {
    printf("failed.\n");
    free(q_res);
    return(a);
  }
  printf("done. a = %lf.\n", a);

  /* MLE estimate q_mle */
  q_mle = dmalloc(n2);
  t_mle = dmalloc(p);
  l_mle = dmalloc(p);
  ini(q_mle, ct, n, p);
  mle2(q_mle, t_mle, l_mle, ct, n, p, 50.0); 
  printf("MLE done.\n");
  free(t_mle); free(l_mle);
  
  /* Computing distances */
  q_blo = dmalloc(n2);
  f = dmalloc(n2);
  blosum_d = dmalloc(n2);
  for (i = 0; i < n_per; i++) {
    BL = blosum(blosum_d, ct, n, p, per[i]);
    if (BL == 1) printf("BLOSUM at %d done.\n", per[i]);
    if (BL == 0) printf("BLOSUM at %d failed.\n", per[i]);
    /* f is the blosum j.d. */
    dblcp(f, blosum_d, n2);
    /* get edt: effective divergence time */
    ftoq(q_blo, blosum_d, &t, n); 
    edt = floor(t);
    if (t - edt > 0.5) edt += 1;
    if (edt > (double)n_F) edt = (double)n_F;
    result[i] = edt;

    /* target j.d. */
    ff = F+(long)edt*n2;
    /* blosum vs true */
    result[n_per+i] = dist_l1(f, ff, n2);
    /* resolvent vs true */
    qtof(f, edt, q_res, n);
    result[2*n_per+i] = dist_l1(f, ff, n2);
    /* mle vs true */
    qtof(f, edt, q_mle, n);
    result[3*n_per+i] = dist_l1(f, ff, n2);
  }  
  free(q_mle); free(q_res); free(q_blo); free(f); free(blosum_d);
  return(a);
}     


/* Multiple simulations */
void msim_BRM(double *mean, double *sd, double *F, int n, int p, int n_F, int n_sim, int *per, int n_per)
{
  int n2 = n*n, i, j, success, *ct, *dmalloci(int);
  double a, m, v, *result, *B_dist, *R_dist, *M_dist, *t;
  double *dmalloc(long);
  double sim_BRM(double *, double *, int, int *, int, int, int *, int);

  ct = dmalloci(n2*p);
  result = dmalloc(4*n_per);
  t = dmalloc(n_per*n_sim);
  B_dist = dmalloc(n_per*n_sim);
  R_dist = dmalloc(n_per*n_sim);
  M_dist = dmalloc(n_per*n_sim);

  success = 0;
  while (success < n_sim) {
    for (j = 0; j < 10; j++) 
      multinomial(5000, ct+j*n2, F+(10*j+9)*n2, n2);
    for (j = 10; j < p; j++) 
      multinomial(5000, ct+j*n2, F+(10*j+9)*n2, n2);
    printf("Start %dth simulation\n", success+1);
    a = sim_BRM(result, F, n_F, ct, n, p, per, n_per);
    if (a > 0) {
      success++;
      i = success - 1;
      printf("%dth simulation done.\n",success);
      for (j = 0; j < n_per; j++) {
        t[j*n_sim+i] = result[j];
        B_dist[j*n_sim+i] = result[n_per+j];
        R_dist[j*n_sim+i] = result[2*n_per+j];
        M_dist[j*n_sim+i] = result[3*n_per+j];
      }
    }
  }
  for (j = 0; j < n_per; j++) {
    desk(t+j*n_sim, n_sim, &m, &v);
    mean[j*4] = m; sd[j*4] = sqrt(v);
    desk(B_dist+j*n_sim, n_sim, &m, &v);
    mean[j*4+1] = m; sd[j*4+1] = sqrt(v);
    desk(R_dist+j*n_sim, n_sim, &m, &v);
    mean[j*4+2] = m; sd[j*4+2] = sqrt(v);;
    desk(M_dist+j*n_sim, n_sim, &m, &v);
    mean[j*4+3] = m; sd[j*4+3] = sqrt(v);;
  }
  free(ct); free(result), free(B_dist); free(M_dist); free(R_dist); free(t);
}       



/* Storage */

/* Loglikelihood of pair data */
double llh2var(double *L, int *ct, double *x, double *t, int n, int p)
{
  int k, n2 = n*n, *aliptr, *ctptr;
  double *d, *e, *sta, *f, *y, a, L0 = 0;
  double *sptr, *lptr, *tptr;
  double *dmalloc(long);
  double gmin(double *, int, int *);
  double multi_llh(int *, double *, int);

  d = dmalloc(n);
  e = dmalloc(n);
  sta = dmalloc(n);
  f = dmalloc(n2);
  y = dmalloc(n2);

  diarama(x, y, d, sta, n);  
  aliptr = ct; tptr = t; 
  for (lptr = L; lptr < L+p; lptr++, tptr++, aliptr += n2) {
    *lptr = 0;
    if (*tptr == 0) {
      ctptr = aliptr; 
      for (sptr = sta; sptr < sta+n; sptr++, ctptr += n+1) {
        *lptr += log(*sptr) * (*ctptr);
      }
    }
    else {
      dqtof(f, *tptr, y, d, sta, n);
      a = gmin(f, n2, &k); 
      if (a <= 0) {
/*        printf("llhp: negative probabilities.\n"); */
        return(1);
      }
      else *lptr += multi_llh(aliptr, f, n2); 
    }
    L0 += *lptr; 
  } 
  free(d); free(e); free(sta); free(f); free(y);
  return(L0);
}

