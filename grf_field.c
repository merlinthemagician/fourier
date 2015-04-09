/*
 *  grf_field.c
 *
 * Data structure for allocating multiple Gaussian Random Fields
 * 
 *
 * Ivo Siekmann, 09/04/2015
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "grf.h"

#define OUT stdout
#define ERR stderr

typedef struct  {
  double **A;
  size_t n;
  int complex;
}grf_field;

/*Allocate real grf_field of size n */
grf_field *grf_field_alloc(size_t n) {
  grf_field *out=malloc(sizeof(*out));

  out->n=n;
  out->A=grf_allocRealMatrix(n);
  out->complex=0;
  return out;
}

/*Allocate complex grf_field of size n */
grf_field *grf_field_complex_alloc(size_t n) {
  grf_field *out=malloc(sizeof(*out));

  out->n=n;
  out->A=grf_allocFourierMatrix(n);
  out->complex=1;
  
  return out;
}

/* Get random field component f_ij */
double grf_field_get(const grf_field *f, int i, int j) {
  if(f->complex) fprintf(ERR, "f is not a real random field\n"), exit(1);
  if( !(i < f->n) || !(j < f->n) ) {
    fprintf(ERR, "f(%i,%i) out of bounds %li\n", i, j, f->n);
    exit(1);
  }
  return f->A[i][j];
}

/* Get complete random field */
double **grf_field_getField(const grf_field *f) {
  return f->A;
}

/* Get real part Re(f_ij) */
double grf_field_getReal(const grf_field *f, int i, int j) {
  double *Ai;
  if(!f->complex) fprintf(ERR, "f is not a complex random field\n"), exit(1);
  if( !(i < f->n) || !(j < f->n) ) {
    fprintf(ERR, "f(%i,%i) out of bounds %li\n", i, j, f->n);
    exit(1);
  }
  Ai=f->A[i];
  return REAL(Ai,j);
}

/* Get imaginary part Im(f_ij) */
double grf_field_getImag(const grf_field *f, int i, int j) {
  double *Ai;
  if(!f->complex) fprintf(ERR, "f is not a complex random field\n"), exit(1);
  if( !(i < f->n) || !(j < f->n) ) {
    fprintf(ERR, "f(%i,%i) out of bounds %li\n", i, j, f->n);
    exit(1);
  }
  Ai=f->A[i];
  return IMAG(Ai,j);
}

/* Print grf_field */
void grf_field_print(FILE *fp, const grf_field *f) {
  if(f->complex) grf_printCompMatrix(fp, (const double **)f->A, f->n);
  else grf_printMatrix(fp, (const double **)f->A, f->n);
}
