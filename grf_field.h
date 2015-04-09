#ifndef GRFFIELD
#define GRFFIELD
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

#include "grf.h"

typedef struct grf_field {
  double **A;
  size_t n;
}grf_field;

/*Allocate grf_field of size n */
grf_field *grf_field_alloc(size_t n);

/*Allocate complex grf_field of size n */
grf_field *grf_field_complex_alloc(size_t n);

/* Get random field component f_ij */
double grf_field_get(const grf_field *f, int i, int j);

/* Get complete random field */
double ** grf_field_getField(const grf_field *f);

/* Get real part Re(f_ij) */
double grf_field_getReal(const grf_field *f, int i, int j);

/* Get imaginary part Im(f_ij) */
double grf_field_getImag(const grf_field *f, int i, int j);

/* Print grf_field */
void grf_field_print(FILE *fp, const grf_field *f);

#endif
