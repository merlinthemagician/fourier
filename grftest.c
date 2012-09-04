/*
 * grftest.c
 *
 *
 * Testet die Erzeugung korrelierter gaussscher Zufallsfelder mithilfe
 * des lgorithmus von Lang (2007).
 * 
 * Aenderung der Frequenzen mithilfe der diskreten Fouriertransformation
 * 
 * 
 *
 * Ivo Siekmann, 23.10.2008
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include "grf.h"

void mat2File(const char *fn, const double **v, int M, int N) {
  int i, j;

  FILE *f = fopen(fn, "w");

  for(i=0; i<M; i++) {
    for(j=0; j<N; j++) {
      fprintf(f, "%g\t", v[i][j]);
    }
    fprintf(f, "\n");
  }
  fflush(f), fclose(f);
}

void mat2LengthFile(const char *fn, const double **v, 
		    double hx, double hy, int M, int N) {
  int i, j;

  FILE *f = fopen(fn, "w");

  for(i=0; i<M; i++) {
    for(j=0; j<N; j++) {
      fprintf(f,"%8.3g %8.3g %g\n", i*hx, j*hy, v[i][j]);
    }
    fprintf(f, "\n");
  }
  fflush(f), fclose(f);
}

void mat2sgn(double **m, int M, int N) {
  int i, j;
  for(i=0; i<M; i++) {
    for(j=0; j<N; j++) {
      m[i][j] = m[i][j] > 0 ?1:0;
    }
  }
}

#define N 100

int main(int argc, char **argv) {
  double seed = 57; 
  double ** noise = malloc(N*sizeof(double*));
  gsl_rng *r;
  int i, j;

  for(i=0; i<N; i++) {
    noise[i]=calloc(N,sizeof(double));
  }

  r = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(r,seed);

  for(i=0; i<N; i++) {
    for(j=0; j<N; j++) {
      noise[i][j] = gsl_ran_gaussian(r,1);
    }
  }

/*   mat2File("whiteNoise.dat", noise, N, N); */
  mat2LengthFile("whiteNoise.dat", noise, 5, 5, N, N);
/*   mat2sgn(noise, N, N); */
/*   mat2File("whiteNoiseSgn.dat", noise, N, N); */
/*   mat2LengthFile("whiteNoise.dat", noise, 5, 5, N, N); */

  /* Initialisiert den Zufallsgenerator */
  initgrf(N);

  /* Aendert Seed des Zufallsgenerators */
  changeSeed(seed);

  /* Rosa Rauschen */
  normgrf2D(1, 1, noise);
/*   mat2File("pinkNoise.dat", noise, N, N); */
  mat2LengthFile("pinkNoise.dat", noise, 5, 5, N, N);
/*   mat2sgn(noise, N, N); */
/*   mat2File("pinkNoiseSgn.dat", noise, N, N); */
/*   mat2LengthFile("pinkNoise.dat", noise, 5, 5, N, N); */

  normgrf2D(2,1,noise);
/*   mat2File("redNoise.dat", noise, N, N); */
  mat2LengthFile("redNoise.dat", noise, 5, 5, N, N);
/*   mat2sgn(noise, N, N); */
/*   mat2File("redNoiseSgn.dat", noise, N, N); */
/*   mat2LengthFile("redNoise.dat", noise, 5, 5, N, N); */

  return 0;  
}
 
