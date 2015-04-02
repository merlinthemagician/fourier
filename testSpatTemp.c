/*
 *  testSpatTemp.c
 *
 * Test spatio-temporal coloured noise
 *
 *  
 *
 * Ivo Siekmann, 27/03/2015
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grf.h"

#define N 64

#define OUT stdout

void printDoubleM(FILE *fp, const double **m) {
  int i, j;
  for(i=0; i<N; i++) {
    for(j=0; j<N-1; j++) {
      fprintf(fp, "%f\t", m[i][j]);
    }
    fprintf(fp, "%f\n", m[i][j]);
    /* fprintf(fp, "\n"); */
  }

}

int main(int argc, char **argv) {
  /* double zeta[N][2*N]; */
  double **zeta, **noise2D;
  double w=5, lambda=3, tau=1;
  int i=0, nIter=2500, j=0, jMax=100;
  FILE *fp;
  /* const char *prefix="results/w5/"; */
  const char *prefix="results/";
  char *fn=malloc(sizeof(*fn)*(strlen(prefix)+strlen("spatTempColNoise_s0000_t0000Fourier.dat")+1));

  zeta=malloc(N*sizeof(*zeta));
  noise2D=malloc(N*sizeof(*noise2D));
  for(i=0; i<N; i++) {
    zeta[i]=malloc(2*N*sizeof(**zeta));
    noise2D[i]=malloc(N*sizeof(**noise2D));
  }

  initgrf(N);

  /* gsl_rng_set(r,323); */
  for(i=1; i<=nIter; i++) {
    j=0;
    changeSeed(i);
    grf_initColouredNoiseFourier(zeta, w, tau, lambda,
			  1, N, 0.05);
    sprintf(fn, "%sspatTempColNoise_s%04i_t%04iFourier.dat", prefix, i,j);
    fp=fopen(fn, "w");
    grf_printCompMatrix(fp, (const double **)zeta, N);
    fclose(fp);
  
    sprintf(fn, "%sspatTempColNoise_s%04i_t%04i.dat", prefix, i,j);
    grf_fourier2Noise(noise2D, (const double **)zeta, N);
    fp=fopen(fn, "w");
    /* grf_printCompMatrix(fp, (const double **)noise2D, N); */
    /* outReal(fp, (const double **)noise2D); */
    printDoubleM(fp, (const double **)noise2D);
    
    fclose(fp);

    for(j=1; j<=jMax; j++) {
      grf_nextColouredNoiseFourier(zeta, w, tau, lambda,
			    1, N, 0.05);
      sprintf(fn, "%sspatTempColNoise_s%04i_t%04iFourier.dat", prefix, i, j);
      fp=fopen(fn, "w");
      grf_printCompMatrix(fp, (const double **)zeta, N);
      fclose(fp);

      sprintf(fn, "%sspatTempColNoise_s%04i_t%04i.dat", prefix, i, j);
      grf_fourier2Noise(noise2D, (const double **)zeta, N);
      fp=fopen(fn, "w");
      /* outReal(fp, (const double **)noise2D); */
      printDoubleM(fp, (const double **)noise2D);
      fclose(fp);
    }
  }
  return 0;  
}
