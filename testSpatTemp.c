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

int main(int argc, char **argv) {
  /* double zeta[N][2*N]; */
  double **zeta, **noise2D;
  double w=0.1, lambda=3, tau=1;
  int i=0, nIter=2500;
  FILE *fp;
  char *fn=malloc(sizeof(*fn)*(strlen("results/w01/spatTempColNoise0000Fourier.dat")+1));

  zeta=malloc(N*sizeof(*zeta));
  noise2D=malloc(N*sizeof(*noise2D));
  for(i=0; i<N; i++) {
    zeta[i]=malloc(2*N*sizeof(**zeta));
    noise2D[i]=malloc(2*N*sizeof(**noise2D));
  }

  initgrf(N);

  /* gsl_rng_set(r,323); */
  i=0;
  grf_initColouredNoise(zeta, w, tau, lambda,
			1, N, 0.05);
  sprintf(fn, "results/w01/spatTempColNoise%04iFourier.dat", i);
  fp=fopen(fn, "w");
  grf_printCompMatrix(fp, (const double **)zeta, N);
  fclose(fp);
  
  sprintf(fn, "results/w01/spatTempColNoise%04i.dat", i);
  grf_fourier2Noise(noise2D, (const double **)zeta, N);
  fp=fopen(fn, "w");
  /* grf_printCompMatrix(fp, (const double **)noise2D, N); */
  outReal(fp, (const double **)noise2D);
  fclose(fp);

  
  for(i=1; i<=nIter; i++) {
    grf_nextColouredNoise(zeta, w, tau, lambda,
  			  1, N, 0.05);
    sprintf(fn, "results/w01/spatTempColNoise%04iFourier.dat", i);
    fp=fopen(fn, "w");
    grf_printCompMatrix(fp, (const double **)zeta, N);
    fclose(fp);

      sprintf(fn, "results/w01/spatTempColNoise%04i.dat", i);
      grf_fourier2Noise(noise2D, (const double **)zeta, N);
      fp=fopen(fn, "w");
      /* grf_printCompMatrix(fp, (const double **)noise2D, N); */
      outReal(fp, (const double **)noise2D);
      fclose(fp);
  }
  
  return 0;  
}
