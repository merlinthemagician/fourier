/*
 *  grf.c
 *
 *
 * Generiert ein Gaussches Feld von Zufallszahlen mit vorgegebener
 * raeumlicher Korrelationsfunktion.
 * 
 * 
 *
 * Ivo Siekmann, 14. Juli 2008
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define OUT stdout
#define ERR stderr

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

static int n = 10;

static gsl_fft_halfcomplex_wavetable *hcwt;
static gsl_fft_real_workspace *work;
static gsl_rng * r;

static gsl_fft_complex_wavetable *cwt;
static gsl_fft_complex_workspace *cwork;
static double **mT;
static double **m;

/* Initialisiert den Zufallsgenerator */
void initgrf(int N) {
  int i;
  if((N%2) != 0) {
    fprintf(ERR, "initgrf(): Coloured noise generation only works for even grid dimension\n");
    exit(1);
  }
  n = N;
  r = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(r,4711);
  hcwt = gsl_fft_halfcomplex_wavetable_alloc(n);
  work = gsl_fft_real_workspace_alloc(n);

  cwt = gsl_fft_complex_wavetable_alloc(n);
  cwork = gsl_fft_complex_workspace_alloc(n);
  mT=malloc(n*sizeof(double*));

  for(i=0; i<n; i++) {
    mT[i] = calloc(2*n, sizeof(double));
  }

  m = malloc(n*sizeof(double*));
  for(i=0; i<n; i++) {
    m[i] = malloc(2*n*sizeof(double));
  }
}

/*
 * Aendert Seed des Zufallsgenerators
 */
void changeSeed(int seed) {
  gsl_rng_set(r, seed);
  fprintf(stdout, "Seed des Zufallsgenerators: %i\n", seed);
}

/*
 * Liefert in v den naechsten mit 1/f^\alpha korrelierten Vektor
 * ab. Zuvor muss die Raumdimension mithilfe von initgrf() festgelegt
 * worden sein.
 */
void grf1D(double alpha, double hx, double *v) {
  int i;
  double d;
#ifdef WHITE
  FILE *fp = fopen("white.dat", "w");
#endif
  d =  gsl_ran_gaussian(r, 1);
  v[0] = d;
#ifdef WHITE
  fprintf(fp, "%f %f\n", 0.0, d);
#endif  
  for(i = 1; i<n; i++) {
    d = gsl_ran_gaussian(r, 1);
#ifdef WHITE
  fprintf(fp, "%f %f\n", i*hx, d);
#endif  
    v[i] = d/sqrt(2);
    v[i] *= exp(-alpha/2* log((i+1)*hx) );
/*     fprintf(OUT, "%i %f\n", i, exp(alpha/2* log(i*hx) )); */
  }
#ifdef WHITE
  fflush(fp), fclose(fp);
#endif
  gsl_fft_halfcomplex_backward(v, 1, n, hcwt, work);  
}

/*Transponiert die Matrix v */
void transpose(double **v, int n) {
  int i, j;
  double d;

  for(i=0; i<n; i++) {
    for(j=i; j<n; j++) {
      d = v[i][j];
      v[i][j]=v[j][i];
      v[j][i] = d;
    }
  }
}


void compTrans(double **res, double**comp, int n) {
  int i, j;

  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      REAL(res[j],i) = REAL(comp[i],j);
      IMAG(res[j],i) = IMAG(comp[i],j);
    }
  }

}

/*Berechne Index*/
int refind(int i) {
  return i == 0? 0: n-i;
}

/*Generiere halbkomplexes zweidimensionales Feld*/
void genHC2D(double alpha, double hx, double **m) {
  int i, j, jmax;
  double x;
  for(i=0; i<=n/2; i++) {
    if( i% (n/2) == 0 ) jmax = n/2;
    else jmax = n-1;
    for(j=0; j<= jmax; j++) {
      if( ( (i==0) && ( (j==0) || (j==n/2)))  || 
	  ( (i==n/2) && ( (j==0) || (j==n/2)))  ) {
	REAL(m[i],j) = gsl_ran_gaussian(r, 1);
	IMAG(m[i],j) = 0;
      }
      else {
	double scl;
	x = i*hx*i*hx + j*hx*j*hx;

	/* Realteil */
	REAL(m[i],j) = gsl_ran_gaussian(r,1)/sqrt(2);
	scl= exp(-alpha/2*log(sqrt(x)));
	/* if(isnan(scl)) { */
	/*   fprintf(ERR, "genHC2D(): For x=%g, scl=nan\n", x); */
	/*   exit(1); */
	/* } */
	REAL(m[i],j) *= scl;

	/* Symmetrie */
   	REAL(m[refind(i)],refind(j)) = REAL(m[i],j);  

	/* Imaginaerteil */
	IMAG(m[i],j) = gsl_ran_gaussian(r,1)/sqrt(2);
	IMAG(m[i],j) *= scl;

	/* Symmetrie */
  	IMAG(m[refind(i)], refind(j)) = -IMAG(m[i],j);
      }
    }
  }
}

/*Generiere halbkomplexes zweidimensionales Feld*/
void genHC2D_noScaling(double **m) {
  int i, j, jmax;
  double x;
  for(i=0; i<=n/2; i++) {
    if( i% (n/2) == 0 ) jmax = n/2;
    else jmax = n-1;
    for(j=0; j<= jmax; j++) {
      if( ( (i==0) && ( (j==0) || (j==n/2)))  || 
	  ( (i==n/2) && ( (j==0) || (j==n/2)))  ) {
	REAL(m[i],j) = gsl_ran_gaussian(r, 1);
	IMAG(m[i],j) = 0;
      }
      else {
	/* Realteil */
	REAL(m[i],j) = gsl_ran_gaussian(r,1)/sqrt(2);

	/* Symmetrie */
   	REAL(m[refind(i)],refind(j)) = REAL(m[i],j);  

	/* Imaginaerteil */
	IMAG(m[i],j) = gsl_ran_gaussian(r,1)/sqrt(2);

	/* Symmetrie */
  	IMAG(m[refind(i)], refind(j)) = -IMAG(m[i],j);
      }
    }
  }
}

/*Generiere halbkomplexes zweidimensionales Feld --- 
  gleichverteilte Zufallszahlen*/
void unirfHC2D(double alpha, double hx, double **m) {
  int i, j, jmax;
  double x;
  for(i=0; i<=n/2; i++) {
    if( i% (n/2) == 0 ) jmax = n/2;
    else jmax = n-1;
    for(j=0; j<= jmax; j++) {
      if( ( (i==0) && ( (j==0) || (j==n/2)))  || 
	  ( (i==n/2) && ( (j==0) || (j==n/2)))  ) {
	REAL(m[i],j) = gsl_rng_uniform(r);
	IMAG(m[i],j) = 0;
      }
      else {
	x = i*hx*i*hx + j*hx*j*hx;

	/* Realteil */
	REAL(m[i],j) = gsl_rng_uniform(r)/sqrt(2);
	REAL(m[i],j) *= exp(-alpha/2*log(sqrt(x)));

	/* Symmetrie */
   	REAL(m[refind(i)],refind(j)) = REAL(m[i],j);  

	/* Imaginaerteil */
	IMAG(m[i],j) = gsl_rng_uniform(r)/sqrt(2);
	IMAG(m[i],j) *= exp(-alpha/2*log(sqrt(x)));

	/* Symmetrie */
  	IMAG(m[refind(i)], refind(j)) = -IMAG(m[i],j);
      }
    }
  }
}


/*Fourier-Transformation fuer quadratische Matrizen*/
void fft2D(double **m) {
  
}

/*Inverse Fourier-Transformation fuer quadratische Matrizen */
void invfft2D(double **m) {
  int i, j;

  for(i=0; i<n; i++) {
    gsl_fft_complex_backward(m[i], 1, n, cwt, cwork);
  }

  compTrans(mT, m, n);

  for(i=0; i<n; i++) {
    gsl_fft_complex_backward(mT[i], 1, n, cwt, cwork);
  }

/*   fprintf(OUT, "Transponiert:\n"); */
/*   for(i=0; i<n; i++) { */
/*     for(j=0; j<n; j++) { */
/*       fprintf(OUT, "(%0.3f,%0.3f) ",  */
/* 	      REAL(mT[i],j), */
/* 	      IMAG(mT[i],j)); */
/*     } */
/*   } */

  compTrans(m, mT, n);
}

void outReal(FILE *fp, double ** m) {
  int i, j;
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      fprintf(fp, "%f ", REAL(m[i],j));
    }
    fprintf(fp, "\n");
  }
}

/*
 * Liefert in v die naechste mit 1/f^\alpha korrelierte Matrix
 * ab. Zuvor muss die Raumdimension mithilfe von initgrf() festgelegt
 * worden sein.
 */
void grf2D(double alpha, double hx, double **v) {
  int i,j;

  genHC2D(alpha, hx, m);
  invfft2D(m);

  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      v[i][j] = REAL(m[i],j);
      /* if(isnan(v[i][j])) { */
      /* 	fprintf(ERR, "grf2D(): v[%i][%i]=nan\n", i, j); */
      /* 	exit(1); */
      /* } */
    }
  }
}

/*
 * Liefert in v die naechste mit 1/f^\alpha korrelierte Matrix
 * ab. Zuvor muss die Raumdimension mithilfe von initgrf() festgelegt
 * worden sein.
 */
void unirf2D(double alpha, double hx, double **v) {
  int i,j;

  unirfHC2D(alpha, hx, m);
  invfft2D(m);

  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      v[i][j] = REAL(m[i],j);
    }
  }
}


/*
 * Liefert in v eine mit 1/f^\alpha korrelierte Matrix. Es wird auf
 * eine Varianz von 1 normiert. Zuvor muss die Raumdimension mithilfe
 * von initgrf() festgelegt worden sein.
 */
void normgrf2D(double alpha, double hx, double **v) {
  int i,j;
  double mean=0, sqrsum=0, variance=0;

  genHC2D(alpha, hx, m);
  invfft2D(m);

  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      v[i][j] = REAL(m[i],j);
      /* if(isnan(v[i][j])) { */
      /* 	fprintf(ERR, "normgrf2D(): v[%i][%i]=nan\n", i, j); */
      /* 	exit(1); */
      /* } */
      mean += v[i][j];
      sqrsum += v[i][j]*v[i][j];
    }
  }
  mean /= n*n;
  /* fprintf(ERR, "normgrf2D(): mean=%g, sqrsum=%g\n", mean, sqrsum); */

  sqrsum /= n*n;

  variance = sqrsum - mean*mean;
  /* fprintf(ERR, "normgrf2D(): variance=%g\n", variance); */
  if(fabs(variance) > 1e-12) {
    for(i=0; i<n; i++) {
      for(j=0; j<n; j++) {
	v[i][j] /= sqrt(variance);
      }
    }
  }
}

/* Generates uncorrelated noise */
void whitegrf2D(double **v) {
  int i, j;
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      v[i][j]=gsl_ran_gaussian(r, 1);
    }
  }
}

#ifdef TEST

#define N (300)

void transferArray(double **m, const double a[N][2*N]) {
  int i, j;

  for(i=0; i<N; i++) {
    for(j=0; j<2*N; j++) {
      m[i][j] = a[i][j];
    }
  }
}

double data[N];
double white[N];
double comp[2*N];

int main(int argc, char **argv) {
  gsl_fft_real_wavetable *rewt = gsl_fft_real_wavetable_alloc(N);
  int i, j;
  FILE *fpW  = fopen("white.dat", "w"), *fpC = fopen("pink.dat", "w");

  double **data2D=malloc(N*sizeof(double *));
  double **comp2D=malloc(N*sizeof(double *));
  double sum=0, sqrsum=0, stddev;

const double fouB[N][2*N] =
{{17.0487, 0, 0.393753, 0.908175, 0.945821, 0.785646, -0.391744, 0, 
  0.945821, -0.785646, 0.393753, -0.908175}, 
 {-1.22802, -3.04739, -0.355645, -1.68993, 1.24848, 0.850071, 
  0.753967, -2.14966, 2.37035, 0.362388, -0.588195,0.106756},
 {0.344115, -1.15042, -0.146674, -0.0961442, -0.386403, 
  -0.537564, -1.38689, -1.00883, -2.83242, 0.988781, -0.750059, 
  1.60249}, 
 {1.393, 0, -0.227459, -0.310577, 2.9467, 0.319474, -3.28235, 0, 
  2.9467, -0.319474, -0.227459, 0.310577}, 
 {0.344115, 1.15042, -0.750059, -1.60249, -2.83242, -0.988781, -1.38689, 
  1.00883, -0.386403, 0.537564, -0.146674, 0.0961442}, 
 {-1.22802, 3.04739, -0.588195, -0.106756, 2.37035, -0.362388, 0.753967, 
  2.14966, 1.24848, -0.850071, -0.355645, 1.68993}};

  for(i=0; i<N; i++) {
    data2D[i] = malloc(N*sizeof(double));
    comp2D[i] = calloc(2*N,sizeof(double));
  }

  initgrf(N);

  gsl_rng_set(r,323);
  normgrf2D(1, 1, data2D);
  for(i=0; i<N; i++) {
    for(j=0; j<N; j++) {
      sum += data2D[i][j];
      sqrsum += data2D[i][j]*data2D[i][j];
      fprintf(fpC, "%f ", data2D[i][j]);
    }
    fprintf(fpC, "\n");
  }
  sum /= N*N;
  sqrsum /= N*N;
  sqrsum -= sum*sum;
  stddev=sqrt(sqrsum);
  fprintf(OUT, "Erwartungswert: mu = %f, Varianz = %f\n", sum, sqrsum);
  exit(0);

  sum=0, sqrsum=0;
  for(i=0; i<N; i++) {
    for(j=0; j<N; j++) {
      data2D[i][j]/=stddev;
      fprintf(fpC, "%f ", data2D[i][j]);
      sum += data2D[i][j];
      sqrsum += data2D[i][j]*data2D[i][j];
    }
    fprintf(fpC, "\n");
  }
  sum /= N*N;
  sqrsum /= N*N;
  sqrsum -= sum*sum;
  fprintf(OUT, "Erwartungswert: mu = %f, Varianz = %f\n", sum, sqrsum);
  exit(0);

  data[0] = gsl_ran_gaussian(r, 1);
  white[0] = gsl_ran_gaussian(r, 1);
  for(i = 1; i<N; i++) {
    white[i] = gsl_ran_gaussian(r, 1)/sqrt(2);
    data[i] = white[i]/sqrt(i);
  }

  fprintf(OUT, "Fourier-Ruecktransformation:\n");
  gsl_fft_halfcomplex_inverse(white, 1, N, hcwt, work);
  gsl_fft_halfcomplex_inverse(data, 1, N, hcwt, work);
  for(i=0; i<N; i++) {
    fprintf(fpW, "%i %f\n", i, white[i]);
    fprintf(fpC, "%i %f\n", i, data[i]);
  }
  fprintf(OUT, "\n");
  fflush(fpW), fclose(fpW);
  fflush(fpC), fclose(fpC);

  fpW=fopen("Fw.dat","w");
  fpC=fopen("Fp.dat","w");
  gsl_fft_real_transform(white, 1, N, rewt, work);
  gsl_fft_real_transform(data, 1, N, rewt, work);
  for(i=0; i<N; i++) {
    fprintf(fpW, "%i %f\n", i, white[i]);
    fprintf(fpC, "%i %f\n", i, data[i]);
  }
  fprintf(OUT, "\n");
  fflush(fpW), fclose(fpW);
  fflush(fpC), fclose(fpC);

  return 0;  
}
#endif
