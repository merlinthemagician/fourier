/*
 *  grf.h
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

/* Initialisiert den Zufallsgenerator */
void initgrf(int N);

/*
 * Aendert Seed des Zufallsgenerators
 */
void changeSeed(int seed);

/*
 * Liefert in v den naechsten mit 1/f^\alpha korrelierten Vektor
 * ab. Zuvor muss die Raumdimension mithilfe von initgrf() festgelegt
 * worden sein.
 */
void grf1D(double alpha, double hx, double *v);

/*
 * Liefert in v die naechste mit 1/f^\alpha korrelierte Matrix
 * ab. Zuvor muss die Raumdimension mithilfe von initgrf() festgelegt
 * worden sein.
 */
void grf2D(double alpha, double hx, double **v);

/*
 * Liefert in v eine mit 1/f^\alpha korrelierte Matrix. Es wird auf
 * eine Varianz von 1 normiert. Zuvor muss die Raumdimension mithilfe
 * von initgrf() festgelegt worden sein.
 */
void normgrf2D(double alpha, double hx, double **v);

/*
 * Liefert in v die naechste mit 1/f^\alpha korrelierte Matrix
 * ab. Zuvor muss die Raumdimension mithilfe von initgrf() festgelegt
 * worden sein.
 */
void unirf2D(double alpha, double hx, double **v);

/* Generates uncorrelated noise */
void whitegrf2D(double **v);

/* Initial condition for spatially coloured noise in Fourier space with:
 - Intensity: w
 - temporal correlation length: tau
 - spatial correlation length: lambda

 Spatial discretisation: hx, n
 Temporal discretisation: ht

 according to Garcia-Ojalvo, Sanchez, (1992)
 */
void grf_initColouredNoise(double **zeta, double w, double tau, double lambda,
			   double hx, int n, double ht);

/* Next time step of spatially coloured noise in Fourier space with:
 - Intensity: w
 - temporal correlation length: tau
 - spatial correlation length: lambda

 Spatial discretisation: hx, n
 Temporal discretisation: ht

 according to Garcia-Ojalvo, Sanchez, (1992)
 */
void grf_nextColouredNoise(double **zeta, double w, double tau, double lambda,
			   double hx, int n, double ht);

/* Print complex matrix zeta to file fp */
void grf_printCompMatrix(FILE *fp, const double **zeta, int n);

/* Convert Fourier-transformed noise back to real */
void grf_fourier2Noise(double **m, const double **zeta, int n);

void outReal(FILE *fp, const double ** m);
