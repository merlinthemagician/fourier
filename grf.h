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
