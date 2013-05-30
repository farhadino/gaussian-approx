#ifndef eig_h
#define eig_h
#include "matrix.h"

extern void eig (Matrix *M, Matrix *Vec, Matrix *Val);

extern void tred2(Matrix M, int n, float d[], float e[]);

extern void tqli(float d[], float e[], int n, Matrix M);

extern float pythag(float a, float b);

extern Matrix addShiftMatrix (Matrix m);

extern Matrix subShiftMatrix (Matrix m);
extern Matrix getValMatrix (float d[], int dim);

#endif
