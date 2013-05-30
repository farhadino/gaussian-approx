#ifndef matrix_h
#define matrix_h

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

typedef float Element;

typedef Element *Line;
typedef Line *Lines;

typedef struct
{
  Lines lines;
  Element *elements;
  int width;
  int height;
} *Matrix;

typedef struct
{
  Lines lines;
  Element *elements;
  int width;
  int height;
} *Matrix2;

// constructors

extern Matrix newMatrix (int h, int w);
extern Matrix zeroMatrix (int h, int w);
extern Matrix unitMatrix (int h, int w);
extern Matrix onesMatrix (int h, int w);

// destructor

extern void freeMatrix (Matrix m);

// matrix operations

extern Matrix sizeOfMatrix (Matrix m);
extern Matrix addMatrix (Matrix m1, Matrix m2);
extern Matrix subMatrix (Matrix m1, Matrix m2);
extern Matrix mulMatrix (Matrix m1, Matrix m2);
extern Matrix transposeMatrix (Matrix m);
extern Matrix choleskyMatrix (Matrix m);
extern Matrix invertCovMatrix (Matrix m);
extern Matrix mulScalarMatrix (float a, Matrix m);
extern Matrix sumMatrix (Matrix m, int dim);
extern Matrix appendMatrix (Matrix m1, int x1, int x2, int y1, int y2, Matrix m2, int q1, int q2, int r1, int r2);
extern void appMatrix (Matrix m1, int x1, int x2, int y1, int y2, Matrix m2, int q1, int q2, int r1, int r2);


// matrix access operations

extern Element elem (Matrix m, int i, int j);
extern void setElem (Matrix m, int i, int j, Element v);

// print whole matrix

extern void printMatrix (Matrix m);

#endif
