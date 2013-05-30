#include "matrix.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tracker.h"

/*----------------------------------------------------------------------------*/

Matrix
newMatrix (int h, int w)
{
  Matrix res = malloc (sizeof (*res));

  if (res)
    {
      res->lines = malloc (h * sizeof (Line));
      res->elements = malloc (h * w * sizeof (Element));
      res->width = w;
      res->height = h;

    {
	Lines lns = res->lines;
	Line l = res->elements;

	if (lns && l)
	  {
	    while (h--)
	      {
		        *lns++ = l;
		         l += w;
	      }

	    return res;
	  }
    }
  }

  // heap overflow
  perror ("newMatrix(): cant allocate matrix");
  exit (1);
}

/*----------------------------------------------------------------------------*/

void
freeMatrix (Matrix m)
{
  free (m->elements);
  free (m->lines);
  free (m);
}

/*----------------------------------------------------------------------------*/

Matrix
zeroMatrix (int h, int w)
{
  Matrix res = newMatrix (h, w);
  Lines l = res->lines;

  int i,j;
  for (i = 0; i < h; ++i)
    {
      for (j = 0; j < w; ++j)
      {
        l[i][j] = 0.0;
      }
    }

  return res;
}

/*----------------------------------------------------------------------------*/

Matrix
unitMatrix (int h, int w)
{
  Matrix res = zeroMatrix (h, w);
  Lines r = res->lines;

  int i;
  for (i = 0; i < w && i < h; ++i)
    {
      r[i][i] = 1.0;
    }

  return res;
}

/*----------------------------------------------------------------------------*/

Matrix
onesMatrix (int h, int w)
{
  Matrix res = newMatrix (h, w);
  Lines r = res->lines;

  int i,j;
  for (i = 0; i < h; ++i)
    {
      for (j = 0; j < w; ++j)
      {
        r[i][j] = 1.0;
      }
    }
  return res;
}

/*----------------------------------------------------------------------------*/

Matrix
sizeOfMatrix (Matrix m)
{
  Matrix res = newMatrix (1, 2);
  Lines r = res->lines;
  r[0][0] = m->height;
  r[0][1] = m->width;  
  
  return res;
}

/*----------------------------------------------------------------------------*/

Matrix
addMatrix (Matrix m1, Matrix m2)
{
  int h = m1->height;
  int w = m1->width;

  assert (w == m2->width && h == m2->height);

    Matrix res = newMatrix (h, w);
    Lines r = res->lines;
    Lines r1 = m1->lines;
    Lines r2 = m2->lines;

    int j;
    for (j = 0; j < h; ++j)
    {
        int i;
        for (i = 0; i < w; ++i)
        {
            r[j][i] = r1[j][i] + r2[j][i];
        }
    }
    return res;
}

/*----------------------------------------------------------------------------*/

Matrix
subMatrix (Matrix m1, Matrix m2)
{
  int h = m1->height;
  int w = m1->width;

  assert (w == m2->width && h == m2->height);

  {
    Matrix res = newMatrix (h, w);
    Lines r = res->lines;
    Lines r1 = m1->lines;
    Lines r2 = m2->lines;

    int j;
    for (j = 0; j < h; ++j)
    {
        int i;
        for (i = 0; i < w; ++i)
        {
            r[j][i] = r1[j][i] - r2[j][i];
        }
    }
    return res;
  }
}

/*----------------------------------------------------------------------------*/

Matrix
mulMatrix (Matrix m1, Matrix m2)
{
  int h1 = m1->height;
  int w1 = m1->width;
  int h2 = m2->height;
  int w2 = m2->width;

  assert (w1 == h2 || h1 == w2);

  {
    Matrix res = newMatrix (h1, w2);
    Lines r = res->lines;
    Lines r1 = m1->lines;
    Lines r2 = m2->lines;

    int i,j,k;
    for (i = 0; i < h1; ++i)
    {
	      for (j = 0; j < w2; ++j)
	      {
            float sum = 0;
            for (k = 0; k < w1; ++k)
            {
                sum = sum + r1[i][k] * r2[k][j];
            }
            r[i][j] = sum;
        }
    }
    return res;
  }
}

/*----------------------------------------------------------------------------*/

Matrix
transposeMatrix (Matrix m)
{
  int h = m->height;
  int w = m->width;

  Matrix res = newMatrix (w, h);
  Lines r = res->lines;
  Lines r1 = m->lines;

  int j;
  for (j = 0; j < h; ++j)
    {
      int i;
      for (i = 0; i < w; ++i)
	    {
	      r[i][j] = r1[j][i];
      }
    }
  return res;
}

/*----------------------------------------------------------------------------*/

Matrix
choleskyMatrix (Matrix m)
{
     int w = m->width;
     int h = m->height;
     
     int i, j, k;
     float v = 0.0 , Cii;
     
     Matrix res = zeroMatrix (h, w);
     Lines Chol = res->lines;
     Lines C = m->lines;
     
     int DIM = w;
     
     assert (w == h);                
     
     // calculate the cholesky decomposition
     // 1st step :
     if (C[0][0] <= 0.0)
     {
        printf ("Cholesky decomposition failed. Matrix is probably not positive definite.\n");
        return 0;
     }
     Chol[0][0] = sqrt(C[0][0]);
     for (j = 1; j < DIM; j++)
     {
         Chol[0][j] = C[j][0] / Chol[0][0];
     }
     // 2nd to n-th step
     for (i = 1; i < DIM; i++)
     {
         v = C[i][i];
         for (k = i-1; k >= 0; k--)
         {
             v -= Chol [k][i]* Chol[k][i];
         }
         if (v <= 0.0)
         {
            printf ("Cholesky decomposition failed. Matrix is probably not positive definite. sum :%f\n", v);
            return 0 ;
         }
         Chol[i][i] = sqrt(v);
         for (j = i+1; j < DIM; j++)
         {
             v = C[j][i];
             for (k = i-1; k >= 0; k--)
             {
                 v -= Chol[k][i] * Chol[k][j];
             }
             Chol[i][j] = v / Chol[i][i];
         }
     }
     // set the values above the diagonal equal to zero
     for (i = 0; i < DIM; i++)
     {
         for (j = i+1; j < DIM; j++)
         {
             Chol[j][i] = 0.0;
         }
     }
  return res;
}

/*----------------------------------------------------------------------------*/

Matrix
invertCovMatrix (Matrix m)
{
     int w = m->width;
     int h = m->height;
     
     int i, j, k;
     float v = 0.0;
     float Cii;
     
     Matrix res = zeroMatrix (h, w);
     Lines CInv = res->lines;
     Lines C = m->lines;
     
     int DIM = w;
     
     assert (w == h);                
     
     // calculate the cholesky decomposition
     // 1st step :
     if (C[0][0] < 0.0)
     {
        printf ("Inverse cholesky failed. Matrix is probably not positive definite.\n");
        return 0;
     }
     CInv[0][0] = sqrt(C[0][0]);
     for (j = 1; j < DIM; j++)
     {
         CInv[0][j] = C[j][0] / CInv[0][0];
     }
     // 2nd to n-th step
     for (i = 1; i < DIM; i++)
     {
         v = C[i][i];
         for (k = i-1; k >= 0; k--)
         {
             v -= CInv [k][i]* CInv[k][i];
         }
         if (v < 0.0)
         {
            printf ("Inverse cholesky failed. Matrix is probably not positive definite. sum :%f\n", v);
            return 0 ;
         }
         CInv[i][i] = sqrt(v);
         for (j = i+1; j < DIM; j++)
         {
             v = C[j][i];
             for (k = i-1; k >= 0; k--)
             {
                 v -= CInv[k][i] * CInv[k][j];
             }
             CInv[i][j] = v / CInv[i][i];
         }
     }
  //
  // Compute inverse of upper triangular matrix .
  //
  for (j = 0; j < DIM; j++)
  {
      if(CInv[j][j] < 0.0)
      {
          return 0;
      }
      CInv[j][j] = 1.0 / CInv[j][j];
    //
    // Compute elements 0:j -1 of j-th column .
    //
      for (i = 0; i < j; i++)
      {
          v = 0.0;
          for (k = i; k < j; k++)
          {
              v += CInv[i][k] * CInv[k][j];
          }
          CInv[i][j] = v * (- CInv[j][j]);
      }
  }
  //
  // InvC = InvL * InvL ’
  //
  for (i = 0; i < DIM; i++)
  {
      Cii = CInv[i][i];
      if( i < DIM-1)
      {
          v = 0.0;
          for (k = i; k < DIM; k++)
          {
              v += CInv[i][k]*CInv[i][k];
          }
          CInv [i][i] = v;
          for (j = 0; j < i; j++)
          {
              v = 0.0;
              for (k = i+1; k < DIM; k++)
              {
                  v += CInv[j][k] * CInv[i][k];
              }
              CInv[j][i] = CInv[j][i] * Cii + v;
              CInv[i][j] = CInv[j][i];
          }
      }
      else
      {
          for (j = 0; j < DIM; j++)
          {
              CInv[j][i] *= Cii;
              CInv[i][j] = CInv[j][i];
          }
      }
  }
  return res;
}

/*----------------------------------------------------------------------------*/

Matrix
mulScalarMatrix (float a, Matrix m)
{
  int w = m->width;
  int h = m->height;
  Matrix res = newMatrix(h, w);
  Lines r2 = res->lines;
  Lines r1 = m->lines;  
  
  int i,j;
  for (i=0; i<h; i++)
  {
      for (j=0; j<w; j++)
      {
          r2[i][j] = a * r1[i][j];
      }
  }
  return res;
}

/*----------------------------------------------------------------------------*/

Matrix
sumMatrix(Matrix m, int dim)
{
  int h = m->height;
  int w = m->width;
  Lines rm = m->lines;
  int i,j;

  if (dim == 1)
  {
     Matrix res = zeroMatrix(1, w);
     Lines sum = res->lines;
     for (i=0; i<w; i++)
     {
         for (j=0; j<h; j++)
         {
           sum[0][i] = sum[0][i] + rm[j][i];
         }
     }
     return res;
  }
  
  else if (dim == 2)
  {
     Matrix res = zeroMatrix(h, 1);
     Lines sum = res->lines;
     for (i=0; i<h; i++)
     {
         for (j=0; j<w; j++)
         {
           sum[0][i] = sum[0][i] + rm[i][j];
         }
     }  
     return res;
  }
  else
  {
      printf("Error: Dim has to be 1 or 2.\n");
  }
  return 0;   
}

/*----------------------------------------------------------------------------*/

Matrix
appendMatrix (Matrix m1, int x1, int x2, int y1, int y2, Matrix m2, int s1, int s2, int t1, int t2)
{
  /* index check */
  assert ( m1->height >= x2 && m1->width >= y2 );
  assert ( m2->height >= s2 && m2->width >= t2 );
  
  Lines r1 = m1->lines;
  Lines r2 = m2->lines;
  
  int i, j, k, l;
  
  k = s1;
  for (i=x1; i<=x2; i++) 
  {
    l = t1;
    for (j=y1; j<=y2; j++)
    {
      r1[i][j] = r2[k][l];
      l++;
    }    
    k++;
  }
  return m1;  
}

/*----------------------------------------------------------------------------*/

void
appMatrix (Matrix m1, int x1, int x2, int y1, int y2, Matrix m2, int s1, int s2, int t1, int t2)
{
  /* index check */
  assert ( m1->height >= x2 && m1->width >= y2 );
  assert ( m2->height >= s2 && m2->width >= t2 );
  
  Lines r1 = m1->lines;
  Lines r2 = m2->lines;
  
  int i, j, k, l;
  
  k = s1;
  for (i=x1; i<=x2; i++) 
  {
    l = t1;
    for (j=y1; j<=y2; j++)
    {
      r1[i][j] = r2[k][l];
      l++;
    }    
    k++;
  }
}

/*----------------------------------------------------------------------------*/

Element
elem (Matrix m, int i, int j)
{
  /* index check */
  assert (0 <= i && i < m->height);
  assert (0 <= j && j < m->width);

  return m->lines[i][j];
}

void
setElem (Matrix m, int i, int j, Element v)
{
  /* index check */
  assert (0 <= i && i < m->height);
  assert (0 <= j && j < m->width);

  m->lines[i][j] = v;
}

/*----------------------------------------------------------------------------*/

void
printMatrix (Matrix m)
{
  int w = m->width;
  int h = m->height;

  Lines r = m->lines;
  
  int i;
  for (i = 0; i < h; ++i)
  {
      int j;
      for (j = 0; j < w; ++j)
      {
          printf("%0.4f ", r[i][j]);
      }
    printf("\n");
  }
  printf("\n");
}

/*----------------------------------------------------------------------------*/


