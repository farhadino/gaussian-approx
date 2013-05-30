#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "noise.h"
#include "gaussianEstimator.h"
#include "tracker.h"

/*----------------------------------------------------------------------------*/

Matrix
gaussianApprox_old ( int L )
{
  Matrix res = newMatrix (1, (L-1) );
  Lines m = res->lines;
  
  FILE *file;
	file = fopen( "GaussianApprox.dat", "r");    
             
  int i = 0;
  int start;
  int eof = 143;
  float diracs[eof]; //143

  while(!feof(file)) { 
    // loop through and store the numbers in the array
    fscanf(file, "%f", &diracs[i]);
    i++;
  }
  
  for (i = 0; i < eof; ++i)
  {
    if (diracs[i] == L)
    {
       start = i + 1;
       break;
    }
  }

  int k = 0;
  for (i = start; i < (start + (L-1)); ++i)
  {
    m[0][k] = diracs[i];
    k++;
  }
   
  return res;
}

/*----------------------------------------------------------------------------*/

Matrix
gaussianApprox ( int L )
{
  Matrix res = newMatrix (1, (L-1) );
  Lines m = res->lines;
  
  if (L == 3) {
        m[0][0] = -1.224744930527570;
        m[0][1] =  1.224744930527570;
  }
  else if (L == 5) {
        m[0][0] = -1.479493430695565;
        m[0][1] = -0.557762663290440;
        m[0][2] =  0.557762663290440;
        m[0][3] =  1.479493430695565;
  }
  else {
       m[0][0] = -1.634559525853960;
       m[0][1] = -0.827490279295782;
       m[0][2] = -0.378780931601150;
       m[0][3] =  0.378780931601150;
       m[0][4] =  0.827490279295782;
       m[0][5] =  1.634559525853960;
  }
  
  return res;
}

/*----------------------------------------------------------------------------*/
