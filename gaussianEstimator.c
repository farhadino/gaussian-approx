#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "noise.h"
#include "gaussianEstimator.h"
#include "eig.h"
#include "tracker.h"

/*----------------------------------------------------------------------------*/

void
gaussianEstimator_Pred_decomp ( Matrix *xEst, Matrix *CEst, Matrix *U, Matrix *Cw, float *dt, Matrix *m_opt)
{
float r;

Matrix sizeMopt;
Matrix xn = zeroMatrix(3,1);
Matrix Cn = zeroMatrix(3,3);
Matrix xl = zeroMatrix(9,1);
Matrix Cl = zeroMatrix(9,9);
Matrix Cnl = zeroMatrix(3,9);
Matrix Cnl_T;
Matrix Cn_i;
Matrix CLN;
Matrix sizeCn;
Matrix Vec;
Matrix Val;
Matrix m1;
Matrix m;
Matrix x;
Matrix A;
Matrix Hi = zeroMatrix(12,9);
Matrix Cy = zeroMatrix(12, 12);
Matrix muy = zeroMatrix(12, 1);
Matrix zeros33 = zeroMatrix(3,3);
Matrix eye33 = unitMatrix(3,3);
Matrix Mat;
Matrix H;
Matrix gi = zeroMatrix(12,1);
Matrix Rot_vec = zeroMatrix(3,1);
Matrix mui;
Matrix muiy;
Matrix Ciy;
Matrix tmp;
Matrix tmp1;
Matrix tmp2;
Matrix tmp3;
Matrix tmp4;
Matrix tmp5;
Matrix tmp6;
Matrix tmp7;
Matrix tmp8;
Matrix tmpHi; 
                              
sizeMopt = sizeOfMatrix(*m_opt);                      //printf("%f\n",*dt);
float D = elem(sizeMopt,0,1)+1;                       //printf("%f\n", D);
freeMatrix(sizeMopt);
float w_opt = 1/D;                                    //printf("%f\n", w_opt);
float Nx = 3;
float d = Nx*(D-1) + 1;                               //printf("%f\n", d);
float w = 1/d;                                        //printf("%f\n", w);



//xn = xEst(4:6); % Rotation vector
appMatrix(xn, 0, 2, 0, 0, *xEst, 3, 5, 0, 0);         //printMatrix(xn);system("PAUSE");

//Cn = CEst(4:6,4:6);
appMatrix(Cn, 0, 2, 0, 2, *CEst, 3, 5, 3, 5);         //printMatrix(Cn);system("PAUSE");

//xl = [xEst(1:3) ; xEst(7:12)]; % Translation, angular velocity, linear velocity
appMatrix(xl, 0, 2, 0, 0, *xEst, 0, 2, 0, 0);
appMatrix(xl, 3, 8, 0, 0, *xEst, 6, 11, 0, 0);         //printMatrix(xl);system("PAUSE");

//Cl = [CEst(1:3,1:3) CEst(1:3,7:12);
//      CEst(7:12,1:3) CEst(7:12,7:12)] ;
appMatrix(Cl, 0, 2, 0, 2, *CEst, 0, 2, 0, 2);
appMatrix(Cl, 0, 2, 3, 8, *CEst, 0, 2, 6, 11);
appMatrix(Cl, 3, 8, 0, 2, *CEst, 6, 11, 0, 2);
appMatrix(Cl, 3, 8, 3, 8, *CEst, 6, 11, 6, 11);        //printMatrix(Cl);system("PAUSE");

//Cnl = [CEst(4:6,1:3) CEst(4:6,7:12)];
appMatrix(Cnl, 0, 2, 0, 2, *CEst, 3, 5, 0, 2);
appMatrix(Cnl, 0, 2, 3, 8, *CEst, 3, 5, 6, 11);      //printMatrix(Cnl);system("PAUSE");

//CLN = Cl - Cnl'*inv(Cn)*Cnl;
Cnl_T = transposeMatrix(Cnl);
                                    // printMatrix(Cn);system("PAUSE");
Cn_i = invertCovMatrix(Cn);         //printMatrix(Cn_i);system("PAUSE");

tmp = mulMatrix( Cnl_T, Cn_i);
tmp7 = mulMatrix(tmp, Cnl);
CLN = subMatrix ( Cl,  tmp7);                //printMatrix(CLN);system("PAUSE");
freeMatrix(tmp);
freeMatrix(tmp7);

// Eigenvectors, Eigenvalues
sizeCn = sizeOfMatrix(Cn);
int dimC = elem ( sizeCn, 0, 0 );
freeMatrix(sizeCn);
Vec = zeroMatrix(dimC, dimC);    
Val = zeroMatrix(dimC, dimC);

eig ( &Cn, &Vec, &Val );    //printMatrix(Cn);printMatrix(Vec);printMatrix(Val);system("PAUSE");

// m1 = vec*sqrtf(val)
int i;
for ( i = 0; i < dimC; ++i )
    setElem(Val, i, i, sqrtf(fabs(elem(Val, i,i))));
m1 = mulMatrix(Vec, Val);           //printMatrix(m1);system("PAUSE");

//  rotate & scale samples: m = m1*S
m = scaledSamplePoints(m1, *m_opt); //printMatrix(m);system("PAUSE");
// x = x*ones(1,d)
x = fillMatrix(xn, d);
// shift samples: m = m + x
tmp = addMatrix(m, x);
appMatrix(m, 0, m->height-1, 0, m->width-1, tmp, 0, tmp->height-1, 0, tmp->width-1 );     //printMatrix(m);system("PAUSE");
freeMatrix(tmp);
//A = [[eye(3,3),t*eye(3,3)];[zeros(3,3),eye(3,3)]];
A = unitMatrix(6,6);
setElem(A, 0, 3, *dt);
setElem(A, 1, 4, *dt);
setElem(A, 2, 5, *dt);                              //printMatrix(A);system("PAUSE");

for (i=0; i<d; i++)
{
    //gi = [zeros(3,1); m(:,i); zeros(6,1)];
    setElem(gi, 3, 0, elem(m, 0, i));
    setElem(gi, 4, 0, elem(m, 1, i));    
    setElem(gi, 5, 0, elem(m, 2, i));               //printMatrix(gi);system("PAUSE");
    //Rot_vec = m(:,i);
    setElem(Rot_vec, 0, 0, elem(m, 0, i));
    setElem(Rot_vec, 1, 0, elem(m, 1, i));
    setElem(Rot_vec, 2, 0, elem(m, 2, i));          //printMatrix(Rot_vec);system("PAUSE");

    //r = norm(Rot_vec);
    r = sqrtf( powf((elem(Rot_vec,0,0)),2) + powf((elem(Rot_vec,1,0)),2) + powf((elem(Rot_vec,2,0)),2) );  //printf("%f\n",r);

    H = zeroMatrix(3,3);

    if (fmod(r, 2*pi) == 0)
       {
         Mat = unitMatrix(3,3);
       }
        
        
    else
       { 
        // build skew symmetric Matrix
        setElem(H, 0, 1, -elem(Rot_vec,2,0));
        setElem(H, 0, 2,  elem(Rot_vec,1,0));
        setElem(H, 1, 0,  elem(Rot_vec,2,0));
        setElem(H, 1, 2, -elem(Rot_vec,0,0));
        setElem(H, 2, 0, -elem(Rot_vec,1,0));
        setElem(H, 2, 1,  elem(Rot_vec,0,0));      //printMatrix(H);system("PAUSE");
        // Bortz equation 
        // Mat = eye(3,3) + 0.5*H + (1- r*sin(r)/( 2*(1-cos(r))))/r^2*H*H;
        // already declared Mat = unitMatrix(3,3);
        tmp1 = mulScalarMatrix(0.5, H);
        tmp4 = addMatrix( eye33 , tmp1 );
        tmp2 = mulMatrix(H, H);
        tmp3 = mulScalarMatrix( (1-(r*sin(r)/(2*(1-cos(r)))))/powf(r,2), tmp2);
        Mat = addMatrix( tmp4, tmp3);
                                               //printMatrix(Mat);system("PAUSE");
        freeMatrix(tmp1);
        freeMatrix(tmp2);
        freeMatrix(tmp3);
        freeMatrix(tmp4);

       }
    
    //Hi = [[A(1:3,1:3) zeros(3,3) A(1:3,4:6)];
    //     [zeros(3,3), t*Mat, zeros(3,3)];
    //     [zeros(3,3), eye(3,3), zeros(3,3)];
    //     [A(4:6,1:3),zeros(3,3), A(4:6,4:6)]];
    
    appMatrix( Hi, 0, 2, 0, 2, A,       0, 2, 0, 2 );
    appMatrix( Hi, 0, 2, 3, 5, zeros33, 0, 2, 0, 2 );
    appMatrix( Hi, 0, 2, 6, 8, A,       0, 2, 3, 5 );

    appMatrix( Hi, 3, 5, 0, 2, zeros33, 0, 2, 0, 2 );
    tmpHi = mulScalarMatrix(*dt, Mat);
    appMatrix( Hi, 3, 5, 3, 5, tmpHi,     0, 2, 0, 2 );
    freeMatrix(tmpHi);
    appMatrix( Hi, 3, 5, 6, 8, zeros33, 0, 2, 0, 2 );
    
    appMatrix( Hi, 6, 8, 0, 2, zeros33, 0, 2, 0, 2 );
    appMatrix( Hi, 6, 8, 3, 5, eye33,   0, 2, 0, 2 );
    appMatrix( Hi, 6, 8, 6, 8, zeros33, 0, 2, 0, 2 );
    
    appMatrix( Hi, 9, 11, 0, 2, A,       3, 5, 0, 2 );
    appMatrix( Hi, 9, 11, 3, 5, zeros33, 0, 2, 0, 2 );
    appMatrix( Hi, 9, 11, 6, 8, A,       3, 5, 3, 5 );     //printMatrix(Hi);system("PAUSE");
    
    // mui = xl + Cnl'*inv(Cn)*(m(:,i)-xn);   //m(:,i) -> Rot_vec
    tmp = mulMatrix(Cnl_T, Cn_i );
    tmp1 = subMatrix(Rot_vec, xn);
    tmp2 = mulMatrix(tmp, tmp1);
    mui = addMatrix(xl, tmp2);
    freeMatrix(tmp);
    freeMatrix(tmp1);
    freeMatrix(tmp2);                                   //printMatrix(mui);system("PAUSE");
        
    // muiy = gi + Hi * mui;
    tmp = mulMatrix(Hi, mui);
    muiy = addMatrix( gi, tmp);     //printMatrix(muiy);system("PAUSE");
    freeMatrix(tmp);
    
    // Ciy = Hi *CLN *Hi';
    tmp1 = mulMatrix(Hi, CLN);
    tmp2 = transposeMatrix(Hi);
    Ciy = mulMatrix( tmp1, tmp2);  //printMatrix(Ciy);system("PAUSE");
    freeMatrix(tmp1);
    freeMatrix(tmp2);
     
    // Cy = Cy + (w*Ciy + w_opt*muiy*muiy');
    tmp3 = mulScalarMatrix(w, Ciy);
    tmp1 = transposeMatrix(muiy);
    tmp2 = mulMatrix(muiy, tmp1);
    tmp4 = mulScalarMatrix( w_opt, tmp2 );
    tmp5 = addMatrix( tmp3, tmp4 );
    tmp6 = addMatrix( Cy, tmp5);
    appMatrix(Cy,0,Cy->height-1,0,Cy->width-1,tmp6, 0,tmp6->height-1,0,tmp6->width-1);  //printMatrix(Cy);system("PAUSE");
    freeMatrix(tmp1);
    freeMatrix(tmp2);
    freeMatrix(tmp3);
    freeMatrix(tmp4);
    freeMatrix(tmp5);
    freeMatrix(tmp6);

    // muy = muy + w*muiy;
    tmp = mulScalarMatrix( w, muiy );
    tmp2 = addMatrix( muy, tmp ); 
    appMatrix(muy,0,muy->height-1,0,muy->width-1, tmp2, 0, tmp2->height-1, 0, tmp2->width-1);  //printMatrix(muy);system("PAUSE");
    freeMatrix(tmp);
    freeMatrix(tmp2);

    freeMatrix(H);
    freeMatrix(Mat);
    freeMatrix(mui);//
    freeMatrix(muiy);//
    freeMatrix(Ciy);
       
} 


appMatrix(*xEst, 0, 11, 0, 0, muy, 0, 11, 0, 0 );                       //printMatrix(muy);system("PAUSE");

//CEst = Cy - muy*muy' * w_opt/w + Cw;
tmp1 = transposeMatrix(muy);
tmp2 = mulMatrix(muy, tmp1);
tmp5 = mulScalarMatrix( w_opt/w, tmp2 );
tmp6 = subMatrix(Cy, tmp5);
tmp8 = addMatrix( tmp6, *Cw);           //printMatrix(*CEst);system("PAUSE");
appMatrix(*CEst,0,11,0,11, tmp8, 0,11,0,11 );                          //printMatrix(tmp8);system("PAUSE");
freeMatrix(tmp1);
freeMatrix(tmp2);
freeMatrix(tmp5);
freeMatrix(tmp6);
freeMatrix(tmp8);

freeMatrix(muy);//
freeMatrix(zeros33);//
freeMatrix(Vec);
freeMatrix(Val);
freeMatrix(Cy);
freeMatrix(xn);
freeMatrix(Cn);
freeMatrix(xl);
freeMatrix(Cl);//
freeMatrix(Cnl);
freeMatrix(Cnl_T);
freeMatrix(Cn_i);
freeMatrix(CLN);//
freeMatrix(m1);
freeMatrix(m);//
freeMatrix(x);
freeMatrix(A);
freeMatrix(eye33);
freeMatrix(Hi);
freeMatrix(gi);
freeMatrix(Rot_vec);



} /* End gaussianPred_decomp */

/*----------------------------------------------------------------------------*/

void
gaussianEstimator_Est (Matrix *xEst, Matrix *CEst, Matrix *y, Matrix *Cv, Matrix (*hfun)(Matrix m), Matrix *m_opt)
{
                      //printMatrix(*xEst);
                      //printMatrix(*CEst);system("PAUSE");
Matrix tmp = sizeOfMatrix(*m_opt);                      
float D = elem(tmp,0,1)+1;       //printf("%f\n", D);
freeMatrix(tmp);
float w_opt = 1/D;                                //printf("%f\n", w_opt);
tmp = sizeOfMatrix(*xEst);
float Nx = elem(tmp,0,0);        // printf("%f\n", Nx);
freeMatrix(tmp);
float d = Nx*(D-1) + 1;                           //printf("%f\n", d);
float w = 1/d;                                   // printf("%f\n", w);system("PAUSE");

// Eigenvectors, Eigenvalues
tmp = sizeOfMatrix(*CEst);
int dimC = elem ( tmp, 0, 0 );
freeMatrix(tmp);
Matrix Vec = zeroMatrix(dimC, dimC);    
Matrix Val = zeroMatrix(dimC, dimC);
eig ( CEst, &Vec, &Val );                   //printMatrix(Vec);printMatrix(Val);system("PAUSE");

// m1 = vec*sqrtf(val)
int i;
for ( i = 0; i < dimC; ++i )
    setElem(Val, i, i, sqrtf(fabs(elem(Val, i,i))));
Matrix m1 = mulMatrix(Vec, Val);                   //printMatrix(m1); system("PAUSE");
freeMatrix(Vec);
freeMatrix(Val);

//*  rotate & scale samples: m = m1*S 
Matrix m = scaledSamplePoints(m1, *m_opt);        // printMatrix(m); system("PAUSE");
Matrix mxDiracs = mulScalarMatrix(1, m);


//* x = x*ones(1,d)
Matrix x = fillMatrix(*xEst, d);

// shift samples: m = m + x
tmp = addMatrix(m, x);
appMatrix(m, 0, m->height-1, 0, m->width-1, tmp, 0, tmp->height-1, 0, tmp->width-1 ) ;                              //printMatrix(m);
freeMatrix(tmp);


//% Predicted Measurements
//* hfun 
// yPredDiracs = feval(hfun, m, [], [], t);
// yPred = w*sum(yPredDiracs, 2);

Matrix yPredDiracs = (*hfun) (m);                  //printMatrix(yPredDiracs );  

Matrix yPredDiracsSum = sumMatrix(yPredDiracs, 2); 
Matrix yPred = mulScalarMatrix(w, yPredDiracsSum); 
// myDiracs = yPredDiracs-repmat(yPred, 1, d);
tmp = fillMatrix(yPred, d);
Matrix myDiracs = subMatrix(yPredDiracs, tmp); 
freeMatrix(tmp);

//* CPred = w_opt*mxDiracs*mxDiracs';     
// Matrix CPred = mulScalarMatrix( w_opt, mulMatrix(mxDiracs, transposeMatrix(mxDiracs)) );
// Matrix CPred = *CEst;


// Cxy = w_opt*mxDiracs*myDiracs';
Matrix tmp1 = transposeMatrix(myDiracs);
Matrix tmp2 = mulMatrix(mxDiracs, tmp1);
Matrix Cxy = mulScalarMatrix( w_opt, tmp2);
freeMatrix(tmp1);
freeMatrix(tmp2);


// Cy  = w_opt*myDiracs*myDiracs'+Cv;
tmp1 = transposeMatrix(myDiracs);
tmp2 = mulMatrix(myDiracs, tmp1);
Matrix tmp3 = mulScalarMatrix( w_opt, tmp2);
Matrix Cy = addMatrix( tmp3 , *Cv );
freeMatrix(tmp1);
freeMatrix(tmp2);
freeMatrix(tmp3);

// K = Cxy / Cy;
tmp = invertCovMatrix(Cy);
Matrix K = mulMatrix( Cxy, tmp);
freeMatrix(tmp);

// I = y - yPred;
Matrix I = subMatrix( *y, yPred );

// xEst = xPred + K*I;
tmp = mulMatrix( K, I  );
Matrix tmp23 = addMatrix( *xEst, tmp);
appMatrix(*xEst,0,5,0,0, tmp23,0,5,0,0);
freeMatrix(tmp);


// CEst = CPred - K*Cy*K';
tmp1 = mulMatrix(K, Cy);
tmp2 = transposeMatrix(K);
tmp3 = mulMatrix( tmp1, tmp2);
Matrix tmp24 = subMatrix(*CEst, tmp3);
appMatrix(*CEst,0,5,0,5, tmp24,0,5,0,5);
freeMatrix(tmp1);
freeMatrix(tmp2);
freeMatrix(tmp3);
freeMatrix(tmp24);

freeMatrix(m1);
freeMatrix(m);
freeMatrix(mxDiracs);
freeMatrix(x);
freeMatrix(yPredDiracs);
freeMatrix(yPredDiracsSum);//
freeMatrix(yPred);//
freeMatrix(myDiracs);
freeMatrix(Cxy);
freeMatrix(Cy);
freeMatrix(K);//
freeMatrix(I);//
freeMatrix(tmp23);


}

/*----------------------------------------------------------------------------*/

Matrix 
scaledSamplePoints(Matrix m, Matrix m_opt)
{
  int w = m_opt->width;
  int h = m->height;
  Matrix res = zeroMatrix(h, h*w + 1);
  
  Lines r1 = m->lines;
  Lines r2 = m_opt->lines;
  Lines r3 = res->lines;
    
  int i, j, k, l;
   
  for (i=0; i<h; i++)
  {
      l = 0;
      j = 1;
      while (j<h*w)
      {
          for (k=0; k<w; k++)
          {
              r3[i][j] = r1[i][l] * r2[0][k];
              j++;
          }
          l++;
      }
  }
  return res;
}

/*----------------------------------------------------------------------------*/

Matrix fillMatrix(Matrix x, int d)
{
  int h = x->height;
  Matrix res = newMatrix(h, d);
  Lines r = res->lines;
  Lines x2 = x->lines;
  
  int i,j;
  for (i=0; i<d; i++)
  {
      for (j=0; j<h; j++)
      {
        r[j][i] = x2[j][0];
      }
  }
   
  return res;
}

/*----------------------------------------------------------------------------*/

void
gaussianEstimator_Pred ( Matrix *xEst, Matrix *CEst, Matrix *U, Matrix *Cw, Matrix (*afun)(Matrix m, float t), float *dt, Matrix *m_opt)
{
float D = elem(sizeOfMatrix(*m_opt),0,1)+1;       //printf("%f\n", D);
float w_opt = 1/D;                                //printf("%f\n", w_opt);
float Nx = elem(sizeOfMatrix(*xEst),0,0);         //printf("%f\n", Nx);
float d = Nx*(D-1) + 1;                           //printf("%f\n", d);
float w = 1/d;                                    //printf("%f\n", w);

/* Eigenvectors, Eigenvalues */
int dimC = elem ( sizeOfMatrix(*CEst), 0, 0 );
Matrix Vec = zeroMatrix(dimC, dimC);    
Matrix Val = zeroMatrix(dimC, dimC);
eig ( CEst, &Vec, &Val );

/* m1 = vec*sqrtf(val) */
int i;
for ( i = 0; i < dimC; ++i )
    setElem(Val, i, i, sqrtf(fabs(elem(Val, i,i))));
Matrix m1 = mulMatrix(Vec, Val);

freeMatrix(Vec);
freeMatrix(Val);

/*  rotate & scale samples: m = m1*S */
Matrix m = scaledSamplePoints(m1, *m_opt);

/* x = x*ones(1,d) */
Matrix x = fillMatrix(*xEst, d);

/* shift samples: m = m + x */
m = addMatrix(m, x);              //printMatrix(m);

/* afun */
/* Prediction: mean
   xPredDiracs = feval(afun,m, [], [], t);
   xPred = w*sum(xPredDiracs, 2);*/
   
Matrix xPredDiracs = (*afun) (m, *dt);                   //printMatrix(xPredDiracs);
Matrix xPredDiracsSum = sumMatrix(xPredDiracs, 2);       //printMatrix(xPredDiracsSum);
Matrix xPred = mulScalarMatrix(w, xPredDiracsSum);       //printMatrix(xPred);

//mxDiracs = xPredDiracs-repmat(xPred, 1, d);
//CPred = w_opt*mxDiracs*mxDiracs';
Matrix mxDiracs = subMatrix(xPredDiracs, fillMatrix(xPred, d));   //printMatrix(mxDiracs);
Matrix CPred = mulScalarMatrix(w_opt, mulMatrix(mxDiracs, transposeMatrix(mxDiracs)));  //printMatrix(CPred);


//RETURN
*xEst = xPred;        //printMatrix(*xEst);
*CEst = CPred;        //printMatrix(*CEst);

freeMatrix(m);
freeMatrix(xPredDiracs);
freeMatrix(xPredDiracsSum);
}

/*----------------------------------------------------------------------------*/






