#ifndef gaussianEstimator_h
#define gaussianEstimator_h


extern void gaussianEstimator_Pred_decomp
 ( Matrix *xEst, Matrix *CEst, Matrix *U, Matrix *Cw, float *dt, Matrix *m_opt);

extern void gaussianEstimator_Pred
 ( Matrix *xEst, Matrix *CEst, Matrix *U, Matrix *Cw, Matrix (*afun)(Matrix m, float t), float *dt, Matrix *m_opt);
 
extern void gaussianEstimator_Est
 ( Matrix *xEst, Matrix *CEst, Matrix *y, Matrix *Cv, Matrix (*hfun)(Matrix m), Matrix *m_opt); 
 
extern Matrix scaledSamplePoints(Matrix m, Matrix m_opt);

extern Matrix fillMatrix(Matrix x, int d);

#endif
