#ifndef GENERAL_PURPOSE_HPP
#define GENERAL_PURPOSE_HPP

#include "stdlib.h"
#include "math.h"
#include "string.h"

#include "lapack.h"
#include "matrix.hpp"

int oedPhi( int dimH, Matrix vecH, int *row, int *col, double *phi, Matrix &dphi, SymMatrix &ddphi, int dmode );
int calcEigenvalues( SymMatrix B, Matrix &ev );
double estimateSmallestEigenvalue( Matrix B );
double fRand( double fMin, double fMax );
double l1VectorNorm( Matrix v );
double l2VectorNorm( Matrix v );
double adotb( Matrix a, Matrix b );
void Atimesb( Matrix A, Matrix b, Matrix &Atimesb );
void sparseAtimesb( double *Anz, int *AIndRow, int *AIndCol, Matrix b, Matrix &Atimesb );
double lInfVectorNorm( Matrix v );
double l1ConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl, Matrix weights );
double l1ConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl );
double l2ConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl );
double lInfConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl );
int symmetricInverse2( Matrix A, Matrix &Ainv );
int symmetricInverse3( SymMatrix A, SymMatrix &Ainv );
void symmetricMultiplication( Matrix A, Matrix B, Matrix &C, char *side, int mode );
int inverse( Matrix A, Matrix &Ainv );
void constructMyInv( int dimH, Matrix measContrib, Matrix Cj, Matrix &myInv, int mode );
void constructMyInv2( int dimH, Matrix CjInv, Matrix measContrib, Matrix &myInv2 );

#endif
