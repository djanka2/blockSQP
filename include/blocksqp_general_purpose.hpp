#ifndef GENERAL_PURPOSE_HPP
#define GENERAL_PURPOSE_HPP

#include "lapack.h"
#include "blocksqp_defs.hpp"
#include "blocksqp_matrix.hpp"

namespace blockSQP
{

double l1VectorNorm( const Matrix &v );
double l2VectorNorm( const Matrix &v );
double lInfVectorNorm( const Matrix &v );

double l1ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl, const Matrix &weights );
double l1ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl );
double l2ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl );
double lInfConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl );

double adotb( const Matrix &a, const Matrix &b );
void Atimesb( const Matrix &A, const Matrix &b, Matrix &Atimesb );
void sparseAtimesb( double *Anz, int *AIndRow, int *AIndCol, const Matrix &b, Matrix &Atimesb );

int calcEigenvalues( const Matrix &B, Matrix &ev );
double estimateSmallestEigenvalue( const Matrix &B );
int inverse( const Matrix &A, Matrix &Ainv );

} // namespace blockSQP

#endif
