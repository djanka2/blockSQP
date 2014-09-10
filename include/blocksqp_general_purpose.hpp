#ifndef GENERAL_PURPOSE_HPP
#define GENERAL_PURPOSE_HPP

#include "lapack.h"
#include "blocksqp_defs.hpp"
#include "blocksqp_matrix.hpp"

namespace blockSQP
{

double l1VectorNorm( Matrix v );
double l2VectorNorm( Matrix v );
double lInfVectorNorm( Matrix v );

double l1ConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl, Matrix weights );
double l1ConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl );
double l2ConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl );
double lInfConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl );

double adotb( Matrix a, Matrix b );
void Atimesb( Matrix A, Matrix b, Matrix &Atimesb );
void sparseAtimesb( double *Anz, int *AIndRow, int *AIndCol, Matrix b, Matrix &Atimesb );

int calcEigenvalues( SymMatrix B, Matrix &ev );
double estimateSmallestEigenvalue( Matrix B );
int inverse( Matrix A, Matrix &Ainv );

} // namespace blockSQP

#endif
