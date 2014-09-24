#include "blocksqp_problemspec.hpp"

namespace blockSQP
{

void Problemspec::evaluate( const Matrix &xi, double *objval, Matrix &constr, int *info )
{
    Matrix lambdaDummy, gradObjDummy;
    SymMatrix *hessDummy;

    int dmode = 0;

    #ifdef QPSOLVER_SPARSE
    double *jacNzDummy;
    int *jacIndRowDummy, *jacIndColDummy;
    evaluate( xi, lambdaDummy, objval, constr, gradObjDummy, jacNzDummy, jacIndRowDummy, jacIndColDummy, hessDummy, dmode, info );
    #else
    Matrix constrJacDummy;
    evaluate( xi, lambdaDummy, objval, constr, gradObjDummy, constrJacDummy, hessDummy, dmode, info );
    #endif
}

} // namespace blockSQP
