#include "blocksqp.hpp"

namespace blockSQP
{

SQPiterate::SQPiterate( Problemspec* prob, SQPoptions* param, bool full )
{
    int maxblocksize = 1;

    // Set nBlocks structure according to if we use block updates or not
    if( param->blockHess == 0 || prob->nBlocks == 1 )
    {
        nBlocks = 1;
        blockIdx = new int[2];
        blockIdx[0] = 0;
        blockIdx[1] = prob->nVar;
        maxblocksize = prob->nVar;
        param->objSecondDerv = 0;
    }
    else if( param->blockHess == 2 && prob->nBlocks > 1 )
    {// hybrid strategy: 1 block for constraints, 1 for objective
        nBlocks = 2;
        blockIdx = new int[3];
        blockIdx[0] = 0;
        blockIdx[1] = prob->blockIdx[prob->nBlocks-1];
        blockIdx[2] = prob->nVar;
    }
    else
    {
        nBlocks = prob->nBlocks;
        blockIdx = new int[nBlocks+1];
        for( int k=0; k<nBlocks+1; k++ )
        {
            blockIdx[k] = prob->blockIdx[k];
            if( k > 0 )
                if( blockIdx[k] - blockIdx[k-1] > maxblocksize )
                    maxblocksize = blockIdx[k] - blockIdx[k-1];
        }
    }

    if( param->hessLimMem && param->hessMemsize == 0 )
        param->hessMemsize = maxblocksize;

    allocMin( prob );

    #ifdef QPSOLVER_DENSE
    constrJac.Dimension( prob->nCon, prob->nVar ).Initialize( 0.0 );
    #else
    jacNz = NULL;
    jacIndCol = NULL;
    jacIndRow = NULL;

    hessNz = NULL;
    hessIndCol = NULL;
    hessIndRow = NULL;
    hessIndLo = NULL;
    #endif

    if( full )
    {
        allocHess();
        allocAlg( prob, param );
    }
}


SQPiterate::SQPiterate( SQPiterate *iter )
{
    int i;

    nBlocks = iter->nBlocks;
    blockIdx = new int[nBlocks+1];
    for( i=0; i<nBlocks+1; i++ )
        blockIdx[i] = iter->blockIdx[i];

    xi = Matrix( iter->xi );
    lambda = Matrix( iter->lambda );
    constr = Matrix( iter->constr );
    gradObj = Matrix( iter->gradObj );
    gradLagrange = Matrix( iter->gradLagrange );

    #ifdef QPSOLVER_DENSE
    constrJac = Matrix( iter->constrJac );
    #else
    int nVar = xi.M();
    int nnz = iter->jacIndCol[nVar];

    jacNz = new double[nnz];
    for( i=0; i<nnz; i++ )
        jacNz[i] = iter->jacNz[i];

    jacIndRow = new int[nnz + (nVar+1) + nVar];
    for( i=0; i<nnz + (nVar+1) + nVar; i++ )
        jacIndRow[i] = iter->jacIndRow[i];
    jacIndCol = jacIndRow + nnz;
    #endif
}


/**
 * Allocate memory for variables
 * required by all optimization
 * algorithms except for the Jacobian
 */
void SQPiterate::allocMin( Problemspec *prob )
{
    // current iterate
    xi.Dimension( prob->nVar ).Initialize( 0.0 );

    // dual variables (for general constraints and variable bounds)
    lambda.Dimension( prob->nVar + prob->nCon ).Initialize( 0.0 );

    // constraint vector with lower and upper bounds
    // (Box constraints are not included in the constraint list)
    constr.Dimension( prob->nCon ).Initialize( 0.0 );

    // gradient of objective
    gradObj.Dimension( prob->nVar ).Initialize( 0.0 );

    // gradient of Lagrangian
    gradLagrange.Dimension( prob->nVar ).Initialize( 0.0 );
}


void SQPiterate::allocHess()
{
    int iBlock, varDim;

    // Create one Matrix for one diagonal block in the Hessian
    hess = new SymMatrix[nBlocks];
    for( iBlock=0; iBlock<nBlocks; iBlock++ )
    {
        varDim = blockIdx[iBlock+1] - blockIdx[iBlock];
        hess[iBlock].Dimension( varDim ).Initialize( 0.0 );
    }
}


/**
 * Convert array *hess to a single symmetric sparse matrix in
 * Harwell-Boeing format (as used by qpOASES)
 */
void SQPiterate::convertHessian( Problemspec *prob )
{
    int iBlock, count, colCountTotal, rowOffset, i, j;
    int nnz, nCols, nRows;

    // 1) count nonzero elements
    nnz = 0;
    for( iBlock=0; iBlock<nBlocks; iBlock++ )
        for( i=0; i<hess[iBlock].N(); i++ )
            for( j=i; j<hess[iBlock].N(); j++ )
                if( hess[iBlock]( i,j ) != 0.0 )
                {
                    nnz++;
                    if( i != j )// off-diagonal elements count twice
                        nnz++;
                }

    if( hessNz != NULL ) delete[] hessNz;
    if( hessIndRow != NULL ) delete[] hessIndRow;

    hessNz = new double[nnz];
    hessIndRow = new int[nnz + (prob->nVar+1) + prob->nVar];
    hessIndCol = hessIndRow + nnz;
    hessIndLo = hessIndCol + (prob->nVar+1);

    // 2) store matrix entries columnwise in hessNz
    count = 0; // runs over all nonzero elements
    colCountTotal = 0; // keep track of position in large matrix
    rowOffset = 0;
    for( iBlock=0; iBlock<nBlocks; iBlock++ )
    {
        nCols = hess[iBlock].N();
        nRows = hess[iBlock].M();

        for( i=0; i<nCols; i++ )
        {
            // column 'colCountTotal' starts at element 'count'
            hessIndCol[colCountTotal] = count;

            for( j=0; j<nRows; j++ )
                if( hess[iBlock]( i,j ) != 0.0 )
                {
                    hessNz[count] = hess[iBlock]( i, j );
                    hessIndRow[count] = j + rowOffset;
                    count++;
                }
            colCountTotal++;
        }

        rowOffset += nRows;
    }
    hessIndCol[colCountTotal] = count;

    // 3) Set reference to lower triangular matrix
    for( j=0; j<prob->nVar; j++ )
    {
        for( i=hessIndCol[j]; i<hessIndCol[j+1] && hessIndRow[i]<j; i++);
        hessIndLo[j] = i;
    }

    if( count != nnz )
         printf( "Error in convertHessian: %i elements processed, should be %i elements!\n", count, nnz );
 }


/**
 * Allocate memory for additional variables
 * needed by the algorithm
 */
void SQPiterate::allocAlg( Problemspec *prob, SQPoptions *param )
{
    int iBlock, varDim, i;
    int hessCount;
    int nVar = prob->nVar;
    int nCon = prob->nCon;

    // current step
    deltaMat.Dimension( nVar, param->hessMemsize, nVar ).Initialize( 0.0 );
    deltaXi.Submatrix( deltaMat, nVar, 1, 0, 0 );
    // trial step (temporary variable, for line search)
    trialXi.Dimension( nVar, 1, nVar ).Initialize( 0.0 );

    // bounds for step (QP subproblem)
    deltaBl.Dimension( nVar+nCon ).Initialize( 0.0 );
    deltaBu.Dimension( nVar+nCon ).Initialize( 0.0 );

    // product of constraint Jacobian with step (deltaXi)
    AdeltaXi.Dimension( nCon ).Initialize( 0.0 );

    // state of each variable (for QP solver)
    istate = new int[nVar+nCon];
    for( i=0; i<nVar+nCon; i++ ) istate[i] = 0;

    // dual variables of QP (simple bounds and general constraints)
    lambdaQP.Dimension( nVar+nCon ).Initialize( 0.0 );

    // line search parameters
    deltaH.Dimension( nBlocks ).Initialize( 0.0 );

    // filter as a set of pairs
    filter = new std::set< std::pair<double,double> >;

    // difference of Lagrangian gradients
    gammaMat.Dimension( nVar, param->hessMemsize, nVar ).Initialize( 0.0 );
    gamma.Submatrix( gammaMat, nVar, 1, 0, 0 );

    // Scalars that are used in various Hessian update procedures
    deltaNorm.Dimension( nBlocks ).Initialize( 1.0 );
    deltaGamma.Dimension( nBlocks ).Initialize( 0.0 );
    deltaBdelta.Dimension( nBlocks ).Initialize( 0.0 );
    noUpdateCounter = new int[nBlocks];
    for( iBlock=0; iBlock<nBlocks; iBlock++ )
        noUpdateCounter[iBlock] = -1;

    // For modified BFGS Updates: for each block save delta^T delta and gamma^T delta
    deltaNormOld.Dimension( nBlocks ).Initialize( 1.0 );
    deltaGammaOld.Dimension( nBlocks ).Initialize( 0.0 );

    updateSequence = new int[param->hessMemsize];
    for( i=0; i<param->hessMemsize; i++ ) updateSequence[i] = 1;
}


void SQPiterate::initIterate( SQPoptions* param )
{
    alpha = 1.0;
    alphaSOC = 0.0;
    reducedStepCount = 0;
    steptype = 0;

    obj = param->inf;
    tol = param->inf;
    cNorm = param->thetaMax;
    gradNorm = param->inf;
    lambdaStepNorm = 0.0;
}

///\todo destructor

} // namespace blockSQP
