#include "blocksqp.hpp"

namespace blockSQP
{

/**
 * Derived class from generic problem specification class that implements all abstract methods
 */
class MyProblem : public Problemspec
{
    public:
        Matrix xi0;

    public:
        MyProblem( int nVar_, int nCon_, int nBlocks_, int *BlockIdx_, const Matrix &bl_, const Matrix &bu_, const Matrix &xi0_ );
        virtual void convertJacobian( const Matrix &constrJac, double *&jacNz, int *&jacIndRow, int *&jacIndCol, bool firstCall = 0 );
        virtual void initialize( Matrix &xi, Matrix &lambda, Matrix &constrJac );
        virtual void initialize( Matrix &xi, Matrix &lambda, double *&jacNz, int *&jacIndRow, int *&jacIndCol );
        virtual void printInfo();

        virtual void evaluate( const Matrix &xi, const Matrix &lambda, double *objval, Matrix &constr,
                               Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                               SymMatrix *&hess, int dmode, int *info );
        virtual void evaluate( const Matrix &xi, const Matrix &lambda, double *objval, Matrix &constr,
                               Matrix &gradObj, Matrix &constrJac, SymMatrix *&hess, int dmode, int *info );
};

void MyProblem::printInfo()
{
    printf( "blockSQP standalone version using\n" );
#ifdef QPSOLVER_QPOPT
    printf( "QPOPT.\n" );
#endif
#ifdef QPSOLVER_QPOASES_DENSE
    printf( "QPOASES (DENSE).\n" );
#endif
#ifdef QPSOLVER_QPOASES_SPARSE
    printf( "QPOASES (w SPARSE matrices).\n" );
#endif
#ifdef QPSOLVER_QPOASES_SCHUR
    printf( "QPOASES (SCHUR COMPLEMENT).\n" );
#endif

    printf( "\nnVar = %i\n", nVar );
    printf( "nCon = %i\n\n", nCon );
}

/**
 * Constraint evaluation methods still operate on dense Jacobian, this is a generic method to convert it
 * to a sparse matrix in Harwell--Boeing (column compressed) format
 */
void MyProblem::convertJacobian( const Matrix &constrJac, double *&jacNz, int *&jacIndRow, int *&jacIndCol, bool firstCall )
{
    int nnz, count, i, j;

    if( firstCall )
    {
        // 1st run: count nonzeros
        nnz = 0;
        for( j=0; j<nVar; j++ )
            for( i=0; i<nCon; i++ )
                if( fabs( constrJac( i, j ) < myInf ) )
                    nnz++;

        if( jacNz != NULL ) delete[] jacNz;
        if( jacIndRow != NULL ) delete[] jacIndRow;

        jacNz = new double[nnz];
        jacIndRow = new int[nnz + (nVar+1) + nVar];
        jacIndCol = jacIndRow + nnz;
    }
    else
    {
        nnz = jacIndCol[nVar];
        /* arrays jacInd* are already allocated! */
    }

    // 2nd run: store matrix entries columnwise in jacNz
    count = 0;
    for( j=0; j<nVar; j++ )
    {
        jacIndCol[j] = count;
        for( i=0; i<nCon; i++ )
            if( fabs( constrJac( i, j ) < myInf ) )
            {
                jacNz[count] = constrJac( i, j );
                jacIndRow[count] = i;
                count++;
            }
    }
    jacIndCol[nVar] = count;
    if( count != nnz )
         printf( "Error in convertJacobian: %i elements processed, should be %i elements!\n", count, nnz );
}

/**
 * Set your initial values for xi (and possibly lambda) here.
 * You can also set the parts of the Jacobian that correspond to purely linear constraints, i.e. that don't change
 */
void MyProblem::initialize( Matrix &xi, Matrix &lambda, Matrix &constrJac )
{
    // set initial values for xi and lambda
    lambda.Initialize( 0.0 );
    for( int i=0; i<nVar; i++ )
        xi( i ) = xi0( i );
}

/**
 * Sparse version of initialize. Sparse Jacobian is allocated here.
 */
void MyProblem::initialize( Matrix &xi, Matrix &lambda, double *&jacNz, int *&jacIndRow, int *&jacIndCol )
{
    Matrix constrDummy, gradObjDummy, constrJac;
    SymMatrix *hessDummy;
    double objvalDummy;
    int info;

    // set initial values for xi and lambda
    lambda.Initialize( 0.0 );
    for( int i=0; i<nVar; i++ )
        xi( i ) = xi0( i );

    // find out Jacobian sparsity pattern
    constrDummy.Dimension( nCon ).Initialize( 0.0 );
    gradObjDummy.Dimension( nVar ).Initialize( 0.0 );
    constrJac.Dimension( nCon, nVar ).Initialize( myInf );
    evaluate( xi, lambda, &objvalDummy, constrDummy, gradObjDummy, constrJac, hessDummy, 1, &info );

    // allocate sparse Jacobian structures
    convertJacobian( constrJac, jacNz, jacIndRow, jacIndCol, 1 );
}


/**
 * Set basic problem data in the constructor:
 * - nVar: number of variables
 * - nCon: number of constraints
 * - bl[nVar+nCon], bu[nVar+nCon]: lower and upper bounds for variables and constraints
 * - nBlocks: number of diagonal blocks in the Hessian
 * - blockIdx[nBlocks+1]: where do blocks start and end?
 *
 * Optional:
 * - conNames[nCon], varNames[nVar]: names for constraints and variables
 */
MyProblem::MyProblem( int nVar_, int nCon_, int nBlocks_, int *blockIdx_, const Matrix &bl_, const Matrix &bu_, const Matrix &xi0_ )
{
    nVar = nVar_;
    nCon = nCon_;

    nBlocks = nBlocks_;
    blockIdx = new int[nBlocks+1];
    if( nBlocks == 1 )
    {
        blockIdx[0] = 0;
        blockIdx[1] = nVar;
    }
    else
    {
        for( int i=0; i<nBlocks+1; i++ )
            blockIdx[i] = blockIdx_[i];
    }

    bl.Dimension( nVar + nCon ).Initialize( -myInf );
    bu.Dimension( nVar + nCon ).Initialize( myInf );

    for( int i=0; i<nVar+nCon; i++ )
    {
        bl( i ) = bl_( i );
        bu( i ) = bu_( i );
    }

    objLo = -myInf;
    objUp = myInf;

    xi0.Dimension( nVar );
    for( int i=0; i<nVar; i++ )
        xi0( i ) = xi0_( i );
}

/**********************************************************************/
/**                     PROBLEM SPECIFIC PART                        **/
/**********************************************************************/

/**
 * Evaluate functions and derivatives (sparse).
 * This is just a wrapper for the dense version.
 */
void MyProblem::evaluate( const Matrix &xi, const Matrix &lambda, double *objval, Matrix &constr,
                          Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                          SymMatrix *&hess, int dmode, int *info )
{
    Matrix constrJac;

    constrJac.Dimension( nCon, nVar ).Initialize( myInf );
    evaluate( xi, lambda, objval, constr, gradObj, constrJac, hess, dmode, info );

    // Convert to sparse format
    if( dmode != 0 )
        convertJacobian( constrJac, jacNz, jacIndRow, jacIndCol, 0 );
}

/**
 * Evaluate functions and derivatives.
 * dmode = 0: objective and constraint evaluation
 * dmode = 1: objective gradient and constraint Jacobian
 */
void MyProblem::evaluate( const Matrix &xi, const Matrix &lambda, double *objval, Matrix &constr,
                          Matrix &gradObj, Matrix &constrJac, SymMatrix *&hess,
                          int dmode, int *info )
{
    *info = 0;
    int i;

    if( dmode >= 0 )
    {
        *objval = xi(0)*xi(0);
        for( i=1; i<nVar; i++ )
            *objval -= i * 0.5*xi(i)*xi(i);

        for( i=1; i<nVar; i++ )
            constr( i-1 ) = xi(0) - sqrt(i*(nVar-1))*xi(i);
    }

    if( dmode > 0 )
    {
        gradObj( 0 ) = 2.0 * xi(0);
        for( i=1; i<nVar; i++ )
            gradObj( i ) = - i * xi(i);

        for( i=1; i<nVar; i++ )
        {
            constrJac( i-1, 0 ) = 1.0;
            constrJac( i-1, i ) = - sqrt(i*(nVar-1));
        }
    }
}

} // namespace blockSQP

int main( int argc, const char* argv[] )
{
    using namespace blockSQP;
    int ret = 0;
    MyProblem *prob;
    SQPmethod *meth;
    SQPoptions *opts;
    SQPstats *stats;
    char outpath[255];
    strcpy( outpath, "./" );


    /*--------------------*/
    /* Setup problem data */
    /*--------------------*/

    int nVar = 2;

    int nCon = nVar-1;

    // Initial values
    Matrix x0;
    x0.Dimension( nVar );
    x0( 0 ) = 10.0;
    for( int i=1; i<nVar; i++ )
        x0( i ) = x0( 0 ) * sqrt(i*(nVar-1));

    // Variable bounds
    Matrix bl, bu;
    bl.Dimension( nVar+nCon ).Initialize( -myInf );
    bu.Dimension( nVar+nCon ).Initialize( myInf );

    // Constraint bounds
    for( int i=nVar; i<nVar+nCon; i++ )
    {
        bl( i ) = bu( i ) = 0.0;
    }

    // Variable partition for block Hessian
    int nBlocks = nVar;
    int blockIdx[nBlocks+1];
    for( int i=0; i<nBlocks+1; i++ )
        blockIdx[i] = i;

    // Create problem evaluation object
    prob = new MyProblem( nVar, nCon, nBlocks, blockIdx, bl, bu, x0 );

    /*------------------------*/
    /* Options for SQP solver */
    /*------------------------*/

    opts = new SQPoptions();
    opts->opttol = 1.0e-10;
    opts->nlinfeastol = 1.0e-10;

    // 0: no globalization, 1: filter line search
    opts->globalization = 1;
    // 0: (scaled) identity, 1: SR1, 2: BFGS
    opts->hessUpdate = 4;
    // 0: initial Hessian is diagonal matrix, 1: scale initial Hessian according to Nocedal p.143,
    // 2: scale initial Hessian with Oren-Luenberger factor 3: scale initial Hessian with geometric mean of 1 and 2
    // 4: scale Hessian in every step with centered Oren-Luenberger sizing according to Tapia paper
    opts->hessScaling = 0;
    // scaling strategy for fallback BFGS update if SR1 and globalization is used
    opts->fallbackScaling = 0;
    // Size of limited memory
    opts->hessMemsize = 200;
    // If too many updates are skipped, reset Hessian
    opts->maxConsecSkippedUpdates = 200;
    // 0: full space Hessian approximation (ignore block structure), 1: blockwise updates
    opts->blockHess = 1;


    /*-------------------------------------------------*/
    /* Create blockSQP method object and run algorithm */
    /*-------------------------------------------------*/

    stats = new SQPstats( outpath );
    meth = new SQPmethod( prob, opts, stats );

    ret = meth->init();
    ret = meth->run( 100 );
    meth->finish();
    if( ret == 1 )
        printf("\033[0;36m***Maximum number of iterations reached.***\n\033[0m");

    printf("\nPrimal solution:\n");
    meth->vars->xi.Print();
    printf("\nDual solution:\n");
    meth->vars->lambda.Print();
    printf("\nHessian approximation at the solution:\n");
    for( int i=0; i<meth->vars->nBlocks; i++ )
        meth->vars->hess[i].Print();
    printf("\nFallback Hessian at the solution:\n");
    for( int i=0; i<meth->vars->nBlocks; i++ )
        meth->vars->hess2[i].Print();

    // Clean up
    delete prob;
    delete stats;
    delete opts;
    delete meth;
}

