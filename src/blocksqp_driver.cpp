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
        Matrix a;

    public:
        MyProblem( int nVar_, int nCon_, int nBlocks_, int *BlockIdx_, Matrix bl_, Matrix bu_, Matrix xi0_ );
        virtual void convertJacobian( Matrix constrJac, double *&jacNz, int *&jacIndRow, int *&jacIndCol, bool firstCall = 0 );
        virtual void initialize( Matrix &xi, Matrix &lambda, Matrix &constrJac );
        virtual void initialize( Matrix &xi, Matrix &lambda, double *&jacNz, int *&jacIndRow, int *&jacIndCol );
        virtual void printInfo();

        virtual void evaluate( Matrix xi, Matrix lambda, double *objval, Matrix &constr,
                               Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                               SymMatrix *&hess, int dmode, int *info );
        virtual void evaluate( Matrix xi, Matrix lambda, double *objval, Matrix &constr,
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
void MyProblem::convertJacobian( Matrix constrJac, double *&jacNz, int *&jacIndRow, int *&jacIndCol, bool firstCall )
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
    lambda.Initialisieren( 0.0 );
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
    lambda.Initialisieren( 0.0 );
    for( int i=0; i<nVar; i++ )
        xi( i ) = xi0( i );

    // find out Jacobian sparsity pattern
    constrDummy.Dimension( nCon ).Initialisieren( 0.0 );
    gradObjDummy.Dimension( nVar ).Initialisieren( 0.0 );
    constrJac.Dimension( nCon, nVar ).Initialisieren( myInf );
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
 * - objLo, objUp: lower and upper bound for objective (should not be active at the solution)
 *
 * Optional:
 * - conNames[nCon], varNames[nVar]: names for constraints and variables
 */
MyProblem::MyProblem( int nVar_, int nCon_, int nBlocks_, int *blockIdx_, Matrix bl_, Matrix bu_, Matrix xi0_ )
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

    bl.Dimension( nVar + nCon ).Initialisieren( -myInf );
    bu.Dimension( nVar + nCon ).Initialisieren( myInf );

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
void MyProblem::evaluate( Matrix xi, Matrix lambda, double *objval, Matrix &constr,
                           Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                           SymMatrix *&hess, int dmode, int *info )
{
    Matrix constrJac;

    constrJac.Dimension( nCon, nVar ).Initialisieren( myInf );
    evaluate( xi, lambda, objval, constr, gradObj, constrJac, hess, dmode, info );

    // Convert to sparse format
    if( dmode != 0 )
        convertJacobian( constrJac, jacNz, jacIndRow, jacIndCol, 0 );
}

#define MINEX
//#define OEDTYPE
//#define OEDTYPE_LIFTED
#define NMESS 10

/**
 * Evaluate functions and derivatives.
 * dmode = 0: objective and constraint evaluation
 * dmode = 1: objective gradient and constraint Jacobian
 * dmode = 2: objective gradient, constraint Jacobian and Hessian of the Lagrangian (not ready to use yet)
 */
void MyProblem::evaluate( Matrix xi, Matrix lambda, double *objval, Matrix &constr,
                           Matrix &gradObj, Matrix &constrJac, SymMatrix *&hess,
                           int dmode, int *info )
{
    *info = 0;

#ifdef OEDTYPE
    if( dmode >= 0 )
    {
        *objval = 0.0;
        for( int i=0; i<NMESS; i++ )
            *objval += ( a(i)*a(i)*xi(i)*xi(i) );
        *objval = 1.0 / *objval;
    }

    if( dmode > 0 )
    {
        for( int i=0; i<NMESS; i++ )
            gradObj( i ) = -(*objval)*(*objval) * 2.0*a(i)*a(i)*xi(i);
    }
#endif

#ifdef OEDTYPE_LIFTED
    double y = xi( nVar-1 );

    if( dmode >= 0 )
    {
        *objval = 1.0 / y;

        constr( 0 ) = 0.0;
        for( int i=0; i<NMESS; i++ )
            constr( 0 ) += ( a(i)*a(i)*xi(i)*xi(i) );
        constr( 0 ) -= y;
    }

    if( dmode > 0 )
    {
        for( int i=0; i<NMESS; i++ )
        {
            gradObj( i ) = 0.0;
            constrJac( 0, i ) = 2.0*a(i)*a(i)*xi(i);
        }
        gradObj( NMESS ) = -1.0 / ( y*y );
        constrJac( 0, NMESS ) = -1.0;
    }
#endif

#ifdef MINEX
    if( dmode >= 0 )
    {
        *objval = xi(0)*xi(0) - 0.5*xi(1)*xi(1);
        constr( 0 ) = xi(0) - xi(1);
    }

    if( dmode > 0 )
    {
        gradObj( 0 ) = 2.0 * xi(0);
        gradObj( 1 ) = - xi(1);

        constrJac( 0, 0 ) = 1.0;
        constrJac( 0, 1 ) = -1.0;
    }
#endif

    if( dmode > 1 )
    {
        printf( "Second derivatives not supported!\n" );
        *info = -1;
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

    // Options for SQP solver
    opts = new SQPoptions();
    opts->opttol = 1.0e-12;
    opts->nlinfeastol = 1.0e-12;
    opts->conSecondDerv = 0;
    opts->objSecondDerv = 0;
    opts->hessUpdate = 2;
    opts->hessScaling = 0;
    opts->fallbackScaling = 0;
    opts->hessLimMem = 1;
    opts->hessMemsize = 200;
    opts->maxConsecSkippedUpdates = 200;
    opts->blockHess = 1;
    opts->restoreFeas = 0;
    opts->globalization = 0;
    opts->iniHessDiag = 1.0;

    ////////////////////////////////////////////////////////////////////

    // Measurement weights
    Matrix a;
    a.Dimension( NMESS );
    for( int i=1; i<NMESS+1; i++ )
        a( i-1 ) = 0.1*i;

#ifdef OEDTYPE
    // Setup problem data
    int nVar = NMESS;
    int nCon = 0;

    // Variable bounds
    Matrix bl, bu;
    bl.Dimension( nVar+nCon ).Initialisieren( 0.5 );
    bu.Dimension( nVar+nCon ).Initialisieren( 1 );

    // Initial values
    Matrix x0;
    x0.Dimension( nVar );
    for( int i=0; i<NMESS; i++ )
        x0( i ) = bl( i );

    // Variable partition for block Hessian
    int nBlocks = 1;
    int blockIdx[nBlocks+1];
    blockIdx[0] = 0;
    blockIdx[1] = nVar;
#endif

#ifdef OEDTYPE_LIFTED
        // Setup problem data
    int nVar = NMESS+1;
    int nCon = 1;

    // Variable bounds
    Matrix bl, bu;
    bl.Dimension( nVar+nCon ).Initialisieren( -myInf );
    bu.Dimension( nVar+nCon ).Initialisieren( myInf );
    for( int i=0; i<NMESS; i++ )
    {
        bl( i ) = 0.5;
        bu( i ) = 1;
    }

    // Initial values
    Matrix x0;
    x0.Dimension( nVar ).Initialisieren( 0.0 );
    for( int i=0; i<NMESS; i++ )
    {
        x0( i ) = bl( i );
        x0( nVar-1 ) += a(i)*a(i)*x0(i)*x0(i);
    }

    // Constraint bounds
    bl( nVar ) = bu( nVar ) = 0.0;

    // Variable partition for block Hessian
    int nBlocks = NMESS+1;
    int blockIdx[nBlocks+1];
    for( int i=0; i<NMESS+1; i++ )
        blockIdx[i] = i;
    blockIdx[NMESS+1] = nVar;
#endif

#ifdef MINEX
    // Setup problem data
    int nVar = 2;
    int nCon = 1;

    // Initial values
    Matrix x0;
    x0.Dimension( nVar );
    x0( 0 ) = 10.0;
    x0( 1 ) = 10.0;

    // Variable bounds
    Matrix bl, bu;
    bl.Dimension( nVar+nCon ).Initialisieren( -myInf );
    bu.Dimension( nVar+nCon ).Initialisieren( myInf );

    // Constraint bounds
    bl( nVar ) = bu( nVar ) = 0.0;

    // Variable partition for block Hessian
    int nBlocks = 2;
    int blockIdx[nBlocks+1];
    blockIdx[0] = 0;
    blockIdx[1] = 1;
    blockIdx[2] = nVar;
#endif
    ////////////////////////////////////////////////////////////////////

    // Create problem (evaluation) object
    prob = new MyProblem( nVar, nCon, nBlocks, blockIdx, bl, bu, x0 );

    prob->a = Matrix( a );

    // Create blockSQP method including memory allocation for iterate
    stats = new SQPstats( outpath );
    meth = new SQPmethod( prob, opts, stats );

    ret = meth->init();
    ret = meth->run( 100 );
    meth->finish();
    if( ret == 1 )
        printf("\033[0;36m***Maximum number of iterations reached.***\n\033[0m");

    printf("\nPrimal solution:\n");
    meth->vars->xi.Ausgabe();
    printf("\nDual solution:\n");
    meth->vars->lambda.Ausgabe();
    printf("\nHessian approximation at the solution:\n");
    for( int i=0; i<meth->vars->nBlocks; i++ )
        meth->vars->hess[i].Ausgabe();

    // Clean up
    delete prob;
    delete meth;
    delete stats;
    delete opts;
}

