//
// THIS IS A DUPLICATE OF THE VPLAN FILE. IT IS INTENDED FOR THE STANDALONE VERSION OF BLOCKSQP.
//

#include "general_purpose.hpp"


/**
 * Compute eigenvalues of a symmetric matrix by DSPEV
 */
int calcEigenvalues( SymMatrix B, Matrix &ev )
{
    int i, j, n;
    SymMatrix temp;
    double *work, *dummy;
    int info, iDummy = 1;

    n = B.M();
    ev.Dimension( n ).Initialisieren( 0.0 );
    work = new double[3*n];

    // copy Matrix, will be overwritten
    temp = SymMatrix( B );

    // DSPEV computes all the eigenvalues and, optionally, eigenvectors of a
    // real symmetric matrix A in packed storage.
    dspev( "N", "L", &n, temp.ARRAY(), ev.ARRAY(), dummy, &iDummy,
            work, &info, strlen("N"), strlen("L") );

    delete[] work;

    return info;
}

/**
 * Estimate the smalles eigenvalue of a sqare matrix
 * with the help og Gershgorin's circle theorem
 */
double estimateSmallestEigenvalue( Matrix B )
{
    int i, j;
    double radius;
    int dim = B.M();
    double lambdaMin = 0.0;

    // For each row, sum up off-diagonal elements
    for( i=0; i<dim; i++ )
    {
        radius = 0.0;
        for( j=0; j<dim; j++ )
            if( j != i )
                radius += fabs( B( i,j ) );

        if( B( i,i ) - radius < lambdaMin )
            lambdaMin = B( i,i ) - radius;
    }

    return lambdaMin;
}


/**
 * Return a pseudo random double between fMin and fMax
 */
double fRand( double fMin, double fMax )
{
    double f = ( double ) rand() / RAND_MAX;
    return fMin + f * ( fMax - fMin );
}

double l1VectorNorm( Matrix v )
{
    double norm = 0.0;

    if( v.N() != 1 )
    {
        printf("v is not a vector!\n");
    }
    else
    {
        for( int k=0; k<v.M(); k++ )
            norm += fabs(v( k ));
    }

    return norm;
}


double l2VectorNorm( Matrix v )
{
    double norm = 0.0;

    if( v.N() != 1 )
    {
        printf("v is not a vector!\n");
    }
    else
    {
        for( int k=0; k<v.M(); k++ )
            norm += v( k )* v( k );
    }

    return sqrt(norm);
}

double adotb( Matrix a, Matrix b )
{
    double norm = 0.0;

    if( a.N() != 1 || b.N() != 1 )
    {
        printf("a or b is not a vector!\n");
    }
    else if( a.M() != b.M() )
    {
        printf("a and b must have the same dimension!\n");
    }
    else
    {
        for( int k=0; k<a.M(); k++ )
            norm += a(k) * b(k);
    }

    return norm;
}

/**
 * Compute the matrix vector product for a column-compressed sparse matrix A
 */
void sparseAtimesb( double *Anz, int *AIndRow, int *AIndCol, Matrix b, Matrix &Atimesb )
{
    int nCol = b.M();
    int nRow = Atimesb.M();
    int i, k;
    int nnz = AIndCol[nCol];

    for( i=0; i<nRow; i++ )
        Atimesb( i ) = 0.0;

    for( i=0; i<nCol; i++ )
    {
        // k runs over all elements in one column
        for( k=AIndCol[i]; k<AIndCol[i+1]; k++ )
            Atimesb( AIndRow[k] ) += Anz[k] * b( i );
    }

}

/**
 * Compute the matrix vector product
 */
void Atimesb( Matrix A, Matrix b, Matrix &Atimesb )
{
    Atimesb.Initialisieren( 0.0 );
    for( int i=0; i<A.M(); i++ )
        for( int k=0; k<A.N(); k++ )
            Atimesb( i ) += A( i, k ) * b( k );
}

double lInfVectorNorm( Matrix v )
{
    double norm = 0.0;

    if( v.N() != 1 )
    {
        printf("v is not a vector!\n");
    }
    else
    {
        for( int k=0; k<v.M(); k++ )
            if( fabs(v( k )) > norm )
                norm = fabs(v( k ));
    }

    return norm;
}


/**
 * Calculate weighted l1 norm of constraint violations
 */
double l1ConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl, Matrix weights )
{
    double norm = 0.0;
    int i;
    int nVar = xi.M();

    if( weights.M() < constr.M() + nVar )
    {
        printf("Weight vector too short!\n");
        return 0.0;
    }

    // Weighted violation of simple bounds
    for( i=0; i<nVar; i++ )
    {
        if( xi( i ) > bu( i ) )
            norm += fabs(weights( i )) * (xi( i ) - bu( i ));
        else if( xi( i ) < bl( i ) )
            norm += fabs(weights( i )) * (bl( i ) - xi( i ));
    }

    // Calculate weighted sum of constraint violations
    for( i=0; i<constr.M(); i++ )
    {
        if( constr( i ) > bu( nVar+i ) )
            norm += fabs(weights( nVar+i )) * (constr( i ) - bu( nVar+i ));
        else if( constr( i ) < bl( nVar+i ) )
            norm += fabs(weights( nVar+i )) * (bl( nVar+i ) - constr( i ));
    }

    return norm;
}


/**
 * Calculate l1 norm of constraint violations
 */
double l1ConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl )
{
    double norm = 0.0;
    int i;
    int nVar = xi.M();

    // Violation of simple bounds
    for( i=0; i<nVar; i++ )
    {
        if( xi( i ) > bu( i ) )
            norm += xi( i ) - bu( i );
        else if( xi( i ) < bl( i ) )
            norm += bl( i ) - xi( i );
    }

    // Calculate sum of constraint violations
    for( i=0; i<constr.M(); i++ )
    {
        if( constr( i ) > bu( nVar+i ) )
            norm += constr( i ) - bu( nVar+i );
        else if( constr( i ) < bl( nVar+i ) )
            norm += bl( nVar+i ) - constr( i );
    }

    return norm;
}


/**
 * Calculate l2 norm of constraint violations
 */
double l2ConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl )
{
    double norm = 0.0;
    int i;
    int nVar = xi.M();

    // Violation of simple bounds
    for( i=0; i<nVar; i++ )
        if( xi( i ) > bu( i ) )
            norm += xi( i ) - bu( i );
        if( xi( i ) < bl( i ) )
            norm += bl( i ) - xi( i );

    // Calculate sum of constraint violations
    for( i=0; i<constr.M(); i++ )
        if( constr( i ) > bu( nVar+i ) )
            norm += pow(constr( i ) - bu( nVar+i ), 2);
        else if( constr( i ) < bl( nVar+i ) )
            norm += pow(bl( nVar+i ) - constr( i ), 2);

    return sqrt(norm);
}


/**
 * Calculate l_Infinity norm of constraint violations
 */
double lInfConstraintNorm( Matrix xi, Matrix constr, Matrix bu, Matrix bl )
{
    double norm = 0.0;
    int i;
    int nVar = xi.M();
    int nCon = constr.M();

    // Violation of simple bounds
    for( i=0; i<xi.M(); i++ )
    {
        if( xi( i ) - bu( i ) > norm )
            norm = xi( i ) - bu( i );
        else if( bl( i ) - xi( i ) > norm )
            norm = bl( i ) - xi( i );
    }

    // Find out the largest constraint violation
    for( i=0; i<constr.M(); i++ )
    {
        if( constr( i ) - bu( nVar+i ) > norm )
            norm = constr( i ) - bu( nVar+i );
        if( bl( nVar+i ) - constr( i ) > norm )
            norm = bl( nVar+i ) - constr( i );
    }

    return norm;
}


/**
 * Compute the inverse of a matrix
 * using LU decomposition
 */
int inverse( Matrix A, Matrix &Ainv )
{
    int i, j;
    int n, ldim, lwork, info = 0;
    int *ipiv;
    double *work;

    for( i=0; i<A.N(); i++ )
        for( j=0; j<A.M(); j++ )
            Ainv( j,i ) = A( j,i );

    n = Ainv.N();
    ldim = Ainv.LDIM();
    ipiv = new int[n];
    lwork = n*n;
    work = new double[lwork];

    // Compute LU factorization
    dgetrf( &n, &n, Ainv.ARRAY(), &ldim, ipiv, &info );
    if ( info != 0 )
        printf( "WARNING: DGETRF returned info=%i\n", info );
    // Compute inverse
    dgetri( &n, Ainv.ARRAY(), &ldim, ipiv, work, &lwork, &info );
    if ( info != 0 )
        printf( "WARNING: DGETRI returned info=%i\n", info );

    return info;
}


/**
 * Compute the inverse of a symmetric
 * indefinite matrix
 */
int symmetricInverse2( Matrix A, Matrix &Ainv )
{
    int i, j;
    int n, ldim, lwork, *ipiv, info = 0;
    double *work;

    for( i=0; i<A.N(); i++ )
        for( j=i; j<A.M(); j++ )
            Ainv( j,i ) = A( j,i );

    n = Ainv.N();
    ldim = Ainv.LDIM();
    lwork = n*n;

    work = new double[lwork];
    ipiv = new int[n];

    // DSYTRF computes the factorization of a real symmetric matrix A using
    // the Bunch-Kaufman diagonal pivoting method.
    dsytrf( "L", &ldim, Ainv.ARRAY(), &ldim, ipiv,
            work, &lwork, &info, strlen("L") );
    if ( info != 0 )
    {
        printf( "WARNING: DSYTRF returned info=%i\n", info );
//         A.Ausgabe();
//         Ainv.Ausgabe();
    }

    // DSYTRI computes the inverse of a real symmetric indefinite matrix
    // A using the factorization A = U*D*U**T or A = L*D*L**T computed by
    // DSYTRF.
    dsytri( "L", &ldim, Ainv.ARRAY(), &ldim, ipiv,
            work, &info, strlen("L") );
    if ( info != 0 )
    {
        printf( "WARNING: DSYTRI returned info=%i\n", info );
//         A.Ausgabe();
//         Ainv.Ausgabe();
    }

    for( i=0; i<A.N(); i++ )
        for( j=i; j<A.M(); j++ )
            Ainv( i,j ) = Ainv( j,i );

    delete[] ipiv;
    delete[] work;

    return info;
}


/**
 * Compute the inverse of a symmetric
 * indefinite matrix in packed storage
 */
int symmetricInverse3( SymMatrix A, SymMatrix &Ainv )
{
    int i, j;
    int n, *ipiv, info = 0;
    double *work;

    for( i=0; i<A.N(); i++ )
        for( j=i; j<A.M(); j++ )
            Ainv( j,i ) = A( j,i );

    n = A.N();

    work = new double[n];
    ipiv = new int[n];

    // DSPTRF computes the factorization of a real symmetric matrix A stored
    // in packed format using the Bunch-Kaufman diagonal pivoting method:
    //
    //    A = U*D*U**T  or  A = L*D*L**T
    //
    // where U (or L) is a product of permutation and unit upper (lower)
    // triangular matrices, and D is symmetric and block diagonal with
    // 1-by-1 and 2-by-2 diagonal blocks.
    dsptrf( "L", &n, Ainv.ARRAY(), ipiv, &info, strlen("L") );
    if ( info != 0 )
    {
        printf( "WARNING: DSPTRF returned info=%i\n", info );
         //A.Ausgabe();
         //Ainv.Ausgabe();
    }

    // DSPTRI computes the inverse of a real symmetric indefinite matrix
    // A in packed storage using the factorization A = U*D*U**T or
    // A = L*D*L**T computed by DSPTRF.
    dsptri( "L", &n, Ainv.ARRAY(), ipiv, work, &info, strlen("L") );
    if ( info != 0 )
    {
        printf( "WARNING: DSPTRI returned info=%i\n", info );
//         A.Ausgabe();
//         Ainv.Ausgabe();
    }

    delete[] ipiv;
    delete[] work;

    return info;
}


/**
 * Input: A as vectorized lower triangular matrix (mode=1)
 *        or as Matrix (mode=2),
 *        B as Matrix, C as Matrix, side=L or R
 * Output: C = A*B + C as Matrix  (side=L)
 *      or C = B*A + C as Matrix  (side=R)
 *
 * Only square matrices allowed!
 */
void symmetricMultiplication( Matrix A, Matrix B, Matrix &C, char *side, int mode )
{
    int dim, ldim, i, j, count;
    Matrix temp1, temp2;
    double alpha, beta;
    alpha = beta = 1.0;

    dim = C.M();
    ldim = C.LDIM();

    if( mode == 0 )
    {
        temp1.Dimension( dim, dim, dim ).Initialisieren( 0.0 );
        count = 0;
        for( i=0; i<dim; i++ )// columns
            for( j=i; j<dim; j++ )// rows
                temp1( i, j ) = temp1( j, i ) = A( count++ );

        temp2.Dimension( dim, dim, dim ).Initialisieren( 0.0 );
        count = 0;
        for( i=0; i<dim; i++ )// columns
            for( j=i; j<dim; j++ )// rows
                temp2( i, j ) = temp2( j, i ) = B( count++ );
    }
    else if( mode == 1 )
    {
        temp1.Dimension( dim, dim, dim ).Initialisieren( 0.0 );
        count = 0;
        for( i=0; i<dim; i++ )// columns
            for( j=i; j<dim; j++ )// rows
                temp1( i, j ) = temp1( j, i ) = A( count++ );

        temp2.Arraymatrix( dim, dim, B.ARRAY() );

    }
    else if( mode == 2 )
    {
        temp1.Arraymatrix( dim, dim, A.ARRAY() );
        temp2.Arraymatrix( dim, dim, B.ARRAY() );
    }

    // C = alpha*temp*B + beta*C;
    dsymm( side, "L", &dim, &dim,
           &alpha, temp1.ARRAY(), &dim,
                   temp2.ARRAY(), &ldim,
           &beta, C.ARRAY(), &ldim,
           1, 1 );
}


/**
 * Input: measContrib and Cj as vectorized
 * lower triangular matrices.
 *
 * Output: myInv = (I+measContrib*Cj)^-1
 * as matrix (in general not symmetric!)
 */
void constructMyInv( int dimH, Matrix measContrib, Matrix Cj, Matrix &myInv, int mode )
{
    int i, j, count;
    Matrix temp;
    double alpha, beta;

    alpha = beta = 1.0;

    temp.Dimension( dimH, dimH, dimH ).Initialisieren( 0.0 );
    count = 0;
    for( i=0; i<dimH; i++ )// columns
        for( j=i; j<dimH; j++ )// rows
            temp( i, j ) = temp( j, i ) = Cj( count++ );

    // Identity matrix
    myInv.Dimension( dimH, dimH , dimH ).Initialisieren( ::delta );

    if( mode == 0 )
        symmetricMultiplication( measContrib, temp, myInv, "L", 1 );
    else if( mode == 1 )
        symmetricMultiplication( measContrib, temp, myInv, "R", 1 );

    inverse( myInv, myInv );

}


/**
 * Input: measContrib as vectorized
 * lower triangular matrix, CjInv as square matrix
 *
 * Output: myInv2 = (Cj^-1+measContrib)^-1
 * as (symmetric) square matrix
 */
void constructMyInv2( int dimH, Matrix CjInv, Matrix measContrib, Matrix &myInv2 )
{
    int i, j, count;
    Matrix temp;
    double alpha, beta;

    alpha = beta = 1.0;

    temp.Dimension( dimH, dimH, dimH ).Initialisieren( 0.0 );
    count = 0;
    for( i=0; i<dimH; i++ )// columns
        for( j=i; j<dimH; j++ )// rows
            temp( i, j ) = temp( j, i ) = CjInv( i,j ) + measContrib( count++ );

    myInv2.Dimension( dimH, dimH , dimH ).Initialisieren( 0.0 );

    symmetricInverse2( temp, myInv2);
}
