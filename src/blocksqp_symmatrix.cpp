#include "blocksqp_matrix.hpp"

namespace blockSQP
{

int SymMatrix::malloc( void )
{
    int len;

    len = m*(m+1)/2.0;

    if ( len == 0 )
       array = NULL;
    else
       if ( ( array = new double[len] ) == NULL )
          Fehler("new erfolglos");

    return 0;
}


 int SymMatrix::free( void )
 {
     if( array != NULL )
         delete[] array;

     return 0;
 }


double &SymMatrix::operator()( int i, int j )
{
    #ifdef MATRIX_DEBUG
    if ( i < 0 || i >= m || j < 0 || j >= n )
    Fehler("Ungueltiger Matrixeintrag");
    #endif

    int pos;

    if( i < j )//reference to upper triangular part
        pos = (int) (j + i*(m - (i+1.0)/2.0));
    else
        pos = (int) (i + j*(m - (j+1.0)/2.0));

    return array[pos];
}


double &SymMatrix::operator()( int i )
{
    Fehler("SymMatrix doesn't support this operation.");
    return array[0];
}


double SymMatrix::a( int i, int j ) const
{
    #ifdef MATRIX_DEBUG
    if ( i < 0 || i >= m || j < 0 || j >= n )
    Fehler("Ungueltiger Matrixeintrag");
    #endif

    int pos;

    if( j > i )//reference to upper triangular part
        pos = (int) (j + i*(m - (i+1.0)/2.0));
    else
        pos = (int) (i + j*(m - (j+1.0)/2.0));

    return array[pos];
}


double SymMatrix::a( int i ) const
{
    Fehler("SymMatrix doesn't support this operation.");
    return 0.0;
}


SymMatrix::SymMatrix( int M )
{
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
}


SymMatrix::SymMatrix( int M, double *ARRAY )
{
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
    array = ARRAY;
}


SymMatrix::SymMatrix( int M, int N, int LDIM )
{
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
}


SymMatrix::SymMatrix( int M, int N, double *ARRAY, int LDIM )
{
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
    array = ARRAY;
}


SymMatrix::SymMatrix( const Matrix &A )
{
    int i, j;

    m = A.M();
    n = A.M();
    ldim = A.M();
    tflag = 0;

    malloc();

    for ( j=0; j<m; j++ )//columns
         for ( i=j; i<m; i++ )//rows
             (*this)(i,j) = A.a(i,j);
}


SymMatrix::SymMatrix( const SymMatrix &A )
{
    int i, j;

    m = A.m;
    n = A.n;
    ldim = A.ldim;
    tflag = 0;

    malloc();

    for ( j=0; j<m; j++ )//columns
         for ( i=j; i<m; i++ )//rows
             (*this)(i,j) = A.a(i,j);
}


SymMatrix::~SymMatrix( void )
{
    Dcount++;

    if( !tflag )
        free();
}


SymMatrix &SymMatrix::operator=( const Matrix &A )
{
    int i, j;

    if( this != &A )
    {
        free();

        m = A.M();
        n = A.N();
        ldim = A.LDIM();

        malloc();

        for ( j=0; j<m; j++ )
            for ( i=j; i<n ; i++ )
                (*this)(i,j) = A.a(i,j);
    }

    return *this;
}


SymMatrix &SymMatrix::Dimension( int M )
{
    free();
    m = M;
    n = M;
    ldim = M;

    malloc();

    return *this;
}


SymMatrix &SymMatrix::Dimension( int M, int N, int LDIM )
{
    free();
    m = M;
    n = M;
    ldim = M;

    malloc();

    return *this;
}


SymMatrix &SymMatrix::set_m( int M )
{
    Fehler("SymMatrix doesn't support this operation.");
    return *this;
}


SymMatrix &SymMatrix::Initialisieren( double (*f)( int, int ) )
{
    int i, j;

    for ( j=0; j<m; j++ )
        for ( i=j; i<n ; i++ )
            (*this)(i,j) = f(i,j);

    return *this;
}


SymMatrix &SymMatrix::Initialisieren( double val )
{
    int i, j;

    for ( j=0; j<m; j++ )
        for ( i=j; i<n ; i++ )
            (*this)(i,j) = val;

    return *this;
}


SymMatrix &SymMatrix::Streiche_Zeilen( int newm )
{
    Fehler("SymMatrix doesn't support this operation.");
    return *this;
}


SymMatrix &SymMatrix::Streiche_Spalten( int newn )
{
    Fehler("SymMatrix doesn't support this operation.");
    return *this;
}


SymMatrix &SymMatrix::Teilmatrix( const Matrix &A, int M, int N, int i0, int j0 )
{
    Fehler("SymMatrix doesn't support this operation.");
    return *this;
}


SymMatrix &SymMatrix::Arraymatrix( int M, double *ARRAY )
{
    if( !tflag )
        free();

    tflag = 1;
    m = M;
    n = M;
    ldim = M;
    array = ARRAY;

    return *this;
}


SymMatrix &SymMatrix::Arraymatrix( int M, int N, double *ARRAY, int LDIM )
{
    if( !tflag )
        free();

    tflag = 1;
    m = M;
    n = M;
    ldim = M;
    array = ARRAY;

    return *this;
}


SymMatrix &SymMatrix::Teilmatrix_del( void )
{
    Fehler("SymMatrix doesn't support this operation.");
    return *this;
}


// Matrix &SymMatrix::convert2Matrix( void )
// {
//     Fehler("...to be implemented.");
// }

} // namespace blockSQP
