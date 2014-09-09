/***************************************************************************
   VPLAN Source Code
   Version: $Id$
   Author: Stefan Koerkel and JRG ExpDesign
 ***************************************************************************/

/**
 * \file matrix.cc
 */

/* Matrix-Datentyp */
/* von Stefan Koerkel 1993-98 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "matrix.hpp"


// #define MATRIX_DEBUG
#define STANDALONE

/* ----------------------------------------------------------------------- */

#ifdef STANDALONE

void Fehler( const char *F )
{    printf("Fehler: %s\n", F );
     exit( 1 );
}

#endif

/* ----------------------------------------------------------------------- */

int Ccount = 0;
int Dcount = 0;
int Ecount = 0;

/* ----------------------------------------------------------------------- */

int Matrix::malloc( void )
{   int len;

    if ( tflag )
       Fehler("malloc auf Teilmatrix");

    if ( ldim < m )
       ldim = m;

    len = ldim*n;

    if ( len == 0 )
       array = NULL;
    else
       if ( ( array = new double[len] ) == NULL )
          Fehler("new erfolglos");

    return 0;
}


int Matrix::free( void )
{   if ( tflag )
       Fehler("free auf Teilmatrix");

    if ( array != NULL )
       delete[] array;

    return 0;
}


double &Matrix::operator()( int i, int j )
{
     #ifdef MATRIX_DEBUG
     if ( i < 0 || i >= m || j < 0 || j >= n )
        Fehler("Ungueltiger Matrixeintrag");
     #endif

     return array[i+j*ldim];
}


double &Matrix::operator()( int i )
{
     #ifdef MATRIX_DEBUG
     if ( i < 0 || i >= m )
        Fehler("Ungueltiger Matrixeintrag");
     #endif

     return array[i];
}


double Matrix::a( int i, int j ) const
{
     #ifdef MATRIX_DEBUG
     if ( i < 0 || i >= m || j < 0 || j >= n )
        Fehler("Ungueltiger Matrixeintrag");
     #endif

     return array[i+j*ldim];
}


double Matrix::a( int i ) const
{
     #ifdef MATRIX_DEBUG
     if ( i < 0 || i >= m )
        Fehler("Ungueltiger Matrixeintrag");
     #endif

     return array[i];
}

/* ----------------------------------------------------------------------- */

Matrix::Matrix( int M, int N, int LDIM )
{    Ccount++;

     m = M;
     n = N;
     ldim = LDIM;
     tflag = 0;

     malloc();
}


Matrix::Matrix( int M, int N, double *ARRAY, int LDIM )
{    Ccount++;

     m = M;
     n = N;
     array = ARRAY;
     ldim = LDIM;
     tflag = 0;

     if ( ldim < m )
        ldim = m;
}


Matrix::Matrix( const Matrix &A )
{    int i, j;

     Ccount++;

     m = A.m;
     n = A.n;
     ldim = A.ldim;
     tflag = 0;

     malloc();

     for ( i = 0; i < m; i++ )
         for ( j = 0; j < n ; j++ )
             (*this)(i,j) = A.a(i,j);
}


Matrix::~Matrix( void )
{    Dcount++;

     if ( !tflag )
        free();
}

/* ----------------------------------------------------------------------- */

int Matrix::M( void ) const
{   return m;
}


int Matrix::N( void ) const
{   return n;
}


int Matrix::LDIM( void ) const
{   return ldim;
}


double *Matrix::ARRAY( void )
{   return array;
}


int Matrix::TFLAG( void ) const
{   return tflag;
}

/* ----------------------------------------------------------------------- */

Matrix &Matrix::operator=( const Matrix &A )
{    int i, j;

     Ecount++;

     if ( this != &A )
        if ( !tflag )
        {
           /* fuer normale Matrizen */

           free();

           m = A.m;
           n = A.n;
           ldim = A.ldim;

           malloc();

           for ( i = 0; i < m; i++ )
               for ( j = 0; j < n ; j++ )
                   (*this)(i,j) = A.a(i,j);
        }
        else
        {
           /* fuer Teilmatrizen */

           if ( m != A.m || n != A.n )
              Fehler("= Operation nicht erlaubt");

           for ( i = 0; i < m; i++ )
               for ( j = 0; j < n ; j++ )
                   (*this)(i,j) = A.a(i,j);
        }

     return *this;
}

/* ----------------------------------------------------------------------- */

Matrix &Matrix::Dimension( int M, int N, int LDIM )
{    if ( M != m || N != n || ( LDIM != ldim && LDIM != -1 ) )
     {  if ( tflag )
           Fehler("Eine Teilmatrix kann nicht redimensioniert werden");
        else
        {  free();

           m = M;
           n = N;
           ldim = LDIM;

           malloc();
        }
     }

     return *this;
}


Matrix &Matrix::set_m( int M )
{    if ( M > ldim )
        Fehler("Neue Dimension ist groesser als Leading Dimension");

     m = M;

     return *this;
}


Matrix &Matrix::Initialisieren( double (*f)( int, int ) )
{    int i, j;

     for ( i = 0; i < m; i++ )
         for ( j = 0; j < n; j++ )
             (*this)(i,j) = f(i,j);

     return *this;
}


Matrix &Matrix::Initialisieren( double val )
{    int i, j;

     for ( i = 0; i < m; i++ )
         for ( j = 0; j < n; j++ )
             (*this)(i,j) = val;

     return *this;
}


Matrix &Matrix::Streiche_Zeilen( int newm )
{    if ( newm > m )
        Fehler("Es koennen keine Zeilen gestrichen werden");

     m = newm;

     return *this;
}


Matrix &Matrix::Streiche_Spalten( int newn )
{    if ( newn > n )
        Fehler("Es koennen keine Spalten gestrichen werden");

     n = newn;

     return *this;
}

/* ----------------------------------------------------------------------- */

Matrix &Matrix::Teilmatrix( const Matrix &A, int M, int N, int i0, int j0 )
{    if ( i0 + M > A.m || j0 + N > A.n )
        Fehler("Teilmatrix kann nicht gebildet werden");

     if ( !tflag )
        free();

     tflag = 1;

     m = M;
     n = N;
     array = &A.array[i0+j0*A.ldim];
     ldim = A.ldim;

     return *this;
}


Matrix &Matrix::Arraymatrix( int M, int N, double *ARRAY, int LDIM )
{    if ( !tflag )
        free();

     tflag = 1;

     m = M;
     n = N;
     array = ARRAY;
     ldim = LDIM;

     if ( ldim < m )
        ldim = m;

     return *this;
}


Matrix &Matrix::Teilmatrix_del( void )
{    if ( !tflag )
        Fehler("Matrix ist keine Teilmatrix");

     tflag = 0;

     m = 1;
     n = 1;
     ldim = -1;

     malloc();

     return *this;
}

/* ----------------------------------------------------------------------- */

Matrix &Matrix::operator+=( const Matrix &A )
{    int i, j;

     if ( m != A.m || n != A.n )
        Fehler("+= Operation nicht erlaubt");

     for ( i = 0; i < m; i++ )
         for ( j = 0; j < n ; j++ )
             (*this)(i,j) += A.a(i,j);

     return *this;
}


Matrix &Matrix::operator-=( const Matrix &A )
{    int i, j;

     if ( m != A.m || n != A.n )
        Fehler("-= Operation nicht erlaubt");

     for ( i = 0; i < m; i++ )
         for ( j = 0; j < n ; j++ )
             (*this)(i,j) -= A.a(i,j);

     return *this;
}


Matrix &Matrix::operator*=( double s )
{    int i, j;

     for ( i = 0; i < m; i++ )
         for ( j = 0; j < n ; j++ )
             (*this)(i,j) *= s;

     return *this;
}

/* ----------------------------------------------------------------------- */

int Matrix::convert2sparse( int *nrows, int *ncols, int *nnonzeros,
                            int **matbeg, double **matval, int **matind,
                            double eps ) const
{   int i, j, cnt;

    *nrows = m;
    *ncols = n;

    /* 1. Durchlauf: Zaehle Nonzeros und allokiere Arrays */

    *nnonzeros = 0;

    for ( j = 0; j < n; j++ )
    	for ( i = 0; i < m; i++ )
    	    if ( fabs(a(i,j)) > eps )
    	       (*nnonzeros)++;


    if ( ( *matbeg = new int[n+1] ) == NULL )
       Fehler("new erfolglos");
    if ( ( *matval = new double[*nnonzeros] ) == NULL )
       Fehler("new erfolglos");
    if ( ( *matind = new int[*nnonzeros] ) == NULL )
       Fehler("new erfolglos");

    /* 2. Durchlauf: Speichere Eintraege */

    cnt = 0;
    for ( j = 0; j < n; j++ )
    {   (*matbeg)[j] = cnt + 1;               /* FORTRAN-Indizierung */
    	for ( i = 0; i < m; i++ )
    	    if ( fabs(a(i,j)) > eps )
    	    {  (*matval)[cnt] = a(i,j) ;
    	       (*matind)[cnt] = i + 1;        /* FORTRAN-Indizierung */
    	       cnt++;
    	    }
    }
    (*matbeg)[n] = cnt + 1;                   /* FORTRAN-Indizierung */

    return 0;
}

/* ----------------------------------------------------------------------- */

Matrix &Matrix::Eingabe( void )
{    int i, j;

     for ( i = 0; i < m; i++ )
     {   for ( j = 0; j < n; j++ )
         {   printf("Eintrag[%d][%d] = ", i, j );
             scanf("%lf", &(*this)(i,j) );
         }
     }

     return *this;
}


Matrix &Matrix::Einlesen( FILE *f )
{    int theM, theN;
     long pos;

     throw "Matrix::Einlesen is not implemented at the moment";

     // Read a matrix from a file

     // Read dimensions

     pos = ftell( f );

     theM = 0;
     theN = 0;

     Dimension( theM, theN );

     // Read values

     fseek( f, pos, SEEK_SET );

     return *this;
}


const Matrix &Matrix::Ausgabe( FILE *f, int DIGITS, int flag ) const
{    int i, j, l;
     double x;
     //const double eps = 1e-6;

     // Flag == 1: Matlab output
     // else: plain output

     if ( flag == 1 )
        fprintf( f, "[" );

     for ( i = 0; i < m; i++ )
     {   for ( j = 0; j < n; j++ )
         {   x = a(i,j);
             //x = (fabs(x) > eps) ? x : 0.0;

	     if ( flag == 1 )
	     {  fprintf( f, j == 0 ? " " : ", " );
        	fprintf( f, "%.*le", DIGITS, x );
             }
	     else
	     {  fprintf( f, j == 0 ? "" : "  " );
        	fprintf( f, "% .*le", DIGITS, x );
             }
         }
	 if ( flag == 1 )
         {  if ( i < m-1 )
	       fprintf( f, ";\n" );
	 }
	 else
         {  if ( i < m-1 )
	       fprintf( f, "\n" );
	 }
     }

     if ( flag == 1 )
     {  fprintf( f, " ];\n" );
     }
     else
     {  fprintf( f, "\n" );
     }

     return *this;
}

/* ----------------------------------------------------------------------- */

double delta( int i, int j )
{    return (i == j) ? 1.0 : 0.0;
}


double nullentry( int, int )
{    return 0.0;
}


double einsentry( int, int )
{    return 1.0;
}


double randomentry( int, int )
{    const double base = 32767.0;

     return 2.0 * ( rand() / base ) - 1.0;
}

/* ----------------------------------------------------------------------- */

int operator==( const Matrix &A, const Matrix &B )
{    int i, j;
     const double eps = 1e-8;

     if ( A.M() != B.M() || A.N() != B.N() )
        return 0;

     for ( i = 0; i < A.M(); i++ )
         for ( j = 0; j < A.N(); j++ )
             if ( fabs( A.a(i,j) - B.a(i,j) ) > eps )
                return 0;

     return 1;
}


int operator!=( const Matrix &A, const Matrix &B )
{    if ( A == B )
        return 0;
     else
        return 1;
}

/* ----------------------------------------------------------------------- */

Matrix Transponierte( const Matrix &A )
{    int i, j;
     double *array;

     if ( ( array = new double[A.N()*A.M()] ) == NULL )
        Fehler("new erfolglos");

     for ( i = 0; i < A.N(); i++ )
         for ( j = 0; j < A.M(); j++ )
             array[i+j*A.N()] = A.a(j,i);

     return Matrix( A.N(), A.M(), array, A.N() );
}


Matrix &Transponierte( const Matrix &A, Matrix &T )
{    int i, j;

     T.Dimension( A.N(), A.M() );

     for ( i = 0; i < A.N(); i++ )
         for ( j = 0; j < A.M(); j++ )
             T(i,j) = A.a(j,i);

     return T;
}


double Spur( const Matrix &A )
{    double s = 0.0;
     int i;

     for ( i = 0; i < A.M() && i < A.N(); i++ )
         s += A.a(i,i);

     return s;
}

/* ----------------------------------------------------------------------- */

Matrix operator*( double s, const Matrix &A )
{    int i, j;
     double *array;

     if ( ( array = new double[A.M()*A.N()] ) == NULL )
        Fehler("new erfolglos");

     for ( i = 0; i < A.M(); i++ )
         for ( j = 0; j < A.N(); j++ )
             array[i+j*A.M()] = s * A.a(i,j);

     return Matrix( A.M(), A.N(), array, A.M() );
}


Matrix operator*( const Matrix &A, double s )
{    return s * A;
}


Matrix operator-( const Matrix &A )
{    return -1.0 * A;
}


Matrix operator+( const Matrix &A, const Matrix &B )
{    int i, j;
     double *array;

     if ( A.M() != B.M() || A.N() != B.N() )
        Fehler("Matrizenaddition nicht definiert");

     if ( ( array = new double[A.M()*A.N()] ) == NULL )
        Fehler("new erfolglos");

     for ( i = 0; i < A.M(); i++ )
         for ( j = 0; j < A.N(); j++ )
             array[i+j*A.M()] = A.a(i,j) + B.a(i,j);

     return Matrix( A.M(), A.N(), array, A.M() );
}


Matrix operator-( const Matrix &A, const Matrix &B )
{    int i, j;
     double *array;

     if ( A.M() != B.M() || A.N() != B.N() )
        Fehler("Matrizensubtraktion nicht definiert");

     if ( ( array = new double[A.M()*A.N()] ) == NULL )
        Fehler("new erfolglos");

     for ( i = 0; i < A.M(); i++ )
         for ( j = 0; j < A.N(); j++ )
             array[i+j*A.M()] = A.a(i,j) - B.a(i,j);

     return Matrix( A.M(), A.N(), array, A.M() );
}


Matrix operator*( const Matrix &A, const Matrix &B )
{    int i, j, k;
     double *array;
     double val;

     if ( A.N() != B.M() )
        Fehler("Matrizenmultiplikation nicht definiert");

     if ( ( array = new double[A.M()*B.N()] ) == NULL )
        Fehler("new erfolglos");

     for ( i = 0; i < A.M(); i++ )
         for ( j = 0; j < B.N(); j++ )
         {   val = 0.0;
             for ( k = 0; k < A.N(); k++ )
                 val += A.a(i,k) * B.a(k,j);
             array[i+A.M()*j] = val;
         }

     return Matrix( A.M(), B.N(), array, A.M() );
}


std::string Matrix::str()
{
  std::string retval;

  retval.reserve(10000);
  for(int row = 0; row < M(); ++row)
  {
      retval += "[";
      for(int col = 0; col < N(); ++col)
      {
          char buffer[256];  // make sure this is big enough!!!
          snprintf(buffer, sizeof(buffer), "%0.7e\t", a(row,col));
          retval += buffer;
      }
      retval += " ]\n";
  }
  return retval;
}
