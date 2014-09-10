/**
 * \file matrix.h
 * \author Stefan Koerkel and JRG ExpDesign
 * \brief Matrix class for handling matrices.
 *
 *  Matrix class for handling matrices. Basic structure for arithmetic
 *  computations is given.
 */

#ifndef BLOCKSQP_MATRIX_HPP
#define BLOCKSQP_MATRIX_HPP

#include <stdio.h>
#include <string>

namespace blockSQP
{

/* ----------------------------------------------------------------------- */

void Fehler( const char *F );

/* ----------------------------------------------------------------------- */

/// Count constructor calls
extern int Ccount;
/// Count destructor calls
extern int Dcount;
/// Count assign operator calls
extern int Ecount;

/* ----------------------------------------------------------------------- */

/** Matrix class
 */

class Matrix
{  private:
      int malloc( void );                           ///< memory allocation
      int free( void );                             ///< memory free

   public:
      int m;                                        ///< internal number of rows
      int n;                                        ///< internal number of columns
      int ldim;                                     ///< internal leading dimension not necesserily equal to m or n
      double *array;                                ///< array of how the matrix is stored in the memory
      int tflag;                                    ///< 1 if it is a Teilmatrix

      Matrix( int = 1, int = 1, int = -1 );         ///< constructor with standard arguments
      Matrix( int, int, double*, int = -1 );
      Matrix( const Matrix& );
      ~Matrix( void );

      int M( void ) const;                          ///< number of rows
      int N( void ) const;                          ///< number of columns
      int LDIM( void ) const;                       ///< leading dimensions
      double *ARRAY( void );                        ///< returns pointer to data array
      int TFLAG( void ) const;                      ///< returns this->tflag (1 if it is a submatrix and does not own the memory and 0 otherwise)

      Matrix &operator=( const Matrix& );
      virtual double &operator()( int i, int j );
      virtual double &operator()( int i );
      virtual double a( int i, int j ) const;
      virtual double a( int i ) const;

      Matrix &Dimension( int, int = 1, int = -1 );
      Matrix &set_m( int );
      Matrix &Initialisieren( double (*)( int, int ) );
      Matrix &Initialisieren( double );
      Matrix &Streiche_Zeilen( int );
      Matrix &Streiche_Spalten( int );

      /// Returns just a pointer to the full matrix
      Matrix& Teilmatrix( const Matrix&, int, int, int = 0, int = 0 );
      /// Matrix that points on <tt>ARRAY</tt>
      Matrix& Arraymatrix( int M, int N, double* ARRAY, int LDIM = -1 );
      Matrix& Teilmatrix_del( void );

      Matrix &operator+=( const Matrix & );
      Matrix &operator-=( const Matrix & );
      Matrix &operator*=( double );

      int convert2sparse( int*, int*, int*,
                          int**, double**, int**,
                          double = 1.0e-12 ) const;

      /** Read in matrix from console input
       */
      Matrix &Eingabe( void );

      /** Read in of a matrix from data file.
       * The format has to be exactly the same which comes out of Ausgabe().
       */
      Matrix &Einlesen( FILE* = stdin );

      /** Flag == 0: bracket output
        * Flag == 1: Matlab output
        * else: plain output
        */
      const Matrix &Ausgabe( FILE* = stdout,   ///< file for output
                             int = 13,         ///< number of digits
                             int = 1           ///< Flag for format
                           ) const;

      std::string str();                ///< returns std::string representation of the matrix



};

class SymMatrix : public Matrix
{
    protected:
        int malloc( void );
        int free( void );

    public:
        SymMatrix( int = 1 );
        SymMatrix( int, double* );
        SymMatrix( int, int, int );
        SymMatrix( int, int, double*, int = -1 );
        SymMatrix( const Matrix& );
        SymMatrix( const SymMatrix& );
//         ~SymMatrix( void );

        SymMatrix &operator=( const Matrix& );
        virtual double &operator()( int i, int j );
        virtual double &operator()( int i );
        virtual double a( int i, int j ) const;
        virtual double a( int i ) const;

        SymMatrix &Dimension( int = 1 );
        SymMatrix &Dimension( int, int, int );
        SymMatrix &set_m( int );
        SymMatrix &Initialisieren( double (*)( int, int ) );
        SymMatrix &Initialisieren( double );
        SymMatrix &Streiche_Zeilen( int );
        SymMatrix &Streiche_Spalten( int );

        SymMatrix& Teilmatrix( const Matrix&, int, int, int = 0, int = 0 );
        SymMatrix& Arraymatrix( int, double* );
        SymMatrix& Arraymatrix( int, int, double*, int = -1 );
        SymMatrix& Teilmatrix_del( void );

//         Matrix &SymMatrix::convert2Matrix( void );
};


/* ----------------------------------------------------------------------- */

double delta( int, int );
double nullentry( int, int );
double einsentry( int, int );
double randomentry( int, int );

int operator==( const Matrix&, const Matrix& );
int operator!=( const Matrix&, const Matrix& );

/// Overwrites \f$ A \f$ with its transpose \f$ A^T \f$
Matrix Transponierte( const Matrix& A);
/// Computes \f$ T = A^T \f$
Matrix &Transponierte( const Matrix &A, Matrix &T );
/// Computes the trace of \f$ A \f$
double Spur( const Matrix& A);

/// Computes \f$ s\cdot  A \f$ with matrix \f$ A\f$ and scalar \f$ s \f$
Matrix operator*( double s, const Matrix& A);
/// Computes \f$ A \cdot s \f$ with matrix \f$ A\f$ and scalar \f$ s \f$
Matrix operator*( const Matrix& A, double s);
/// Computes \f$ - A \f$
Matrix operator-( const Matrix& A);
/// Computes \f$ A + B \f$ for matrices \f$ A\f$ and \f$ B\f$
Matrix operator+( const Matrix& A, const Matrix& B);
/// Computes \f$ A - B \f$ for matrices \f$ A\f$ and \f$ B\f$
Matrix operator-( const Matrix& A, const Matrix& B);
/** Computes \f$ A * B \f$ for matrices
 * \f$ A \in R^{m\times n }\f$ and \f$ B\in R^{n \times k}\f$
 */
Matrix operator*( const Matrix& A, const Matrix& B);

/* ----------------------------------------------------------------------- */

} // namespace blockSQP

#endif
