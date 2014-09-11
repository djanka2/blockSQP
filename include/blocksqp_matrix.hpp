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

#include "blocksqp_defs.hpp"

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

      virtual double &operator()( int i, int j );
      virtual double &operator()( int i );
      virtual double a( int i, int j ) const;
      virtual double a( int i ) const;

      Matrix &Dimension( int, int = 1, int = -1 );
      Matrix &Initialize( double (*)( int, int ) );
      Matrix &Initialize( double );

      /// Returns just a pointer to the full matrix
      Matrix& Submatrix( const Matrix&, int, int, int = 0, int = 0 );
      /// Matrix that points on <tt>ARRAY</tt>
      Matrix& Arraymatrix( int M, int N, double* ARRAY, int LDIM = -1 );

      //int convert2sparse( int*, int*, int*,
                          //int**, double**, int**,
                          //double = 1.0e-12 ) const;

      /** Flag == 0: bracket output
        * Flag == 1: Matlab output
        * else: plain output
        */
      const Matrix &Print( FILE* = stdout,   ///< file for output
                             int = 13,         ///< number of digits
                             int = 1           ///< Flag for format
                           ) const;
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
        ~SymMatrix( void );

        virtual double &operator()( int i, int j );
        virtual double &operator()( int i );
        virtual double a( int i, int j ) const;
        virtual double a( int i ) const;

        SymMatrix &Dimension( int = 1 );
        SymMatrix &Dimension( int, int, int );
        SymMatrix &Initialize( double (*)( int, int ) );
        SymMatrix &Initialize( double );

        SymMatrix& Submatrix( const Matrix&, int, int, int = 0, int = 0 );
        SymMatrix& Arraymatrix( int, double* );
        SymMatrix& Arraymatrix( int, int, double*, int = -1 );

//         Matrix &SymMatrix::convert2Matrix( void );
};

double delta( int, int );
/// Overwrites \f$ A \f$ with its transpose \f$ A^T \f$
Matrix Transpose( const Matrix& A);
/// Computes \f$ T = A^T \f$
Matrix &Transpose( const Matrix &A, Matrix &T );

} // namespace blockSQP

#endif
