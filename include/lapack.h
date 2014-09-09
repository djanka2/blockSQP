/***************************************************************************
   VPLAN Source Code
   Version: $Id: lapack.h 992 2012-05-14 14:36:48Z skoerkel $
   Author: Stefan Koerkel and JRG ExpDesign
 ***************************************************************************/
/**
 * \file lapack.h
 * \author Stefan Koerkel
 * \brief C interface sto fortran lapack routines
 */

#ifndef LAPACK_H
#define LAPACK_H

#include "fornorm.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* ----------------------------------------------------------------------- */

#define dsyevx __FORNORML(dsyevx, DSYEVX)
#define dsyev __FORNORML(dsyev, DSYEV)
#define dsytrf __FORNORML(dsytrf, DSYTRF)
#define dsytri __FORNORML(dsytri, DSYTRI)
#define dspev __FORNORML(dspev, DSPEV)
#define dsptrf __FORNORML(dsptrf, DSPTRF)
#define dsptri __FORNORML(dsptri, DSPTRI)
//#define dpotrf __FORNORML(dpotrf, DPOTRF)
#define dpotri __FORNORML(dpotri, DPOTRI)
#define dlamch __FORNORML(dlamch, DLAMCH)
#define dgeqp3 __FORNORML(dgeqp3, DGEQP3)
#define dgeqrf __FORNORML(dgeqrf, DGEQRF)
#define dormqr __FORNORML(dormqr, DORMQR)
#define dgetrf __FORNORML(dgetrf, DGETRF)
#define dgetri __FORNORML(dgetri, DGETRI)
#define dcopy __FORNORMB(dcopy, DCOPY)
//#define dgemm __FORNORMB(dgemm, DGEMM)
#define dsymm __FORNORMB(dsymm, DSYMM)
#define dgesvd __FORNORMB(dgesvd, DGESVD)

/* ----------------------------------------------------------------------- */

/// see LAPACK documentation
void dsyevx( char *jobz, char *range, char *uplo,
             int *n, double *a, int *lda,
             double *vl, double *vu, int *il, int *iu, double *abstol,
             int *m, double *w, double *z, int *ldz,
             double *work, int *lwork, int *iwork, int *ifail, int *info,
             int strlen_jobz, int strlen_range, int strlen_uplo );

void dsyev( char *jobz, char *uplo, int *n, double *a, int *lda,
            double *w, double *work, int *lwork, int *info,
            int strlen_jobz, int strlen_uplo );

void dsytrf( char *uplo, int *n, double *a, int *lda, int *ipiv,
             double *work, int *lwork, int *info, int strlen_uplo );

void dsytri( char *uplo, int *n, double *a, int *lda, int *ipiv,
             double *work, int *info, int strlen_uplo );

void dspev( char *jobz, char *uplo, int *n, double *ap, double *w, double *z, int *ldz,
            double *work, int *info, int strlen_jobz, int strlen_uplo );

void dsptrf( char *uplo, int *n, double *ap, int *ipiv,
             int *info, int strlen_uplo );

void dsptri( char *uplo, int *n, double *ap, int *ipiv,
             double *work, int *info, int strlen_uplo );

void dpotrf( char *uplo, int *n, double *a, int *lda, int *info,
             int strlen_uplo );

void dpotri( char *uplo, int *n, double *a, int *lda, int *info,
             int strlen_uplo );

double dlamch( char *cmach, int strlen_cmach );

void dgeqp3( int *m, int *n, double *a, int *lda, int *jpvt, double *tau,
             double *work, int *lwork, int *info );

void dgeqrf( int *m, int *n, double *a, int *lda, double *tau,
             double *work, int *lwork, int *info );

void dormqr( char *side, char *trans,
             int *m, int *n, int *k, double *a, int *lda,
             double *tau, double *c, int *ldc,
             double *work, int *lwork, int *info,
             int strlen_side, int strlen_trans );

void dgetrf( int *m, int *n, double *a, int *lda, int *ipiv, int *info );

void dgetri( int *n, double *a, int *lda,
             int *ipiv, double *work, int *lwork, int *info );

void dcopy( const int *n, const double *x, const int *incx,
            double *y, const int *incy );

void dgesvd ( char *jobu, char *jobvt,
              int *m, int *n, double *A, int *lda,
              double *S,
              double *U, int *ldu,
              double *VT, int *ldvt,
              double *work, int *lwork, int *info);

//void dgemm( char *transa, char *transb, int *m, int *n, int *k,
            //double *alpha, double *a, int *lda, double *b, int *ldb,
            //double *beta, double *c, int *ldc,
            //int strlen_transa, int strlen_transb );

void dsymm( char *side, char *uplo, int *m, int *n, double *alpha,
            double *a, int *lda,
            double *b, int *ldb, double *beta,
            double *c, int *ldc,
            int strlen_side, int strlen_uplo );

/* ----------------------------------------------------------------------- */

#ifdef __cplusplus
}
#endif

#endif
