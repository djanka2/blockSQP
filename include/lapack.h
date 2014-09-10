#ifndef LAPACK_H
#define LAPACK_H

#ifdef __cplusplus
extern "C"
{
#endif

void dsyev_( char *jobz, char *uplo, int *n, double *a, int *lda,
            double *w, double *work, int *lwork, int *info,
            int strlen_jobz, int strlen_uplo );

void dspev_( char *jobz, char *uplo, int *n, double *ap, double *w, double *z, int *ldz,
            double *work, int *info, int strlen_jobz, int strlen_uplo );

void dgetrf_( int *m, int *n, double *a, int *lda, int *ipiv, int *info );

void dgetri_( int *n, double *a, int *lda,
             int *ipiv, double *work, int *lwork, int *info );

#ifdef __cplusplus
}
#endif

#endif
