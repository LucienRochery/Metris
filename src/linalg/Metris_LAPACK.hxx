// lapack_declarations.h
#ifndef METRIS_LAPACK
#define METRIS_LAPACK

// Unfortunately, some LAPACKE routines call malloc, and we use them for very
// small matrices (2x2, 3x3...)
#ifdef METRIS_USE_LAPACK

extern "C" {

// LU factorization for a general m×n matrix
void dgetrf_(int* m, int* n, double* a, int* lda,
             int* ipiv, int* info);

// Compute eigenvalues and (optionally) eigenvectors of a real symmetric matrix
void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda,
            double* w, double* work, int* lwork, int* info);

// Compute inverse of a symmetric positive-definite packed matrix using its Cholesky factor
void dpptri_(char* uplo, int* n, double* ap, int* info);

// Cholesky factorization of symmetric positive-definite packed matrix
void dpptrf_(char* uplo, int* n, double* ap, int* info);

// Bunch–Kaufman factorization of symmetric packed matrix
void dsptrf_(char* uplo, int* n, double* ap, int* ipiv, int* info);

// Compute inverse using previously computed factorization from dsptrf
void dsptri_(char* uplo, int* n, double* ap, int* ipiv, double* work, int* info);

// Compute matrix inverse from its LU factorization (after dgetrf)
void dgetri_(int* n, double* a, int* lda, int* ipiv,
             double* work, int* lwork, int* info);

} // extern "C"

#endif

#endif