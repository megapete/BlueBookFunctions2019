//
//  BLAS_Matrix.h
//  BlueBookFunctions
//
//  Created by Peter Huber on 2018-01-11.
//  Copyright © 2018 Peter Huber. All rights reserved.
//

// NOTE: Only double and double_complex matrices are modeled. This file is intened to make it easier to create and use BLAS and LAPACK matrices and functions. It is written in C to allow easier portability between Macs and Windows machines. There are some preprocessor directives for correct compilation but the underlying projects will need to load the correct framework/library.

// NOTE: All buffers are assumed to be in COLUMN-MAJOR order.

// NOTE: We use "UPLO = U" symantics for symmetric, positive-definite, and hermitian matrices

#ifndef BLAS_Matrix_h
#define BLAS_Matrix_h

#include <stdio.h>
#include <stdlib.h>

#if defined(__APPLE__)

#include <Accelerate/Accelerate.h>

#endif

// The different matrix types that we can handle. Note that the idea is that eventually, all functions will apply to all matrix types, but this will probably be an evolutionary "if I need it" sort of thing.
typedef enum {generalMatrix, bandedMatrix, lowerTriangularMatrix, upperTriangularMatrix, diagonalMatrix, symmetricMatrix, hermitianMatrix, positiveDefiniteMatrix, vectorMatrix} MatrixType;

// This is 2018. We only work with double (and double complex) types
typedef enum {doublePrecisionMatrix, complexPrecisionMatrix} MatrixPrecision;

// Error codes that can be returned by the Matrix creation routines
typedef enum {NO_MATRIX_ERROR, ZERO_DIMENSION_MATRIX_ERROR, VECTOR_DIMENSION_MATRIX_ERROR, EXPECTED_SQUARE_MATRIX_ERROR, NONCOMPLEX_HERMITIAN_MATRIX_ERROR, BUFFER_ALLOCATION_MATRIX_ERROR, ILLEGAL_PRECISION_MATRIX_ERROR, BUFFERSIZE_MATRIX_ERROR, ILLEGAL_TYPE_MATRIX_ERROR, PRECISION_MISMATCH_MATRIX_ERROR, DIMENSION_MISMATCH_MATRIX_ERROR, LAPACK_ROUTINE_ERROR, SOLVE_MATRIX_ERROR, ILLEGAL_INDEX_MATRIX_ERROR, UNREACHABLE_INDEX_MATRIX_ERROR, ILLEGAL_FACTORING_SOLVE_ERROR} MatrixStatus;

// Since this is C, we don't have access to classes. So we define a struct instead. Most of the routines in this library take a BLAS_Matrix struct pointer as their first argument (I suppose I could pass the structs as values). Routines that call the matrix creation routines should check the value of status before using the Matrix. In the event of an LAPACK call error where the "info" argument holds the error code, the matrix routine will return LAPACK_ROUTINE_ERROR with the code in lapack_info field. If an error occurs during one of the Solve...LinearSystem calls, the status will hold SOLVE_MATRIX_ERROR with the underlying LAPACK info error in lapack_info.
typedef struct _tagMatrix {
    
    MatrixType type;
    MatrixPrecision precision;
    unsigned int numRows;
    unsigned int numCols;
    unsigned int numSuperDiags;
    unsigned int numSubDiags;
    
    MatrixStatus status;
    int lapack_info;
    
    size_t bufferSize; // in bytes
    __CLPK_doublereal *buffer;
    
} BLAS_Matrix;


// The general solve routine can return all sorts of information. We create a struct to hold it all. There are more simple routines (which ultimately call the most general) that simply return results. See the LAPACK documentation for the meaning of most of the members of the struct. Note that all of the pointers in this class are allocated on the heap. Call the ReleaseSolveResult to free all the memory that has been allocated. Any routine that returns a SolveResult will set the resultIsValid field before returning. If resultIsValid is false, then it is guaranteed that all the pointers in the struct have been released and set to NULL.
typedef struct _tagSolveResult {
    
    __CLPK_integer lda; // always N for an NxN coefficient matrix
    void *A;
    __CLPK_integer ldaf; // always the same as lda
    void *AF;
    __CLPK_integer *ipiv;
    char equed;
    void *R;
    void *C;
    
    void *B;
    __CLPK_integer ldb;
    
    __CLPK_doublereal rCond;
    __CLPK_doublereal *fErr;
    __CLPK_doublereal *bErr;
    
    __CLPK_integer info;
    
    void *X;
    __CLPK_integer ldx;
    
    MatrixStatus status;
    
} SolveResult;

// Enums for the input flags for the Solve routines
typedef enum {BLAS_NoFactor, BLAS_Factored/* UNIMPLEMENTED */, BLAS_Equilibrate} BLAS_Factoring;
typedef enum {BLAS_NoTranspose, BLAS_Transpose, BLAS_ConjugateTranspose} BLAS_Transposing;
typedef enum {BLAS_NoEquilibration, BLAS_RowEquilibration, BLAS_ColumnEquilibration, BLAS_BothEquilibration} BLAS_Equilibrating;

// Error codes that can be returned by the Get...Value routines
#define GET_VALUE_ERROR DBL_MAX

// Get the value stored at the given row/column. In the event of an error, the functions return GET_VALUE_ERROR (for the Complex type, this is returned in the real field) and status is set to the error code. If there is no error, status will hold NO_MATRIX_ERROR. Note that the status argument can safely be set to NULL if the caller is not interested in error codes.
__CLPK_doublereal GetDoubleValue(const BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, MatrixStatus *status);
__CLPK_doublecomplex GetComplexValue(const BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, MatrixStatus *status);

// Set the value at the given row/column. If an error occurs, the error code is returned (NO_MATRIX_ERROR if no error occurs).
MatrixStatus SetDoubleValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, __CLPK_doublereal value);
MatrixStatus SetComplexValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, __CLPK_doublecomplex value);

// Create a new vector (convenience routine)
BLAS_Matrix *CreateVector(MatrixPrecision precision, unsigned int numElements);
// Create a new matrix
BLAS_Matrix *CreateMatrix(MatrixType type, MatrixPrecision precision,  unsigned int rows, unsigned int columns, unsigned int subDiagonals, unsigned int superDiagonals);

// Delete a matrix using this function to avoid memory leaks. Note that this function should only be used if the matrix was created using either CreateVector() or CreateMatrix() or if the matrix struct itself was created using malloc(). That is, do not call this function if you create the matrix on the stack (in which case the calling rouitne is responsible for freeing the buffer's memory)
void ReleaseMatrix(BLAS_Matrix *matrix);

// Copy an existing matrix (or vector)
BLAS_Matrix *CopyMatrix(const BLAS_Matrix *srcMatrix);

// Get the transpose of a matrix
BLAS_Matrix *TransposeMatrix(const BLAS_Matrix *srcMatrix);

// Matrix multiplication. The function is basically a wrapper over the BLAS routine, which is defined as follows:
// D := alpha * op(A) * op(B) + beta * C (note that 'D' does not exist in the BLAS, with the result clobbering C instead)
// op(A) and op(B) depend on the transA and transB arguments passed in
// Note that as far as I can tell, only general matrices are treated by the BLAS
// A and B cannot be NULL
// C can be NULL, but in the event that it is not, it will not be clobbered by the answer (as is the case in BLAS)
BLAS_Matrix *MultiplyDoubleMatrices(__CLPK_doublereal alpha, int transA, BLAS_Matrix *A, int transB, BLAS_Matrix *B, __CLPK_doublereal beta, BLAS_Matrix *C);

BLAS_Matrix *MultiplyComplexMatrices(__CLPK_doublecomplex alpha, int transA, BLAS_Matrix *A, int transB, BLAS_Matrix *B, __CLPK_doublecomplex beta, BLAS_Matrix *C);

// Get the inverse of a matrix, if it is invertible (otherwise, return NULL)
BLAS_Matrix *Inverse(const BLAS_Matrix *srcMatrix);

// The most general solve routine available. Calling routines should probably use one of the more specific solver routines that this library makes available.
// Notes for the call: B _MUST_ be a "General Matrix" type. A can be anything.
SolveResult Solve(bool useExtendedMethod, BLAS_Matrix *A, BLAS_Matrix *B, BLAS_Factoring factor, BLAS_Transposing transpose, BLAS_Equilibrating equilibrate);

// Call this to safely release all allocated memory used in creating the solve result
void ReleaseSolveResult(SolveResult sResult);

// Given matrices A and B, find X for the system AX = B. Uses the simplest of the available solve routines. The calling routine should check the status field of the returned matrix (if it is not equal to NO_MATRIX_ERROR, then the matrix is not valid),
BLAS_Matrix *SolveSimpleLinearSystem(BLAS_Matrix *A, BLAS_Matrix *B);

// Simple function to get the string representation of the matrix. Note that this should probably only be used for small matrices (ie: during debugging)
char *MatrixAsString(BLAS_Matrix *theMatrix);

#endif /* BLAS_Matrix_h */
