//
//  BLAS_Matrix.h
//  BlueBookFunctions
//
//  Created by Peter Huber on 2018-01-11.
//  Copyright © 2018 Peter Huber. All rights reserved.
//

// NOTE: Only double and double_complex matrices are modeled. This file is intened to make it easier to create and use BLAS and LAPACK matrices and functions. It is written in C to allow easier portability between Macs and Windows machines. There are some preprocessor directives for correct compilation but the underlying projects will need to load the correct framework/library.

// NOTE: All buffers are assumed to be in COLUMN-MAJOR order

// NOTE: We "UPLO = U" symantics for symmetric matrices

#ifndef BLAS_Matrix_h
#define BLAS_Matrix_h

#include <stdio.h>
#include <stdlib.h>

#if defined(__APPLE__)

#include <Accelerate/Accelerate.h>

#endif

// The different matrix types that we can handle. Note that the idea is that eventually, all functions will apply to all matrix types, but this will probably be an evolutionary "if I need it" sort of thing.
typedef enum {generalMatrix, bandedMatrix, lowerTriangularMatrix, upperTriangularMatrix, diagonalMatrix, symmetricMatrix, hermitianMatrix, vectorMatrix} MatrixType;

// This is 2018. We only work with double (and double complex) types
typedef enum {doublePrecisionMatrix, complexPrecisionMatrix} MatrixPrecision;

// Since this is C, we don't have access to classes. So we define a struct instead. Most of the routines in this library take a BLAS_Matrix struct pointer as their first argument. I suppose I could pass the structs as values.
typedef struct _tagMatrix {
    
    MatrixType type;
    MatrixPrecision precision;
    unsigned int numRows;
    unsigned int numCols;
    unsigned int numSuperDiags;
    unsigned int numSubDiags;
    
    size_t bufferSize; // in bytes
    __CLPK_doublereal *buffer;
    
} BLAS_Matrix;

// Error codes that can be returned by some routines
#define GET_VALUE_ERROR DBL_MAX

// Get the value stored at the given row/column
__CLPK_doublereal GetDoubleValue(const BLAS_Matrix *theMatrix, unsigned int row, unsigned int col);
__CLPK_doublecomplex GetComplexValue(const BLAS_Matrix *theMatrix, unsigned int row, unsigned int col);

// Set the value at the given row/column
bool SetDoubleValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, __CLPK_doublereal value);
bool SetComplexValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, __CLPK_doublecomplex value);

// Create a new vector (convenience routine)
BLAS_Matrix *CreateVector(MatrixPrecision precision, unsigned int numElements);
// Create a new matrix
BLAS_Matrix *CreateMatrix(MatrixType type, MatrixPrecision precision,  unsigned int rows, unsigned int columns, unsigned int subDiagonals, unsigned int superDiagonals);
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


// Simple function to get the string representation of the matrix. Note that this should probably only be used for small matrices (ie: during debugging)
char *MatrixAsString(BLAS_Matrix *theMatrix);

#endif /* BLAS_Matrix_h */
