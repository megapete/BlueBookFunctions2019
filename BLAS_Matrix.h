//
//  BLAS_Matrix.h
//  BlueBookFunctions
//
//  Created by Peter Huber on 2018-01-11.
//  Copyright © 2018 Peter Huber. All rights reserved.
//

// NOTE: Only double and double_complex matrices are modeled. This file is intened to make it easier to create and use BLAS and LAPACK matrices and functions. It is written in C to allow easier portability between Macs and Windows machines. There are a bunch of preprocesser directives for correct compilation but the underlying projects will need to load the correct framework/library.

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
    
    __CLPK_doublereal *buffer;
    
} BLAS_Matrix;

// Error codes that can be returned by some routines
#define GET_VALUE_ERROR DBL_MAX

// Get the value stored at the given row/column
__CLPK_doublereal GetDoubleValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col);
__CLPK_doublecomplex GetComplexValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col);

// Set the value at the given row/column
bool SetDoubleValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, __CLPK_doublereal value);
bool SetComplexValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, __CLPK_doublecomplex value);)

// Create a new vector (convenience routine)
BLAS_Matrix *CreateVector(MatrixType type, MatrixPrecision precision, unsigned int numElements);
// Create a new matrix
BLAS_Matrix *CreateMatrix(MatrixType type, MatrixPrecision precision,  unsigned int rows, unsigned int columns, unsigned int subDiagonals, unsigned int superDiagonals);



#endif /* BLAS_Matrix_h */
