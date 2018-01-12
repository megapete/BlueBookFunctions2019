//
//  BLAS_Matrix.h
//  BlueBookFunctions
//
//  Created by Peter Huber on 2018-01-11.
//  Copyright © 2018 Peter Huber. All rights reserved.
//

// Only double and double_complex matrices are modeled. This file is intened to make it easier to create and use BLAS and LAPACK matrices and functions. It is written in C to allow easier portability between Macs and Windows machines. There are a bunch of preprocesser directives for correct compilation but the underlying projects will need to load the correct framework/library.

#ifndef BLAS_Matrix_h
#define BLAS_Matrix_h

#include <stdio.h>

#if defined(__APPLE__)

#include <Accelerate/Accelerate.h>

#endif

typedef enum {generalMatrix, bandedMatrix, triangularMatrix, diagonalMatrix, symmetricMatrix, hermitianMatrix, vectorMatrix} MatrixType;
typedef enum {doublePrecisionMatrix, complexPrecisionMatrix} MatrixPrecision;

__CLPK_doublereal *AllocateMatrix(MatrixType type, MatrixPrecision precision, int rows, int columns);

#endif /* BLAS_Matrix_h */
